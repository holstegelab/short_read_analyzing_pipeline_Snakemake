configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")
gatk = config['gatk']
samtools = config['samtools']
bcftools = config['bcftools']
dragmap = config['dragmap']
verifybamid2 = config['verifybamid2']

ref = config['RES'] + config['ref']


wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+"


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO = load_samplefiles('.', config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

rule gVCF_all:
    input:
        expand("{gvcfs}/{chr}/{sample}.{chr}.g.vcf.gz", gvcfs=config['gVCF'], chr = main_chrs, sample = sample_names),
        rules.Aligner_all.input
    default_target: True



# calibrate model
# step for HaplotypeCaller in dragen mode
# for better resolution in complex regions
rule CalibrateDragstrModel:
    input:
        bam = rules.markdup.output.mdbams
    output:
        dragstr_model = config['BAM'] + "/{sample}-dragstr.txt"
    priority: 16
    params:
        str_ref = config['RES'] + config['str_ref']
    conda: "envs/preprocess.yaml"
    log: config['LOG'] + '/' + "{sample}_calibratedragstr.log"
    benchmark: config['BENCH'] + "/{sample}_calibrate_dragstr.txt"
    shell:
        "{gatk} CalibrateDragstrModel -R {ref} -I {input} -O {output} -str {params.str_ref} 2>{log}"

# verifybamid
# verifybamid has some bugs and require samtools lower version
# to avoid any bugs verifybamid step runs in different conda enviroment
#
rule verifybamid:
    input:
        rules.markdup.output.mdbams
    output:
        VBID_stat = config['STAT'] + '/contam/{sample}_verifybamid.pca2.selfSM'
    # end of this file hardcoded in Haplotypecaller and read_contam_w
    threads: config['verifybamid']['n']
    priority: 35
    params:
        VBID_prefix = config['STAT'] + '/contam/{sample}_verifybamid.pca2',
        SVD = config['RES'] + config['verifybamid_exome']
    conda: 'envs/verifybamid.yaml'
    shell:
        """
        verifybamid2 --BamFile {input} --SVDPrefix {params.SVD} --Reference {ref} --DisableSanityCheck --NumThread {threads} --Output {params.VBID_prefix}
        """

def get_chrom_capture_kit(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    chr = wildcards.chr
    capture_kit_chr_path = config['RES'] + config['kit_folder'] + capture_kit + '_hg38/' + capture_kit + '_hg38_' + chr + '.interval_list'
    return capture_kit_chr_path
    # return os.path.join(RESOURCES, 'capture_kits', capture_kit, capture_kit + '.chr' + wildcards.chrom + '.bed')


def read_contam_w(wildcards):
    # filename = rules.verifybamid.output.VBID_stat
    filename = os.path.join(config['STAT'], 'contam', wildcards['sample'] + '_verifybamid.pca2.selfSM')
    with open(filename,'r', encoding='utf-8') as f:
        c = csv.reader(f, delimiter='\t')
        c = iter(c)
        headers = c.__next__()
        data = c.__next__()
        freemix = data[6]
    return freemix

#find SNPs from bams
rule HaplotypeCaller:
    input:
        # check what bam file we need to use (with or without additional cleanup)
        bams = rules.markdup.output.mdbams,
        model = rules.CalibrateDragstrModel.output.dragstr_model,
        contam= rules.verifybamid.output.VBID_stat,
        interval= get_chrom_capture_kit,
        # command to get path to capture_kit interval list from SAMPLEFILE
    output:
        gvcf= config['gVCF'] + "/{chr}/{sample}.{chr}.g.vcf.gz"
    log:
        HaplotypeCaller=config['LOG'] + "/{sample}_{chr}_haplotypecaller.log"
    benchmark:
        config['BENCH'] + "/{sample}_{chr}_haplotypecaller.txt"
    conda: "envs/preprocess.yaml"
    params:
        dbsnp = config['RES'] + config['dbsnp'],
        padding=100,  # extend intervals to this bp
        contam_frac = read_contam_w
    priority: 25
    shell:
        "{gatk} HaplotypeCaller \
                 -R {ref} -L {input.interval} -ip {params.padding} -D {params.dbsnp} -ERC GVCF --contamination {params.contam_frac} \
                 -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
                 -I {input.bams} -O {output.gvcf}  \
                  --dragen-mode true --dragstr-params-path {input.model} 2> {log.HaplotypeCaller}"