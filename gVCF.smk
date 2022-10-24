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
        expand("{gvcfs}/reblock/{chr}/{sample}.{chr}.g.vcf.gz", gvcfs=config['gVCF'], chr = main_chrs, sample = sample_names),
        rules.Aligner_all.input
    default_target: True

def get_mem_mb_CalibrateDragstrModel(wildcrads, attempt):
    return attempt * int(config['CalibrateDragstrModel']['mem'])

# calibrate model
# step for HaplotypeCaller in dragen mode
# for better resolution in complex regions
rule CalibrateDragstrModel:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output:
        dragstr_model = config['BAM'] + "/{sample}-dragstr.txt"
    priority: 26
    params:
        str_ref = config['RES'] + config['str_ref']
    resources: mem_mb = get_mem_mb_CalibrateDragstrModel
    conda: "envs/preprocess.yaml"
    log: config['LOG'] + '/' + "{sample}_calibratedragstr.log"
    benchmark: config['BENCH'] + "/{sample}_calibrate_dragstr.txt"
    shell:
        "{gatk} CalibrateDragstrModel -R {ref} -I {input.bam} -O {output} -str {params.str_ref} 2>{log}"

# verifybamid
# verifybamid has some bugs and require samtools lower version
# to avoid any bugs verifybamid step runs in different conda enviroment
#
rule verifybamid:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output:
        VBID_stat = config['STAT'] + '/contam/{sample}_verifybamid.pca2.selfSM'
    # end of this file hardcoded in Haplotypecaller and read_contam_w
    threads: config['verifybamid']['n']
    benchmark: config['BENCH'] + "/{sample}_verifybamid.txt"
    priority: 27
    params:
        VBID_prefix = config['STAT'] + '/contam/{sample}_verifybamid.pca2',
        SVD = config['RES'] + config['verifybamid_exome']
    conda: 'envs/verifybamid.yaml'
    shell:
        """
        verifybamid2 --BamFile {input.bam} --SVDPrefix {params.SVD} --Reference {ref} --DisableSanityCheck --NumThread {threads} --Output {params.VBID_prefix}
        """

def get_chrom_capture_kit(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    chr = wildcards.chr
    if SAMPLEINFO[wildcards['sample']]['sample_type'].startswith('illumina_wgs'):
        capture_kit_chr_path = chr
    else:
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

def get_mem_mb_HaplotypeCaller(wildcrads, attempt):
    return attempt * int(config['HaplotypeCaller']['mem'])


#find SNPs from bams
# change RAM too function
rule HaplotypeCaller:
    input:
        # check what bam file we need to use (with or without additional cleanup)
        bams = rules.markdup.output.mdbams,
        model = rules.CalibrateDragstrModel.output.dragstr_model,
        contam= rules.verifybamid.output.VBID_stat,
        bai = rules.markdup_index.output.mdbams_bai
    output:
        gvcf= config['gVCF'] + "/{chr}/{sample}.{chr}.g.vcf.gz"
    log:
        HaplotypeCaller=config['LOG'] + "/{sample}_{chr}_haplotypecaller.log"
    benchmark:
        config['BENCH'] + "/{sample}_{chr}_haplotypecaller.txt"
    conda: "envs/preprocess.yaml"
    resources: mem_mb = get_mem_mb_HaplotypeCaller
    threads: config['HaplotypeCaller']['n']
    params:
        dbsnp = config['RES'] + config['dbsnp'],
        padding=100,  # extend intervals to this bp
        contam_frac = read_contam_w, # get contamination fraction per sample
        # command to get path to capture_kit interval list from SAMPLEFILE
        interval= get_chrom_capture_kit,
    priority: 28
    shell:
        "{gatk} HaplotypeCaller \
                 -R {ref} -L {params.interval} -ip {params.padding} -D {params.dbsnp} -ERC GVCF --contamination {params.contam_frac} \
                 -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
                 -I {input.bams} -O {output.gvcf}  --native-pair-hmm-threads {threads} \
                  --dragen-mode true --dragstr-params-path {input.model} 2> {log.HaplotypeCaller}"

rule reblock_gvcf:
    input:
        gvcf = rules.HaplotypeCaller.output.gvcf,
    output: gvcf_reblock = config['gVCF'] + "/reblock/{chr}/{sample}.{chr}.g.vcf.gz"
    log: Reblock=config['LOG'] + "/{sample}_{chr}_reblock.log"
    benchmark:
        config['BENCH'] + "/{sample}_{chr}_reblock.txt"
    conda: "envs/preprocess.yaml"
    priority: 29
    params:
        dbsnp=config['RES'] + config['dbsnp'],
    shell:
        "{gatk} ReblockGVCF --keep-all-alts -D {params.dbsnp} -R {ref} -V {input.gvcf} -O {output.gvcf_reblock} -G StandardAnnotation -G AS_StandardAnnotation 2> {log}"
