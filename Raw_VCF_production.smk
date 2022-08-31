configfile: "Snakefile.cluster.json"
configfile: "Snakefile.paths.yaml"
gatk = config['miniconda'] + config['gatk']
verifybamid2 = config['miniconda'] + config['verifybamid2']
ref = config['RES'] + config['ref']

wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+"


from read_samples import *
from common import *

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

rule all:
    input:
        expand("{vcf}/ALL_chrs.vcf.gz",vcf=config['VCF']),

module Aligner:
    snakefile: 'Aligner.smk'


# check amount of supplementary reads
# if value higher than 0.5% - run additional cleanup steps
# from checkpoint step
# trigger this additional steps only in case if these steps necessary
def check_supp(wildcards):
    with checkpoints.bamstats_all.get(sample=wildcards.sample).output[0].open() as f:
        lines = f.readlines()
        # 4th column (3rd if 0-based) is column with supplementary fraction
        if float((lines[1].split()[3])) >= float(0.005):
            # return cleaned bam in case if clean up necessary
            return rules.sort_back.output.ready_bams
            # return os.path.join(config['BAM'] + '/' + str(wildcards) + '.DeClipped.bam')
        else:
            # return original MD bam if clean up not necessary
            return rules.markdup.output.mdbams
            # return os.path.join(config['BAM'] + '/' + str(wildcards) + '.markdup.bam')

# calibrate model
# step for HaplotypeCaller in dragen mode
# for better resolution in complex regions
rule CalibrateDragstrModel:
    input:
        bam = check_supp
    output:
        dragstr_model = config['BAM'] + "/{sample}-dragstr.txt"
    priority: 16
    params:
        str_ref = config['RES'] + config['str_ref']
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
        check_supp
    output:
        VBID_stat = config['STAT'] + '/contam/{sample}_verifybamid.pca2.selfSM'
    # end of this file hardcoded in Haplotypecaller and read_contam_w
    threads: config['verifybamid']['n']
    priority: 35
    params:
        VBID_prefix = config['STAT'] + '/contam/{sample}_verifybamid.pca2',
        SVD = config['RES'] + config['verifybamid_exome']
    conda: config['CONDA_VERIFYBAMID']
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
    with open(filename,'r') as f:
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
        bams = check_supp,
        model = rules.CalibrateDragstrModel.output.dragstr_model,
        contam= rules.verifybamid.output.VBID_stat,
        interval= get_chrom_capture_kit,
        # command to get path to capture_kit interval list from SAMPLEFILE
    output:
        gvcf="gvcfs/{chr}/{sample}.{chr}.g.vcf.gz"
    log:
        HaplotypeCaller=config['LOG'] + "/{sample}_{chr}_haplotypecaller.log"
    benchmark:
        config['BENCH'] + "/{sample}_{chr}_haplotypecaller.txt"
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

#Genomics DBImport instead CombineGVCFs
rule GenomicDBImport:
    input:
        expand("gvcfs/{chr}/{sample}.{chr}.g.vcf.gz", sample = sample_names, chr = main_chrs)
    log: config['LOG'] + '/' + "GenomicDBImport.{chr}.log"
    benchmark: config['BENCH'] + "/GenomicDBImport.{chr}.txt"
    output:
        dbi=directory("genomicsdb_{chr}"),
        gvcf_list = temp("{chr}_gvcfs.list")
    threads: config['GenomicDBImport']['n']
    # params:
        # N_intervals=5,
        # threads=16,
        # padding = 100
    priority: 30
    shell:
        "ls gvcfs/{wildcards.chr}/*.g.vcf.gz > {output.gvcf_list} && {gatk} GenomicsDBImport --reader-threads {threads} \
        -V {wildcards.chr}_gvcfs.list --intervals {wildcards.chr}  -R {ref} --genomicsdb-workspace-path {output.dbi} \
         --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader 2> {log}"
        # "ls gvcfs/*.g.vcf > gvcfs.list && {gatk} GenomicsDBImport -V gvcfs.list --intervals {chrs}  -R {ref} --genomicsdb-workspace-path {output} \
        #      --max-num-intervals-to-import-in-parallel {params.N_intervals} --reader-threads {params.threads}"

# genotype
rule GenotypeDBI:
    input:
        rules.GenomicDBImport.output.dbi
    output:
        raw_vcfDBI=config['VCF'] + "/Merged_raw_DBI_{chr}.vcf.gz"
    log: config['LOG'] + '/' + "GenotypeDBI.{chr}.log"
    benchmark: config['BENCH'] + "/GenotypeDBI.{chr}.txt"
    params:
        dbsnp = config['RES'] + config['dbsnp']
    priority: 40
    shell:
            "{gatk} GenotypeGVCFs -R {ref} -V gendb://{input} -O {output} -D {params.dbsnp} --intervals {wildcards.chr} 2> {log}"


rule Mergechrs:
    input:
        expand(config['VCF'] + "/Merged_raw_DBI_{chr}.vcf.gz", chr = chr)
    params:
        vcfs = list(map("-I {}/Merged_raw_DBI_{}.vcf.gz".format, config['VCF'], chr))
    log: config['LOG'] + '/' + "Mergechrs.log"
    benchmark: config['BENCH'] + "/Mergechrs.txt"
    output:
        vcf = config['VCF'] + "/ALL_chrs.vcf.gz"
    priority: 45
    shell:
        "{gatk} GatherVcfs {params.vcfs} -O {output} -R {ref} 2> {log} && {gatk} IndexFeatureFile -I {output} "

