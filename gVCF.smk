import pandas as pd
import read_stats
import os
import getpass

configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")
gatk = config['gatk']
samtools = config['samtools']
bcftools = config['bcftools']
dragmap = config['dragmap']
verifybamid2 = config['verifybamid2']

ref = config['RES'] + config['ref']

tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())
tmpdir_alternative = os.path.join(config['tmpdir'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

wildcard_constraints:
    sample="[\w\d_\-@]+",
    # readgroup="[\w\d_\-@]+"


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

module Stat:
    snakefile: 'Stat.smk'
    config: config
use rule verifybamid from Stat

mode = config.get("computing_mode", "WES")
cur_dir = os.getcwd()

def generate_gvcf(wildcards):
    """Generate gvcf file name."""
    res = []
    for sample in sample_names:
        sinfo = SAMPLEINFO[sample]
        if sinfo['sex'] == 'F':
            for chrom in main_chrs_ploidy_female:
                res.append(os.path.join(cur_dir, config['gVCF'], chrom, sample + '.' + chrom + '.female.g.vcf.gz'))
        else:
            for chrom in main_chrs_ploidy_female:
                res.append(os.path.join(cur_dir, config['gVCF'], chrom, sample + '.' + chrom + '.male.g.vcf.gz'))
    return res

rule gVCF_all:
    input:
        generate_gvcf,
        rules.Aligner_all.input
    default_target: True

def get_mem_mb_CalibrateDragstrModel(wildcards, attempt):
    return attempt * int(7100)


rule CalibrateDragstrModel:
    """CalibrateDragstrModel. Estimates the parameters of the dragstr model from a set of aligned reads.
    Provides better STR variant calling.
    """
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output:
        dragstr_model = config['BAM'] + "/{sample}-{sex}-dragstr.txt"
    priority: 26
    params:
        str_ref = os.path.join(config['RES'], config['str_ref']),
        java_options=config['DEFAULT_JAVA_OPTIONS']
    resources: 
        n = 2,
        mem_mb = get_mem_mb_CalibrateDragstrModel
    conda: "envs/preprocess.yaml"
    log: config['LOG'] + '/' + "{sample}_{sex}_calibratedragstr.log"
    benchmark: config['BENCH'] + "/{sample}_{sex}_calibrate_dragstr.txt"
    shell:
        """{gatk} CalibrateDragstrModel --java-options "-Xmx{resources.mem_mb}m {params.java_options}"  -R {ref} -I {input.bam} -O {output} -str {params.str_ref} 2>{log}"""

# verifybamid
# verifybamid has some bugs and require samtools lower version
# to avoid any bugs verifybamid step runs in different conda enviroment



def get_chrom_capture_kit(wildcards):
    """Get capture kit in case of WES, or otherwise chromosome in case of WGS."""
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    chrom = wildcards.chrom
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        capture_kit_chr_path = chrom
    else:        
        capture_kit_chr_path = os.path.join(config['RES'], config['kit_folder'], capture_kit + '_hg38',  capture_kit + '_hg38_' + chrom + '.interval_list')

    return capture_kit_chr_path



def read_contam_w(wildcards):
    """Read contamination from verifybamid output file."""
    filename = os.path.join(config['STAT'], 'contam', wildcards['sample'] + '.' + wildcards['sex'] + '.verifybamid.pca2.selfSM')
    with open(filename,'r', encoding='utf-8') as f:
        c = csv.reader(f, delimiter='\t')
        c = iter(c)
        headers = c.__next__()
        data = c.__next__()
        freemix = data[6]
    return freemix

def get_mem_mb_HaplotypeCaller(wildcrads, attempt):
    """Get memory for HaplotypeCaller."""
    return attempt * int(2500)


def get_chrom_merged_capture_kit(wildcards):
    chrom = wildcards.chrom
    chrom = chrom.replace('XH', 'X').replace('YH', 'Y')
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        capture_kit_chr_path = chrom
    else:
        capture_kit_chr_path = os.path.join(config['RES'], config['kit_folder'], config['MERGED_CAPTURE_KIT'] + '_hg38', config['MERGED_CAPTURE_KIT'] + '_hg38_' + chrom + '.interval_list')
    
    return capture_kit_chr_path

def get_chrom_ploidy(wildcards):
    if wildcards['chrom'] in ['chrXH', 'chrYH']:
        ploidy = 1
    else:
        ploidy = 2
    return ploidy

rule HaplotypeCaller:
    """HaplotypeCaller. Call SNPs and indels for each sample."""
    input:
        # check what bam file we need to use (with or without additional cleanup)
        bams = rules.markdup.output.mdbams,
        model = rules.CalibrateDragstrModel.output.dragstr_model,
        contam= rules.verifybamid.output.VBID_stat,
        bai = rules.markdup.output.mdbams_bai
    output:
        gvcf= ensure(os.path.join(cur_dir, config['gVCF'], "{chrom}/{sample}.{chrom}.{sex}.g.vcf.gz"), non_empty=True),
        tbi = ensure(os.path.join(cur_dir, config['gVCF'], "{chrom}/{sample}.{chrom}.{sex}.g.vcf.gz.tbi"), non_empty=True),
    log:
        HaplotypeCaller=config['LOG'] + "/{sample}_{chrom}_{sex}_haplotypecaller.log"
    benchmark:
        config['BENCH'] + "/{sample}_{chrom}_{sex}_haplotypecaller.txt"
    conda: "envs/preprocess.yaml"
    resources: 
               n=2,
               mem_mb = get_mem_mb_HaplotypeCaller,
               tmpdir = tmpdir_alternative
    params:
        dbsnp = config['RES'] + config['dbsnp'],
        padding=300,  # extend intervals to this bp
        contam_frac = read_contam_w, # get contamination fraction per sample
        # command to get path to capture_kit interval list from SAMPLEFILE
        interval= get_chrom_merged_capture_kit,
        ploidy = get_chrom_ploidy,
        java_options=config['DEFAULT_JAVA_OPTIONS']
    priority: 28
    shell:
        """ 
                {gatk} HaplotypeCaller  --java-options "-Xmx{resources.mem_mb}M  {params.java_options}"   \
                 -R {ref} -L {params.interval} -ip {params.padding} -D {params.dbsnp} -ERC GVCF --contamination {params.contam_frac} \
                 --ploidy {params.ploidy} -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
                 -I {input.bams} -O {output.gvcf}  --native-pair-hmm-threads {resources.n}  --create-output-variant-index true\
                  --dragen-mode true --dragstr-params-path {input.model} 2> {log.HaplotypeCaller}"""


def get_mem_mb_reblock_gvcf(wildcrads, attempt):
    return attempt * int(2500)

rule reblock_gvcf:
    input:
        gvcf = rules.HaplotypeCaller.output.gvcf,
        idx = rules.HaplotypeCaller.output.tbi
    output: gvcf_reblock = ensure(os.path.join(cur_dir, config['gVCF'], "reblock/{chrom}/{sample}.{chrom}.{sex}.g.vcf.gz"), non_empty=True),
            tbi = ensure(os.path.join(cur_dir, config['gVCF'], "reblock/{chrom}/{sample}.{chrom}.{sex}.g.vcf.gz.tbi"), non_empty=True)
    log: Reblock=config['LOG'] + "/{sample}_{chrom}_{sex}_reblock.log"
    benchmark:
        config['BENCH'] + "/{sample}_{chrom}_{sex}_reblock.txt"
    conda: "envs/preprocess.yaml"
    priority: 29
    params:
        dbsnp=config['RES'] + config['dbsnp'],
        java_options=config['DEFAULT_JAVA_OPTIONS']
    resources: 
        n=2,
        mem_mb = get_mem_mb_reblock_gvcf
    shell:
        """{gatk} ReblockGVCF  --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" --keep-all-alts --create-output-variant-index true -D {params.dbsnp} -R {ref} -V {input.gvcf} -O {output.gvcf_reblock} -GQB 3 -GQB 5 -GQB 8 -GQB 10 -GQB 15 -GQB 20 -GQB 30 -GQB 50 -GQB 70 -GQB 100 -G StandardAnnotation -G AS_StandardAnnotation 2> {log}"""

