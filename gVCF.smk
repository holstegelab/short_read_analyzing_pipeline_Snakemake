import pandas as pd
import read_stats
import os
import getpass
import yaml

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

module Tools:
    snakefile: 'Tools.smk'
    config: config
use rule * from Tools

mode = config.get("computing_mode", "WES")




rule gVCF_all:
    input:
        expand("{gvcf}/{sample}.done",sample=sample_names, gvcf = config['gVCF']),
        rules.Aligner_all.input
    default_target: True

def get_gvcf_files(wildcards):
    sample = wildcards['sample']
    res = []
    for chrom in main_chrs_ploidy_male:
        if 'Y' in chrom or 'H' in chrom:
            res.append(os.path.join(config['gVCF'], chrom, sample + '.' + chrom + '.g.vcf.gz'))
        else:        
            res.append(os.path.join(config['gVCF'], chrom, sample + '.' + chrom + '.wg.vcf.gz'))
    return res

rule gvcf_sample_done:
    input:
        get_gvcf_files
    output:
        touch(os.path.join(config['gVCF'], "{sample}.done"))




rule ComposeSTRTableFile:
    input:
        ref_fasta = "{path}.fa"
    output:
        str_table = "{path}.str.zip"
    params:
        java_options=config['DEFAULT_JAVA_OPTIONS']
    resources: 
        n = 1,
        mem_mb = 3000
    conda: "envs/gatk.yaml"        
    shell:
        """gatk ComposeSTRTableFile --java-options "-Xmx{resources.mem_mb}m {params.java_options}" -R {input.ref_fasta} -O {output.str_table}"""



 
def get_strref_by_validated_sex(wildcards, input):
    sex = get_validated_sex_file(input)
    if sex == 'female':
        ref=os.path.join(config['RES'],config['ref_female_str'])
    else:
        ref=os.path.join(config['RES'],config['ref_male_str'])

    return ref

def get_ref_by_validated_sex(wildcards, input):
    sex = get_validated_sex_file(input)
    if sex == 'female':
        ref=os.path.join(config['RES'],config['ref_female'])
    else:
        ref=os.path.join(config['RES'],config['ref_male'])

    return ref

 
def get_mem_mb_CalibrateDragstrModel(wildcards, attempt):
    return attempt * int(1750)

rule CalibrateDragstrModel:
    """CalibrateDragstrModel. Estimates the parameters of the dragstr model from a set of aligned reads.
    Provides better STR variant calling.
    """
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
        validated_sex = rules.get_validated_sex.output.yaml       
    output:
        dragstr_model = config['BAM'] + "/{sample}-dragstr.txt"
    priority: 26
    params:
        ref=get_ref_by_validated_sex,
        str_ref = get_strref_by_validated_sex,
        java_options=config['DEFAULT_JAVA_OPTIONS']
    resources: 
        n = 1,
        mem_mb = get_mem_mb_CalibrateDragstrModel
    log: config['LOG'] + '/' + "{sample}_calibratedragstr.log"
    benchmark: config['BENCH'] + "/{sample}_calibrate_dragstr.txt"
    conda: 'envs/gatk.yaml'
    shell:
        """{gatk} CalibrateDragstrModel --java-options \
                    "-Xmx{resources.mem_mb}m {params.java_options}"  -R {params.ref} -I {input.bam} \
                    -O {output} -str {params.str_ref} 2>{log}"""

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
        capture_kit_chr_path = os.path.join(config['RES'], config['kit_folder'], capture_kit,  capture_kit + '_chrom_' + chrom + '.interval_list')

    return capture_kit_chr_path



def read_contam_w(wildcards):
    """Read contamination from verifybamid output file."""
    filename = os.path.join(config['STAT'], 'contam', wildcards['sample'] +  '.verifybamid.pca2.selfSM')
    with open(filename,'r', encoding='utf-8') as f:
        c = csv.reader(f, delimiter='\t')
        c = iter(c)
        headers = c.__next__()
        data = c.__next__()
        freemix = float(data[6])
    if freemix < 0.01: #only remove contamination if it is more than 1%
        freemix = 0.0        
    return freemix

def get_mem_mb_HaplotypeCaller(wildcards, attempt):
    """Get memory for HaplotypeCaller."""
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        res = 1750
    else:
        res = 1400
    return attempt * res


def get_chrom_merged_capture_kit(wildcards):
    chrom = wildcards.chrom
    chrom = chrom.replace('XH', 'X').replace('YH', 'Y')
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        capture_kit_chr_path = []
    else:
        capture_kit_chr_path = [os.path.join(config['RES'], config['kit_folder'], config['MERGED_CAPTURE_KIT'], config['MERGED_CAPTURE_KIT'] + '_chrom_' + chrom + '.interval_list')]
    
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
        bai = rules.markdup.output.mdbams_bai,
        model = rules.CalibrateDragstrModel.output.dragstr_model,
        contam= rules.verifybamid.output.VBID_stat,
        interval= get_chrom_merged_capture_kit,
        validated_sex = rules.get_validated_sex.output.yaml
    output:
        gvcf= ensure(os.path.join(config['gVCF'], "{chrom}/{sample}.{chrom}.g.vcf.gz"), non_empty=True),
        tbi = ensure(os.path.join(config['gVCF'], "{chrom}/{sample}.{chrom}.g.vcf.gz.tbi"), non_empty=True),
    log:
        HaplotypeCaller=config['LOG'] + "/{sample}_{chrom}_haplotypecaller.log"
    benchmark:
        config['BENCH'] + "/{sample}_{chrom}_haplotypecaller.txt"
    conda: "envs/gatk.yaml"
    resources: 
               n=1, #average 1.3 cores
               mem_mb = get_mem_mb_HaplotypeCaller,
               tmpdir = tmpdir_alternative
    params:
        interval = lambda wildcards, input: wildcards['chrom'].replace('YH', 'Y').replace('XH','X') if len(input.interval) == 0 else input.interval[0],
        dbsnp = config['RES'] + config['dbsnp'],
        padding=300,  # extend intervals to this bp
        contam_frac = read_contam_w, # get contamination fraction per sample
        # command to get path to capture_kit interval list from SAMPLEFILE
        ploidy = get_chrom_ploidy,
        java_options=config['DEFAULT_JAVA_OPTIONS'],
        ref=get_ref_by_validated_sex,
        dragen_mode = lambda wildcards: '--dragen-mode true' if not 'H' in wildcards['chrom'] else ''
    priority: 28
    #Overview of available annotations (version 4.4)
    #AS_StandardAnnotation: AS_BaseQualityRankSumTest, AS_FisherStrand, AS_InbreedingCoeff, AS_MappingQualityRankSumTest, AS_QualByDepth, AS_RMSMappingQuality, AS_ReadPosRankSumTest, AS_StrandOddsRatio (all allele specific annotations)
    #StandardHCAnnotation: DepthPerSampleHC
    #StandardANnotation: Coverage, ChromosomeCounts, BaseQualityRankSumTest, MappingQualityRankSumTest, FisherStrand, QualByDepth, InbreedingCoeff, DepthPerAlleleBySample, ReadPosRankSumTest, ExcessHet
    #                      StrandOddsRatio, RMSMappingQuality, 
    #not included in standard annotation: FragmentDepthPerAlleleBySample,StrandBiasBySample, AssemblyComplexity,  FragmentLength, OrientationBiasReadCounts, TandemRepeat, UniqueAltReadCount
    
    #findings:
    #FragmentDepthPerAlleleBySample, OrientationBiasReadCounts, TandemRepeat --> do not output anything in haplotypecaller, omitted

    #regarding linked-de-bruijn-graph, see comments in https://gatk.broadinstitute.org/hc/en-us/community/posts/360077647812-Why-do-a-clear-expected-variant-not-show-up-in-the-Mutect2-vcf-file
    #also: https://gatk.broadinstitute.org/hc/en-us/articles/360043491652-When-HaplotypeCaller-and-Mutect2-do-not-call-an-expected-variant
    #still, seems quite experimental, so not using it for now
    
    #related: --recover-all-dangling-branches
    #to consider: --population-callset in combination with gnomad,  --use-posteriors-to-calculate-qual

    #max-mnp-distance: not compatible with multiple sample gVCF calling, omitted.
    shell:
        """ 
                {gatk} --java-options "-Xmx6000M  {params.java_options}" HaplotypeCaller     \
                 -R {params.ref} -L {params.interval} -ip {params.padding} -D {params.dbsnp} -ERC GVCF --contamination {params.contam_frac} \
                 --ploidy {params.ploidy} -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
                 --annotate-with-num-discovered-alleles --adaptive-pruning \
                 -A StrandBiasBySample -A AssemblyComplexity -A FragmentLength \
                 -I {input.bams} -O {output.gvcf}  --native-pair-hmm-threads 2  --create-output-variant-index true\
                  {params.dragen_mode} --dragstr-params-path {input.model} 2> {log.HaplotypeCaller}"""


rule PrepareWhatshap:
    """Prepare input for whatshap phasing."""
    input:
        gvcf = rules.HaplotypeCaller.output.gvcf,
        tbi = rules.HaplotypeCaller.output.tbi,
        validated_sex = rules.get_validated_sex.output.yaml
    output:
        vcf = temp(os.path.join(config['gVCF'], "{chrom}/{sample}.{chrom}.vcf.gz")),
        vcf_tbi = temp(os.path.join(config['gVCF'], "{chrom}/{sample}.{chrom}.vcf.gz.tbi"))
    params:
        java_options=config['DEFAULT_JAVA_OPTIONS'],
        ref=get_ref_by_validated_sex,
    resources: 
        n=1,
        mem_mb = 2500
    conda: "envs/gatk.yaml"        
    shell: """
        {gatk} --java-options "-Xmx4000M  {params.java_options}" GenotypeGVCFs \
            -R {ref} -V {input.gvcf} -O {output.vcf} 
        
        """        

rule WhatshapPhasingMerge:
    """Phase VCF with Whatshap and merge into the gVCF"""
    input:
        vcf = rules.PrepareWhatshap.output.vcf,
        vcf_tbi = rules.PrepareWhatshap.output.vcf_tbi,
        gvcf = rules.HaplotypeCaller.output.gvcf,
        gvcf_tbi = rules.HaplotypeCaller.output.tbi,
        bams = rules.markdup.output.mdbams,
        bai = rules.markdup.output.mdbams_bai,
        validated_sex = rules.get_validated_sex.output.yaml
    output:
        tmp_vcf = temp(os.path.join(config['gVCF'], "{chrom}/{sample}.{chrom}.tmp.vcf.gz")),
        vcf = temp(os.path.join(config['gVCF'], "{chrom}/{sample}.{chrom}.w.vcf.gz")),
        wstats = os.path.join(config['STAT'], "whatshap_phasing/{sample}.{chrom}.stats"),
        mwstats = os.path.join(config['STAT'], "whatshap_phasing/{sample}.{chrom}.merge_stats"),
        tmp_gvcf= temp(os.path.join(config['gVCF'], "{chrom}/{sample}.{chrom}.wg.vcf")), 
        gvcf= os.path.join(config['gVCF'], "{chrom}/{sample}.{chrom}.wg.vcf.gz"), 
        gvcf_tbi = os.path.join(config['gVCF'], "{chrom}/{sample}.{chrom}.wg.vcf.gz.tbi"), 
    params:
        merge_script=srcdir("scripts/merge_phasing.py"),
        ref=get_ref_by_validated_sex
    resources: 
        n=1,
        mem_mb = 2500
    conda: "envs/whatshap.yaml"
    shell: """
        mkdir -p `dirname {output.wstats}`
        whatshap unphase  {input.vcf} >  {output.tmp_vcf}
        whatshap phase  --ignore-read-groups --reference {params.ref} {output.tmp_vcf} {input.bams} -o {output.vcf}
        whatshap stats {input.vcf} > {output.wstats}
        whatshap stats {output.vcf} >> {output.wstats}
        python {params.merge_script} {input.gvcf} {output.vcf} {output.tmp_gvcf} {output.mwstats}
        bcftools annotate -x 'FORMAT/PGT,FORMAT/PS,FORMAT/PID' {output.tmp_gvcf} -o {output.gvcf}
        bcftools index --tbi {output.gvcf}
        """        


def get_mem_mb_reblock_gvcf(wildcrads, attempt):
    return attempt * int(2500)

rule reblock_gvcf:
    input:
        gvcf = rules.HaplotypeCaller.output.gvcf,
        idx = rules.HaplotypeCaller.output.tbi
    output: gvcf_reblock = ensure(os.path.join(config['gVCF'], "reblock/{chrom}/{sample}.{chrom}.g.vcf.gz"), non_empty=True),
            tbi = ensure(os.path.join(config['gVCF'], "reblock/{chrom}/{sample}.{chrom}.g.vcf.gz.tbi"), non_empty=True)
    log: Reblock=config['LOG'] + "/{sample}_{chrom}_reblock.log"
    benchmark:
        config['BENCH'] + "/{sample}_{chrom}_reblock.txt"
    conda: "envs/gatk.yaml"
    priority: 29
    params:
        dbsnp=config['RES'] + config['dbsnp'],
        java_options=config['DEFAULT_JAVA_OPTIONS']
    resources: 
        n=2,
        mem_mb = get_mem_mb_reblock_gvcf
    shell:
        """{gatk} --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" ReblockGVCF   --keep-all-alts --create-output-variant-index true -D {params.dbsnp} -R {ref} -V {input.gvcf} -O {output.gvcf_reblock} -GQB 3 -GQB 5 -GQB 8 -GQB 10 -GQB 15 -GQB 20 -GQB 30 -GQB 50 -GQB 70 -GQB 100 -G StandardAnnotation -G AS_StandardAnnotation 2> {log}"""

