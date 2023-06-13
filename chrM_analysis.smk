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

ref = config['RES'] + config['ref']

tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())
tmpdir_alternative = os.path.join(config['tmpdir'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

wildcard_constraints:
    sample="[\w\d_\-@]+",


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

mode = config.get("computing_mode", "WES")
cur_dir = os.getcwd()

def get_refdir_by_sex(wildcards):
    if wildcards['sex'] == 'female':
        ref_dir=os.path.join(config['RES'],config['ref_female_dir'])
    else:
        ref_dir=os.path.join(config['RES'],config['ref_male_dir'])

    return ref_dir

sample_names = SAMPLEINFO.keys()


rule chrM_analysis_all:
    input:
        rules.Aligner_all.input,
    # os.path.join(config['chrM'],'variants','{sample}.chrM_filtred.vcf.gz'
        expand("{chrM}/variants/{sample}.chrM_filtred.vcf.gz", chrM = config['chrM'], sample=sample_names)
    default_target: True

rule extract_chrM_reads:
    input: rules.markdup.output.mdbams
    output: bam = ensure(os.path.join(config['chrM'], '{sample}_chrM.reads.bam'), non_empty = True),
            bai = os.path.join(config['chrM'], '{sample}_chrM.reads.bai')
    benchmark: os.path.join(config['BENCH'], '{sample}.printreads_chrM.txt')
    log: os.path.join(config['LOG'],"{sample}.printreads_chrM.log")
    conda: "envs/preprocess.yaml"
    shell:
        "gatk PrintReads -I {input} -L chrM -O {output.bam} --read-filter NotDuplicateReadFilter 2> {log}"

rule sort_by_name:
    input: rules.extract_chrM_reads.output.bam
    output: ensure(temp(os.path.join(config['chrM'], '{sample}_chrM_namesorted.reads.bam')), non_empty = True)
    benchmark: os.path.join(config['BENCH'], '{sample}.namesort_chrM.txt')
    log: os.path.join(config['LOG'],"{sample}.namesort_chrM.log")
    conda: "envs/preprocess.yaml"
    shell: "samtools sort -n {input} > {output} 2> {log}"


rule realign_to_shifted_ref:
    input: rules.sort_by_name.output
    output: bam_shifted = temp(os.path.join(config['chrM'], '{sample}_chrM_shifted.reads.bam')),
            bai_shifted= temp(os.path.join(config['chrM'],'{sample}_chrM_shifted.reads.bai'))
    conda: "envs/preprocess.yaml"
    params:
        ref_dir=os.path.join(config['RES'], config['SHIFTED_MT']),
        dragmap=os.path.join(config['RES'], config['SOFTWARE'],'dragen-os'),
    log: os.path.join(config['LOG'],"{sample}.shiftedchrM_align.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.shiftedchrM_align.txt')
    resources: n = 12,
                mem_mb = 22000
    shell: "{params.dragmap} -r {params.ref_dir} -b {input} --interleaved | samtools view -@ {resources.n} -o {output.bam_shifted} 2> {log}"

rule mutect_orig:
    input: rules.extract_chrM_reads.output.bam
    output: vcf = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.chrM_orig.vcf.gz')), non_empty = True),
            idx = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.chrM_orig.vcf.gz.idx')), non_empty = True)
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.mutect_orig.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.mutect_orig.txt')
    params: mt_ref = os.path.join(config['RES'], config['ORIG_MT_fa'])
    shell: "gatk Mutect2 -R {params.mt_ref} -L chrM --mitochondria-mode -I {input} -O {output.vcf}"

rule mutect_shifted:
    input: rules.realign_to_shifted_ref.output.bam_shifted
    output: vcf = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.chrM_shifted.vcf.gz')), non_empty = True),
            idx = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.chrM_shifted.vcf.gz.idx')), non_empty = True),
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.mutect_shift.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.mutect_shift.txt')
    params: mt_ref_shift = os.path.join(config['RES'], config['SHIFTED_MT_fa'])
    shell: "gatk Mutect2 -R {params.mt_ref_shift} -L chrM --mitochondria-mode -I {input} -O {output.vcf}"


rule shift_back:
    input: rules.mutect_shifted.output.vcf
    output: vcf = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.chrM_shifted_backshifted.vcf.gz')), non_empty = True),
            idx = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.chrM_shifted_backshifted.vcf.gz.idx')), non_empty = True),
    params: chain = config['MT_CHAIN'],
            mt_ref= os.path.join(config['RES'],config['ORIG_MT_fa'])
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.mutect_shift_back.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.mutect_shift_back.txt')
    shell: "gatk LiftoverVcf -I {input} -O {output.vcf} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null "

rule merge_vcfs:
    input: o_vcf = rules.mutect_orig.output.vcf,
            sb_vcf = rules.shift_back.output.vcf
    output: merged_vcf = ensure(os.path.join(config['chrM'], 'variants', '{sample}.chrM_merged.vcf.gz'), non_empty = True),
            merged_idx = ensure(os.path.join(config['chrM'], 'variants', '{sample}.chrM_merged.vcf.gz.idx'), non_empty = True),
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.vcf_merge.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.vcf_merge.txt')
    shell: "gatk MergeVcfs -I {input.o_vcf} -I {input.sb_vcf} -O {output}"

rule filter_mutect_calls:
    input: rules.merge_vcfs.output.merged_vcf
    output: filtred_vcf = ensure(os.path.join(config['chrM'], 'variants', '{sample}.chrM_filtred.vcf.gz'), non_empty = True),
            filtred_idx = ensure(os.path.join(config['chrM'], 'variants', '{sample}.chrM_filtred.vcf.gz.idx'), non_empty = True),
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.filtr.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.filtr.txt')
    params: mt_ref= os.path.join(config['RES'],config['ORIG_MT_fa'])
    shell: "gatk FilterMutectCalls -V {input} -R {params.mt_ref} --mitochondria-mode True -O {output.filtred_vcf}"

######################################
#####   Process NUMTs regions   ######
######################################
    
rule extract_NUMTs_reads:
    input: rules.markdup.output.mdbams
    output: bam = ensure(os.path.join(config['chrM'], '{sample}_NUMTs.reads.bam'), non_empty = True),
            bai = os.path.join(config['chrM'], '{sample}_NUMTs.reads.bai')
    benchmark: os.path.join(config['BENCH'], '{sample}.printreads_NUMTs.txt')
    log: os.path.join(config['LOG'],"{sample}.printreads_NUMTs.log")
    conda: "envs/preprocess.yaml"
    params: NUMTs_bed = os.path.join(config['RES'], config ['NUMTs'])
    shell:
        "gatk PrintReads -I {input} -L {params.NUMTs_bed} -L chrM -O {output.bam} --read-filter NotDuplicateReadFilter 2> {log}"

rule sort_by_name_NUMT:
    input: rules.extract_NUMTs_reads.output.bam
    output: ensure(temp(os.path.join(config['chrM'], "NUMT", '{sample}_NUMTs_namesorted.reads.bam')), non_empty = True)
    benchmark: os.path.join(config['BENCH'], '{sample}.namesort_NUMTs.txt')
    log: os.path.join(config['LOG'],"{sample}.namesort_NUMTs.log")
    conda: "envs/preprocess.yaml"
    shell: "samtools sort -n {input} > {output} 2> {log}"


rule align_NUMT_to_chrM:
    input: rules.sort_by_name_NUMT.output
    output: bam = temp(os.path.join(config['chrM'], "NUMT", '{sample}_NUMT_shifted.reads.bam')),
            bai = temp(os.path.join(config['chrM'], "NUMT", '{sample}_NUMTS_shifted.reads.bai'))
    conda: "envs/preprocess.yaml"
    params:
        ref_dir=os.path.join(config['RES'], config['ORIG_MT']),
        dragmap=os.path.join(config['RES'], config['SOFTWARE'],'dragen-os'),
    log: os.path.join(config['LOG'],"{sample}.origchrM_NUMT_align.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.origchrM_NUMT_align.txt')
    resources: n = 12,
                mem_mb = 22000
    shell: "{params.dragmap} -r {params.ref_dir} -b {input} --interleaved | samtools view -@ {resources.n} -o {output.bam_shifted} 2> {log}"

rule realign_to_shifted_ref_NUMT:
    input: rules.sort_by_name_NUMT.output
    output: bam_shifted = temp(os.path.join(config['chrM'], "NUMT", '{sample}_NUMT_shifted.reads.bam')),
            bai_shifted= temp(os.path.join(config['chrM'], "NUMT", '{sample}_NUMTS_shifted.reads.bai'))
    conda: "envs/preprocess.yaml"
    params:
        ref_dir=os.path.join(config['RES'], config['SHIFTED_MT']),
        dragmap=os.path.join(config['RES'], config['SOFTWARE'],'dragen-os'),
    log: os.path.join(config['LOG'],"{sample}.shiftedchrM_NUMT_align.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.shiftedchrM_NUMT_align.txt')
    resources: n = 12,
                mem_mb = 22000
    shell: "{params.dragmap} -r {params.ref_dir} -b {input} --interleaved | samtools view -@ {resources.n} -o {output.bam_shifted} 2> {log}"

rule mutect_orig_NUMT:
    input: rules.align_NUMT_to_chrM.output.bam
    output: vcf = ensure(temp(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_NUMT_orig.vcf.gz')), non_empty = True),
            idx = ensure(temp(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_NUMT_orig.vcf.gz.idx')), non_empty = True)
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.mutect_orig_NUMT.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.mutect_orig_NUMT.txt')
    params: mt_ref = os.path.join(config['RES'], config['ORIG_MT_fa'])
    shell: "gatk Mutect2 -R {params.mt_ref} -L chrM --mitochondria-mode -I {input} -O {output.vcf}"

rule mutect_shifted_NUMT:
    input: rules.realign_to_shifted_ref_NUMT.output.bam_shifted
    output: vcf = ensure(temp(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_NUMT_shifted.vcf.gz')), non_empty = True),
            idx = ensure(temp(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_NUMT_shifted.vcf.gz.idx')), non_empty = True),
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.mutect_shift_NUMT.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.mutect_shift_NUMT.txt')
    params: mt_ref_shift = os.path.join(config['RES'], config['SHIFTED_MT_fa'])
    shell: "gatk Mutect2 -R {params.mt_ref_shift} -L chrM --mitochondria-mode -I {input} -O {output.vcf}"


rule shift_back_NUMT:
    input: rules.mutect_shifted_NUMT.output.vcf
    output: vcf = ensure(temp(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_NUMT_shifted_backshifted.vcf.gz')), non_empty = True),
            idx = ensure(temp(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_NUMT_shifted_backshifted.vcf.gz.idx')), non_empty = True),
    params: chain = config['MT_CHAIN'],
            mt_ref= os.path.join(config['RES'],config['ORIG_MT_fa'])
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.mutect_shift_back_NUMT.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.mutect_shift_back_NUMT.txt')
    shell: "gatk LiftoverVcf -I {input} -O {output.vcf} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null "

rule merge_vcfs_NUMT:
    input: o_vcf = rules.mutect_orig_NUMT.output.vcf,
            sb_vcf = rules.shift_back_NUMT.output.vcf
    output: merged_vcf = ensure(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_NUMT_merged.vcf.gz'), non_empty = True),
            merged_idx = ensure(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_NUMT_merged.vcf.gz.idx'), non_empty = True),
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.vcf_merge_NUMT.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.vcf_merge_NUMT.txt')
    shell: "gatk MergeVcfs -I {input.o_vcf} -I {input.sb_vcf} -O {output}"

rule filter_mutect_calls_NUMT:
    input: rules.merge_vcfs_NUMT.output.merged_vcf
    output: filtred_vcf = ensure(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_NUMTs_filtred.vcf.gz'), non_empty = True),
            filtred_idx = ensure(os.path.join(config['chrM'], 'variants', 'NUMTs', '{sample}.chrM_filtred.vcf.gz.idx'), non_empty = True),
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.filtr_NUMT.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.filtr_NUMT.txt')
    params: mt_ref= os.path.join(config['RES'],config['ORIG_MT_fa'])
    shell: "gatk FilterMutectCalls -V {input} -R {params.mt_ref} --mitochondria-mode True -O {output.filtred_vcf}"



#
# CN estimation is basically done by calculation coverage across chrM and autosomal chr.
#     MT_CN(wgs) = 2*(mean mtDNA cov/haploid autosomal coverage)
#     mtDNA cov could be calculated by bedtools
#     autosomal cov calculated by GATK (collectwgsmetrics?)



