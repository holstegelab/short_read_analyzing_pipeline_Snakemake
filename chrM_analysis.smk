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
sample_sex_names = get_validated_sex(SAMPLEINFO)

rule chrM_analysis_all:
    input:
        rules.Aligner_all.input,
    # os.path.join(config['chrM'],'variants','{sample}.{sex}.chrM_filtred.vcf.gz'
        expand("{chrM}/variants/{sample_sex}.chrM_filtred.vcf.gz", chrM = config['chrM'], sample_sex=sample_sex_names)
    default_target: True

rule extract_chrM_reads:
    input: rules.markdup.output.mdbams
    output: bam = ensure(os.path.join(config['chrM'], '{sample}.{sex}_chrM.reads.bam'), non_empty = True),
            bai = os.path.join(config['chrM'], '{sample}.{sex}_chrM.reads.bai')
    benchmark: os.path.join(config['BENCH'], '{sample}.{sex}.printreads_chrM.txt')
    log: os.path.join(config['LOG'],"{sample}.{sex}.printreads_chrM.log")
    conda: "envs/preprocess.yaml"
    shell:
        "gatk PrintReads -I {input} -L chrM -O {output.bam} --read-filter NotDuplicateReadFilter 2> {log}"

rule sort_by_name:
    input: rules.extract_chrM_reads.output.bam
    output: ensure(temp(os.path.join(config['chrM'], '{sample}.{sex}_chrM_namesorted.reads.bam')), non_empty = True)
    benchmark: os.path.join(config['BENCH'], '{sample}.{sex}.namesort_chrM.txt')
    log: os.path.join(config['LOG'],"{sample}.{sex}.namesort_chrM.log")
    conda: "envs/preprocess.yaml"
    shell: "samtools sort -n {input} > {output} 2> {log}"


rule realign_to_shifted_ref:
    input: rules.sort_by_name.output
    output: bam_shifted = temp(os.path.join(config['chrM'], '{sample}.{sex}_chrM_shifted.reads.bam')),
            bai_shifted= temp(os.path.join(config['chrM'],'{sample}.{sex}_chrM_shifted.reads.bai'))
    conda: "envs/preprocess.yaml"
    params:
        ref_dir=os.path.join(config['RES'], config['SHIFTED_MT']),
        dragmap=os.path.join(config['RES'], config['SOFTWARE'],'dragen-os'),
    log: os.path.join(config['LOG'],"{sample}.{sex}.shiftedchrM_align.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.{sex}.shiftedchrM_align.txt')
    resources: n = 12,
                mem_mb = 22000
    shell: "{params.dragmap} -r {params.ref_dir} -b {input} --interleaved | samtools view -@ {resources.n} -o {output.bam_shifted} 2> {log}"

rule mutect_orig:
    input: rules.extract_chrM_reads.output.bam
    output: vcf = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_orig.vcf.gz')), non_empty = True),
            idx = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_orig.vcf.gz.idx')), non_empty = True)
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.{sex}.mutect_orig.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.{sex}.mutect_orig.txt')
    params: mt_ref = os.path.join(config['RES'], config['ORIG_MT_fa'])
    shell: "gatk Mutect2 -R {params.mt_ref} -L chrM --mitochondria-mode -I {input} -O {output.vcf}"

rule mutect_shifted:
    input: rules.realign_to_shifted_ref.output.bam_shifted
    output: vcf = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_shifted.vcf.gz')), non_empty = True),
            idx = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_shifted.vcf.gz.idx')), non_empty = True),
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.{sex}.mutect_shift.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.{sex}.mutect_shift.txt')
    params: mt_ref_shift = os.path.join(config['RES'], config['SHIFTED_MT_fa'])
    shell: "gatk Mutect2 -R {params.mt_ref_shift} -L chrM --mitochondria-mode -I {input} -O {output.vcf}"


rule shift_back:
    input: rules.mutect_shifted.output.vcf
    output: vcf = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_shifted_backshifted.vcf.gz')), non_empty = True),
            idx = ensure(temp(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_shifted_backshifted.vcf.gz.idx')), non_empty = True),
    params: chain = config['MT_CHAIN'],
            mt_ref= os.path.join(config['RES'],config['ORIG_MT_fa'])
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.{sex}.mutect_shift_back.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.{sex}.mutect_shift_back.txt')
    shell: "gatk LiftoverVcf -I {input} -O {output.vcf} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null "

rule merge_vcfs:
    input: o_vcf = rules.mutect_orig.output.vcf,
            sb_vcf = rules.shift_back.output.vcf
    output: merged_vcf = ensure(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_merged.vcf.gz'), non_empty = True),
            merged_idx = ensure(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_merged.vcf.gz.idx'), non_empty = True),
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.{sex}.vcf_merge.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.{sex}.vcf_merge.txt')
    shell: "gatk MergeVcfs -I {input.o_vcf} -I {input.sb_vcf} -O {output}"

rule filter_mutect_calls:
    input: rules.merge_vcfs.output.merged_vcf
    output: filtred_vcf = ensure(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_filtred.vcf.gz'), non_empty = True),
            filtred_idx = ensure(os.path.join(config['chrM'], 'variants', '{sample}.{sex}.chrM_filtred.vcf.gz.idx'), non_empty = True),
    conda: "envs/preprocess.yaml"
    log: os.path.join(config['LOG'],"{sample}.{sex}.filtr.log")
    benchmark: os.path.join(config['BENCH'], '{sample}.{sex}.filtr.txt')
    params: mt_ref= os.path.join(config['RES'],config['ORIG_MT_fa'])
    shell: "gatk FilterMutectCalls -V {input} -R {params.mt_ref} --mitochondria-mode True -O {output.filtred_vcf}"



