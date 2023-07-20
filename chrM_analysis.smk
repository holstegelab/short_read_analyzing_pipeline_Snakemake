import pandas as pd
import read_stats
import os
import getpass
import utils
from common import *


configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")


tmpdir = pj(config['TMPDIR'],getpass.getuser())
tmpdir_alternative = pj(config['tmpdir'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

wildcard_constraints:
    sample="[\w\d_\-@]+",



SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner
module Reference_preparation:
    snakefile: "Reference_preparation.smk"
    config: config
# use rule * from Reference_preparation

mode = config.get("computing_mode", "WES")
cur_dir = os.getcwd()


sample_names = SAMPLEINFO.keys()


rule chrM_analysis_all:
    input:
        rules.Aligner_all.input,
    # pj(chrM,'variants','{sample}.chrM_filtred.vcf.gz'
        expand("{chrM}/variants/{sample}.chrM_filtred_NORM.vcf.gz", chrM = chrM, sample=sample_names),
        expand("{chrM}/variants/gvcf/{sample}.chrM_merged_BP_annotated.g.vcf.gz", chrM = chrM, sample=sample_names),
        expand("{chrM}/variants/NUMTs/{sample}.chrM_NUMTs_filtred_NORM.vcf.gz", chrM = chrM, sample=sample_names),
        expand("{chrM}/variants/NUMTs/gVCF/{sample}.chrM_NUMT_merged_with_anno.g.vcf.gz", chrM = chrM, sample=sample_names),
    default_target: True

rule extract_chrM_reads:
    input: rules.markdup.output.mdbams
    output: bam = ensure(temp(pj(chrM, '{sample}_chrM.reads.bam')), non_empty = True),
            bai = ensure(temp(pj(chrM, '{sample}_chrM.reads.bai')), non_empty = True)
    benchmark: pj(BENCH, '{sample}.printreads_chrM.txt')
    log: pj(LOG,"{sample}.printreads_chrM.log")
    conda: "envs/gatk.yaml"
    shell:
        "gatk PrintReads -I {input} -L chrM -O {output.bam} --read-filter NotDuplicateReadFilter 2> {log}"

rule sort_by_name:
    input: rules.extract_chrM_reads.output.bam
    output: ensure(temp(pj(chrM, '{sample}_chrM_namesorted.reads.bam')), non_empty = True)
    benchmark: pj(BENCH, '{sample}.namesort_chrM.txt')
    log: pj(LOG,"{sample}.namesort_chrM.log")
    conda: "envs/preprocess.yaml"
    shell: "samtools sort -n {input} > {output} 2> {log}"

rule realign_to_orig_ref:
    input: rules.sort_by_name.output
    output: bam = temp(pj(chrM, '{sample}_chrM_orig.reads.bam')),
            bai = temp(pj(chrM,'{sample}_chrM_orig.reads.bai')),
            bam_shifted= temp(pj(chrM,'{sample}_chrM_shifted.reads.bam')),
            bai_shifted=temp(pj(chrM,'{sample}_chrM_shifted.reads.bai'))
    conda: "envs/dragmap_chrm.yaml"
    params:
        ref_dir=pj(ORIG_MT),
        threads_per_task = 4
    log: pj(LOG,"{sample}.orig_mt_align.log")
    benchmark: pj(BENCH, '{sample}.orig_mt_align.txt')
    resources: n = 8,
                mem_mb = 14000
    shell:
        "(dragen-os -r {params.ref_dir} -b {input} --interleaved | samtools sort -O bam -@ {resources.n} -o {output.bam} && samtools index -@ {params.threads_per_task} -o {output.bai} {output.bam} 2> {log}) &"
        "dragen-os -r {params.ref_dir} -b {input} --interleaved | samtools sort -O bam -@ {resources.n} -o {output.bam_shifted} && samtools index -@ {params.threads_per_task} -o {output.bai_shifted} {output.bam_shifted} 2> {log}"

rule mutect_orig:
    input: bam = rules.realign_to_orig_ref.output.bam,
            bai = rules.realign_to_orig_ref.output.bai,
            bam_shifted = rules.realign_to_orig_ref.output.bam_shifted,
            bai_shifted = rules.realign_to_orig_ref.output.bai_shifted
    output: vcf = ensure(temp(pj(chrM, 'variants', '{sample}.chrM_orig.vcf.gz')), non_empty = True),
            tbi = ensure(temp(pj(chrM, 'variants', '{sample}.chrM_orig.vcf.gz.tbi')), non_empty = True),
            stat = ensure(temp(pj(chrM, 'variants', '{sample}.chrM_orig.vcf.gz.stats')), non_empty = True),
            vcf_shift = ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted.vcf.gz')),non_empty=True),
            tbi_shift = ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted.vcf.gz.tbi')),non_empty=True),
            stat_shift = ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted.vcf.gz.stats')),non_empty=True),
            vcf_shift_back = ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted_backshifted.vcf.gz')),non_empty=True),
            tbi_shift_back = ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted_backshifted.vcf.gz.tbi')),non_empty=True),
    conda: "envs/gatk.yaml"
    log: pj(LOG,"{sample}.mutect_orig.log")
    benchmark: pj(BENCH, '{sample}.mutect_orig.txt')
    params:
            mt_ref = pj(ORIG_MT_fa),
            mt_ref_shift= pj(SHIFTED_MT_fa),
            chain= pj(MT_CHAIN),
            variants_dir = pj(chrM, 'variants')
    resources: n = 2
    shell:
            """
                mkdir -p {params.variants_dir}
             ((gatk Mutect2 -R {params.mt_ref} -L chrM --mitochondria-mode -I {input.bam} -O {output.vcf}  2> {log})   
             gatk Mutect2 -R {params.mt_ref_shift} -L chrM:1-500 -L chrM:16069-16569 --mitochondria-mode -I {input.bam_shifted} -O {output.vcf_shift} 2>> {log}) &&
             gatk LiftoverVcf -I {output.vcf_shift} -O {output.vcf_shift_back} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null 2>> {log} 
            """
# -L chrM:500-16069
# -L chrM:1-650 -L chrM:15919-16069


rule merge_vcfs:
    input: o_vcf = rules.mutect_orig.output.vcf,
            o_tbi = rules.mutect_orig.output.tbi,
            sb_vcf = rules.mutect_orig.output.vcf_shift_back,
            sb_tbi = rules.mutect_orig.output.tbi_shift_back,
            orig= rules.mutect_orig.output.stat,
            shift  = rules.mutect_orig.output.stat_shift
    output: merged_vcf = ensure(pj(chrM, 'variants', '{sample}.chrM_merged.vcf.gz'), non_empty = True),
            merged_tbi = ensure(pj(chrM, 'variants', '{sample}.chrM_merged.vcf.gz.tbi'), non_empty = True),
            merged_stat= ensure(pj(chrM,'variants','{sample}.chrM_merged.vcf.gz.stats'),non_empty=True),
            filtred_vcf= ensure(pj(chrM,'variants','{sample}.chrM_filtred.vcf.gz'),non_empty=True),
            filtred_tbi=ensure(pj(chrM,'variants','{sample}.chrM_filtred.vcf.gz.tbi'),non_empty=True),
            filtred_norm_vcf = ensure(pj(chrM,'variants','{sample}.chrM_filtred_NORM.vcf'),non_empty=True),
            filtred_norm_vcf_gz = ensure(pj(chrM,'variants','{sample}.chrM_filtred_NORM.vcf.gz'),non_empty=True),
    conda: "envs/gatk.yaml"
    log: pj(LOG,"{sample}.vcf_merge.log")
    benchmark: pj(BENCH, '{sample}.vcf_merge.txt')
    params: mt_ref=pj(ORIG_MT_fa)
    resources: n = 2
    shell:
            """
                 gatk MergeMutectStats --stats {input.orig} --stats {input.shift} -O {output.merged_stat} 2> {log}
               gatk MergeVcfs -I {input.sb_vcf} -I {input.o_vcf} -O {output.merged_vcf} 2>> {log} && 
               gatk FilterMutectCalls  -OVI true -V {output.merged_vcf} -R {params.mt_ref} --mitochondria-mode True -O {output.filtred_vcf} 2>> {log} && 
               bcftools norm -d exact -O v -o {output.filtred_norm_vcf} {output.filtred_vcf} 2>> {log} &&
               bgzip {output.filtred_norm_vcf} && tabix {output.filtred_norm_vcf_gz}
            """
rule mutect_orig_bp_resolut:
    input: bam = rules.realign_to_orig_ref.output.bam,
            bai = rules.realign_to_orig_ref.output.bai,
            bam_shifted = rules.realign_to_orig_ref.output.bam_shifted,
            bai_shifted = rules.realign_to_orig_ref.output.bai_shifted,
            anno_file = rules.merge_vcfs.output.filtred_norm_vcf_gz
    output: vcf_BP = ensure(temp(pj(chrM, 'variants', 'gvcf', '{sample}.chrM_orig_BP.g.vcf.gz')), non_empty = True),
            tbi_BP = ensure(temp(pj(chrM, 'variants', 'gvcf', '{sample}.chrM_orig_BP.g.vcf.gz.tbi')), non_empty = True),
            stat_BP = ensure(temp(pj(chrM, 'variants', 'gvcf', '{sample}.chrM_orig_BP.g.vcf.gz.stats')), non_empty = True),
            vcf_shift_BP = ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_BP.g.vcf.gz')),non_empty=True),
            tbi_shift_BP = ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_BP.g.vcf.gz.tbi')),non_empty=True),
            stat_shift_BP = ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_BP.g.vcf.gz.stats')),non_empty=True),
            vcf_shift_back_BP = ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_backshifted_BP.g.vcf.gz')),non_empty=True),
            tbi_shift_back_BP = ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_backshifted_BP.g.vcf.gz.tbi')),non_empty=True),
            merged_vcf_BP = ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_merged_BP.g.vcf.gz')),non_empty=True),
            merged_tbi_BP = ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_merged_BP.g.vcf.gz.tbi')),non_empty=True),
            merged_vcf_BP_with_anno = ensure(pj(chrM,'variants','gvcf','{sample}.chrM_merged_BP_annotated.g.vcf.gz'),non_empty=True),
    conda: "envs/gatk.yaml"
    log: pj(LOG,"{sample}.mutect_orig_BP_res.log")
    benchmark: pj(BENCH, '{sample}.mutect_orig_BP_res.txt')
    params:
            mt_ref = pj(ORIG_MT_fa),
            mt_ref_shift= pj(SHIFTED_MT_fa),
            chain= pj(MT_CHAIN),
            variants_dir = pj(chrM, 'variants')
    resources: n = 2
    shell:
            """
                mkdir -p {params.variants_dir}
             ((gatk Mutect2 -ERC BP_RESOLUTION -R {params.mt_ref} -L chrM --mitochondria-mode -I {input.bam} -O {output.vcf_BP}  2> {log})   
              gatk Mutect2 -ERC BP_RESOLUTION -R {params.mt_ref_shift} -L chrM:1-500 -L chrM:16069-16569 --mitochondria-mode -I {input.bam_shifted} -O {output.vcf_shift_BP} 2>> {log}) &&
              gatk LiftoverVcf -I {output.vcf_shift_BP} -O {output.vcf_shift_back_BP} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null 2>> {log} && 
              gatk MergeVcfs -I {output.vcf_shift_back_BP} -I {output.vcf_BP} -O {output.merged_vcf_BP} 2>> {log} && 
              bcftools norm -d exact -O v {output.merged_vcf_BP} | bcftools annotate -a {input.anno_file} -c FILTER -O z -o {output.merged_vcf_BP_with_anno}  - 2>> {log}
            """

# (gatk MergeMutectStats --stats {input.orig} --stats {input.shift} -O {output.merged_stat} 2> {log})


######################################
#####   Process NUMTs regions   ######
######################################

rule extract_NUMTs_reads:
    input: rules.markdup.output.mdbams
    output: bam = ensure(pj(chrM, "NUMTs", '{sample}_NUMTs.reads.bam'), non_empty = True),
            bai = pj(chrM, "NUMTs", '{sample}_NUMTs.reads.bai')
    benchmark: pj(BENCH, '{sample}.printreads_NUMTs.txt')
    log: pj(LOG,"{sample}.printreads_NUMTs.log")
    conda: "envs/gatk.yaml"
    params: NUMTs_bed = NUMTs
    shell:
        "gatk PrintReads -I {input} -L {params.NUMTs_bed} -L chrM -O {output.bam} --read-filter NotDuplicateReadFilter 2> {log}"

rule sort_by_name_NUMT:
    input: rules.extract_NUMTs_reads.output.bam
    output: ensure(temp(pj(chrM, "NUMTs", '{sample}_NUMTs_namesorted.reads.bam')), non_empty = True)
    benchmark: pj(BENCH, '{sample}.namesort_NUMTs.txt')
    log: pj(LOG,"{sample}.namesort_NUMTs.log")
    conda: "envs/preprocess.yaml"
    shell: "samtools sort -n {input} > {output} 2> {log}"


rule align_NUMT_to_chrM:
    input: rules.sort_by_name_NUMT.output
    output: bam = temp(pj(chrM, "NUMTs", '{sample}_NUMTs.realign.bam')),
            bai = temp(pj(chrM, "NUMTs", '{sample}_NUMTs.realign.bai')),
            bam_shifted= temp(pj(chrM,"NUMTs",'{sample}_NUMTs_shifted.reads.bam')),
            bai_shifted=temp(pj(chrM,"NUMTs",'{sample}_NUMTs_shifted.reads.bai'))
    conda: "envs/dragmap_chrm.yaml"
    params:
        ref_dir=pj(ORIG_MT),
        threads_per_task = 4,

    log: pj(LOG,"{sample}.origchrM_NUMT_align.log")
    benchmark: pj(BENCH, '{sample}.origchrM_NUMT_align.txt')
    resources: n = 8,
                mem_mb = 14000
    shell:
            """
            dragen-os -r {params.ref_dir} -b {input} --interleaved | samtools sort -O bam -@ {params.threads_per_task} -o {output.bam} && samtools index -@ {resources.n} -o {output.bai} {output.bam} 2> {log} &
             dragen-os -r {params.ref_dir} -b {input} --interleaved | samtools sort -O bam -@ {params.threads_per_task} -o {output.bam_shifted} && samtools index -@ {resources.n} -o {output.bai_shifted} {output.bam_shifted}  2>> {log}
            """
rule mutect_orig_NUMT:
    input: bam = rules.align_NUMT_to_chrM.output.bam,
            bai = rules.align_NUMT_to_chrM.output.bai,
            bam_shifted = rules.align_NUMT_to_chrM.output.bam_shifted,
            bai_shifted = rules.align_NUMT_to_chrM.output.bai_shifted
    output: vcf = ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_orig.vcf.gz')), non_empty = True),
            tbi = ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_orig.vcf.gz.tbi')), non_empty = True),
            stat = ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_orig.vcf.gz.stats')), non_empty = True),
            vcf_shifted = ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted.vcf.gz')),non_empty=True),
            tbi_shifted = ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted.vcf.gz.tbi')),non_empty=True),
            stat_shifted = ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted.vcf.gz.stats')),non_empty=True),
            vcf_shiftback = ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted_backshifted.vcf.gz')),non_empty=True),
            tbi_shiftback = ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted_backshifted.vcf.gz.tbi')),non_empty=True),
    conda: "envs/gatk.yaml"
    log: pj(LOG,"{sample}.mutect_orig_NUMT.log")
    benchmark: pj(BENCH, '{sample}.mutect_orig_NUMT.txt')
    params: mt_ref = pj(ORIG_MT_fa),
            mt_ref_shift= pj(SHIFTED_MT_fa),
            chain= pj(MT_CHAIN),
            results_dir = pj(chrM, 'variants' , 'NUMTs')
    shell:
            """
             mkdir -p {params.results_dir}
             gatk Mutect2 -R {params.mt_ref} -L chrM --mitochondria-mode -I {input.bam} -O {output.vcf} 2> {log}   
             gatk Mutect2 -R {params.mt_ref_shift} -L chrM:1-500 -L chrM:16069-16569 --mitochondria-mode -I {input.bam_shifted} -O {output.vcf_shifted} 2>>{log} &&
             gatk LiftoverVcf -I {output.vcf_shifted} -O {output.vcf_shiftback} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null 2>> {log}
            """

rule merge_vcfs_NUMT:
    input: o_vcf = rules.mutect_orig_NUMT.output.vcf,
            o_tbi = rules.mutect_orig_NUMT.output.tbi,
            sb_vcf = rules.mutect_orig_NUMT.output.vcf_shiftback,
            sb_tbi = rules.mutect_orig_NUMT.output.tbi_shiftback,
            stats = rules.mutect_orig_NUMT.output.stat,
            stats_shifted = rules.mutect_orig_NUMT.output.stat_shifted
    output: merged_vcf = ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_merged.vcf.gz')), non_empty = True),
            merged_tbi = ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_merged.vcf.gz.tbi')), non_empty = True),
            merged_stat= ensure(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_merged.vcf.gz.stats')),
            filtred_vcf= ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMTs_filtred.vcf.gz')),non_empty=True),
            filtred_tbi = ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMTs_filtred.vcf.gz.tbi')),non_empty=True),
            filtred_norm_vcf= ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMTs_filtred_NORM.vcf')),non_empty=True),
            filtred_norm_vcf_gz = ensure(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMTs_filtred_NORM.vcf.gz'),non_empty=True),
    conda: "envs/gatk.yaml"
    log: pj(LOG,"{sample}.vcf_merge_NUMT.log")
    benchmark: pj(BENCH, '{sample}.vcf_merge_NUMT.txt')
    params: mt_ref= pj(ORIG_MT_fa)
    shell:
            """
            gatk MergeVcfs -I {input.sb_vcf} -I {input.o_vcf} -O {output.merged_vcf} 2> {log} 
             gatk MergeMutectStats --stats {input.stats} --stats {input.stats_shifted} -O {output.merged_stat} 2>> {log} &&
             gatk FilterMutectCalls  -OVI true -V {output.merged_vcf} -R {params.mt_ref} --mitochondria-mode True -O {output.filtred_vcf} 2>> {log} &&
               bcftools norm -d exact -O v -o {output.filtred_norm_vcf} {output.filtred_vcf} 2>> {log} &&
               bgzip {output.filtred_norm_vcf_gz} && tabix {output.filtred_norm_vcf_gz}
            """

rule mutect_orig_NUMT_BP_resolution:
    input: bam = rules.align_NUMT_to_chrM.output.bam,
            bai = rules.align_NUMT_to_chrM.output.bai,
            bam_shifted = rules.align_NUMT_to_chrM.output.bam_shifted,
            bai_shifted = rules.align_NUMT_to_chrM.output.bai_shifted,
            anno_file = rules.merge_vcfs_NUMT.output.filtred_norm_vcf_gz
    output: vcf = ensure(temp(pj(chrM, 'variants', 'NUMTs', 'gVCF', '{sample}.chrM_NUMT_orig_BP_res.g.vcf.gz')), non_empty = True),
            tbi = ensure(temp(pj(chrM, 'variants', 'NUMTs', 'gVCF', '{sample}.chrM_NUMT_orig_BP_res.g.vcf.gz.tbi')), non_empty = True),
            stat = ensure(temp(pj(chrM, 'variants', 'NUMTs', 'gVCF', '{sample}.chrM_NUMT_orig_BP_res.g.vcf.gz.stats')), non_empty = True),
            vcf_shifted = ensure(temp(pj(chrM,'variants','NUMTs', 'gVCF','{sample}.chrM_NUMT_shifted_BP_res.g.vcf.gz')),non_empty=True),
            tbi_shifted = ensure(temp(pj(chrM,'variants','NUMTs', 'gVCF','{sample}.chrM_NUMT_shifted_BP_res.g.vcf.gz.tbi')),non_empty=True),
            stat_shifted = ensure(temp(pj(chrM,'variants','NUMTs', 'gVCF','{sample}.chrM_NUMT_shifted_BP_res.g.vcf.gz.stats')),non_empty=True),
            vcf_shiftback = ensure(temp(pj(chrM,'variants','NUMTs', 'gVCF','{sample}.chrM_NUMT_shifted_backshifted_BP_res.g.vcf.gz')),non_empty=True),
            tbi_shiftback = ensure(temp(pj(chrM,'variants','NUMTs', 'gVCF','{sample}.chrM_NUMT_shifted_backshifted_BP_res.g.vcf.gz.tbi')),non_empty=True),
            merged_vcf = ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_merged.g.vcf.gz')),non_empty=True),
            merged_vcf_with_anno = ensure(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_merged_with_anno.g.vcf.gz'),non_empty=True),
    conda: "envs/gatk.yaml"
    log: pj(LOG,"{sample}.mutect_orig_NUMT_BP_res.log")
    benchmark: pj(BENCH, '{sample}.mutect_orig_NUMT_BP_res.txt')
    params: mt_ref = pj(ORIG_MT_fa),
            mt_ref_shift= pj(SHIFTED_MT_fa),
            chain= pj(MT_CHAIN),
            results_dir = pj(chrM, 'variants' , 'NUMTs')
    shell:
            """
             mkdir -p {params.results_dir}
             gatk Mutect2  -ERC BP_RESOLUTION -R {params.mt_ref} -L chrM --mitochondria-mode -I {input.bam} -O {output.vcf} 2> {log}   
             gatk Mutect2  -ERC BP_RESOLUTION -R {params.mt_ref_shift} -L chrM:1-500 -L chrM:16069-16569 --mitochondria-mode -I {input.bam_shifted} -O {output.vcf_shifted} 2>>{log} &&
             gatk LiftoverVcf -I {output.vcf_shifted} -O {output.vcf_shiftback} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null 2>> {log} && 
              gatk MergeVcfs -I {output.vcf_shiftback} -I {output.vcf} -O {output.merged_vcf} 2>> {log} && 
              bcftools norm -d exact -O v {output.merged_vcf} | bcftools annotate -a {input.anno_file} -c FILTER -O z -o {output.merged_vcf_with_anno}  - 2>> {log}
            """



#
# CN estimation is basically done by calculation coverage across chrM and autosomal chr.
#     MT_CN(wgs) = 2*(mean mtDNA cov/haploid autosomal coverage)
#     mtDNA cov could be calculated by bedtools
#     autosomal cov calculated by GATK (collectwgsmetrics?)
#
#
#
