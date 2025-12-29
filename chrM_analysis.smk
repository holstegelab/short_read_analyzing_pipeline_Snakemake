import pandas as pd
import read_stats
import os
import getpass
import utils
from shlex import quote
from common import *
onsuccess: shell("rm -fr logs/chrM/*")

wildcard_constraints:
    sample=r"[\w\d_\-@]+",



module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner
module Reference_preparation:
    snakefile: "Reference_preparation.smk"
    config: config

mode = config.get("computing_mode", "WES")
cur_dir = os.getcwd()


rule chrM_analysis_all:
    input:
        rules.Aligner_all.input,
        expand("{chrM}/variants/gvcf/{sample}.chrM_merged_BP_annotated.g.vcf.gz", chrM = chrM, sample=sample_names),
        pj(chrM, "chrM_tar_uploads.done"),


def chrM_gvcf_inputs(wildcards):
    files = []
    for sample in sample_names:
        base = pj(chrM, "variants", "gvcf", f"{sample}.chrM_merged_BP_annotated.g.vcf.gz")
        files.append(base)
        files.append(base + ".tbi")
    return files


rule tar_chrM_gvcfs:
    input:
        chrM_gvcf_inputs
    output:
        tar=temp(pj(chrM, "tar", "chrM_gvcfs.tar.gz"))
    params:
        files=lambda wildcards, input: " ".join(quote(path) for path in input),
        outdir=lambda wildcards, output: quote(os.path.dirname(output.tar))
    resources:
        mem_mb=2000,
        n="0.5"
    shell:
        """
        mkdir -p {params.outdir}
        tar -czf {output.tar:q} {params.files}
        """


rule copy_chrM_gvcfs_to_dcache:
    input:
        tar=pj(chrM, "tar", "chrM_gvcfs.tar.gz")
    output:
        copied=temp(pj(chrM, "tar", "chrM_gvcfs.tar.copied")),
        checksum=temp(pj(chrM, "tar", "chrM_gvcfs.tar.ADLER32"))
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1",
        dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        if not SAMPLE_FILES:
            raise ValueError("No sample files available to resolve remote destination for chrM gVCFs")
        samplefile = next(iter(SAMPLE_FILES))
        remote_dir = os.path.join(remote_base_for_samplefile(samplefile), "chrM")
        remote_name = os.path.basename(input.tar)
        copy_with_checksum(str(input.tar), remote_dir, remote_name, str(output.checksum), AGH_DCACHE_CONFIG, params.ada_script)
        shell(f"touch {quote(str(output.copied))}")


rule chrM_tar_all:
    input:
        pj(chrM, "tar", "chrM_gvcfs.tar.copied")
    output:
        done=touch(pj(chrM, "chrM_tar_uploads.done"))

rule chrM_sample_done:
    input:
        gvcf=pj(chrM, "variants", "gvcf", "{sample}.chrM_merged_BP_annotated.g.vcf.gz"),
        numt=pj(chrM, "variants", "NUMTs", "gVCF", "{sample}.chrM_NUMT_merged_with_anno.g.vcf.gz")
    output:
        done=touch(pj(chrM, "{sample}.done"))

rule extract_chrM_reads:
    input: rules.markdup.output.mdbams
    output:
        fq1 = ensure(temp(pj(chrM, '{sample}_chrM.R1.fastq.gz')), non_empty = True),
        fq2 = ensure(temp(pj(chrM, '{sample}_chrM.R2.fastq.gz')), non_empty = True)
    conda: CONDA_VCF
    resources:
        mem_mb=1000
    params:
        threads=2
    shell:
        """
        tmp_dir=$(mktemp -d)
        trap 'rm -rf "$tmp_dir"' EXIT
        samtools view -@ {params.threads} -b -o "$tmp_dir/{wildcards.sample}.chrM.raw.bam" {input} chrM
        mkdir -p $(dirname {output.fq1})
        samtools sort -n -@ {params.threads} -o "$tmp_dir/{wildcards.sample}.chrM.sorted.bam" "$tmp_dir/{wildcards.sample}.chrM.raw.bam"
        samtools collate -@ {params.threads} -o "$tmp_dir/{wildcards.sample}.chrM.collated.bam" "$tmp_dir/{wildcards.sample}.chrM.sorted.bam"
        samtools fastq -O -N -@ {params.threads} -0 /dev/null -s /dev/null -1 {output.fq1} -2 {output.fq2} "$tmp_dir/{wildcards.sample}.chrM.collated.bam"
        """
 
rule extract_NUMTs_reads:
    input: rules.markdup.output.mdbams
    output:
        fq1 = ensure(temp(pj(chrM, "NUMTs", '{sample}_NUMTs.R1.fastq.gz')), non_empty = True),
        fq2 = ensure(temp(pj(chrM, "NUMTs", '{sample}_NUMTs.R2.fastq.gz')), non_empty = True)
    conda: CONDA_VCF
    params:
        NUMTs_bed = NUMTs,
        threads = 2
    resources:
        mem_mb=1000
    shell:
        """
        tmp_dir=$(mktemp -d)
        trap 'rm -rf "$tmp_dir"' EXIT
        tmpbed="$tmp_dir/NUMTs_plus_chrM.bed"
        cat {params.NUMTs_bed} > "$tmpbed"
        echo -e "chrM\t1\t999999999" >> "$tmpbed"
        samtools view -@ {params.threads} -b -L "$tmpbed" -o "$tmp_dir/{wildcards.sample}.NUMTs.raw.bam" {input}
        mkdir -p $(dirname {output.fq1})
        samtools sort -n -@ {params.threads} -o "$tmp_dir/{wildcards.sample}.NUMTs.sorted.bam" "$tmp_dir/{wildcards.sample}.NUMTs.raw.bam"
        samtools collate -@ {params.threads} -o "$tmp_dir/{wildcards.sample}.NUMTs.collated.bam" "$tmp_dir/{wildcards.sample}.NUMTs.sorted.bam"
        samtools fastq -O -N -@ {params.threads} -0 /dev/null -s /dev/null -1 {output.fq1} -2 {output.fq2} "$tmp_dir/{wildcards.sample}.NUMTs.collated.bam"
        """

rule align_chrM_and_NUMTs:
    input:
        chrM_fq1=rules.extract_chrM_reads.output.fq1,
        chrM_fq2=rules.extract_chrM_reads.output.fq2,
        numt_fq1=rules.extract_NUMTs_reads.output.fq1,
        numt_fq2=rules.extract_NUMTs_reads.output.fq2
    output:
        bam_chrM=temp(pj(chrM, '{sample}_chrM_orig.reads.bam')),
        bai_chrM=temp(pj(chrM,'{sample}_chrM_orig.reads.bai')),
        bam_shifted_chrM=temp(pj(chrM,'{sample}_chrM_shifted.reads.bam')),
        bai_shifted_chrM=temp(pj(chrM,'{sample}_chrM_shifted.reads.bai')),
        bam_NUMTs=temp(pj(chrM, 'NUMTs', '{sample}_NUMTs.realign.bam')),
        bai_NUMTs=temp(pj(chrM, 'NUMTs', '{sample}_NUMTs.realign.bai')),
        bam_shifted_NUMTs=temp(pj(chrM,'NUMTs','{sample}_NUMTs_shifted.reads.bam')),
        bai_shifted_NUMTs=temp(pj(chrM,'NUMTs','{sample}_NUMTs_shifted.reads.bai'))
    conda: CONDA_MAIN
    params:
        mt_ref=pj(ORIG_MT_fa),
        mt_ref_shift=pj(SHIFTED_MT_fa),
        threads_per_task=4,
        rg=lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}"
    log:
        chrM_log=pj(LOG,'chrM','{sample}.orig_mt_align.log'),
        numt_log=pj(LOG,'chrM','{sample}.origchrM_NUMT_align.log')
    resources:
        n="3",
        mem_mb=750
    shell:
        """
        mkdir -p $(dirname {log.chrM_log})
        mkdir -p $(dirname {log.numt_log})
        bwa mem -t 4 -R "{params.rg}" {params.mt_ref} {input.chrM_fq1} {input.chrM_fq2} | samtools sort -O bam -@ {params.threads_per_task} -o {output.bam_chrM}
        samtools index -@ {params.threads_per_task} -o {output.bai_chrM} {output.bam_chrM}
        bwa mem -t 4 -R "{params.rg}" {params.mt_ref_shift} {input.chrM_fq1} {input.chrM_fq2} | samtools sort -O bam -@ {params.threads_per_task} -o {output.bam_shifted_chrM}
        samtools index -@ {params.threads_per_task} -o {output.bai_shifted_chrM} {output.bam_shifted_chrM}

        bwa mem -t 4 -R "{params.rg}" {params.mt_ref} {input.numt_fq1} {input.numt_fq2} | samtools sort -O bam -@ {params.threads_per_task} -o {output.bam_NUMTs}
        samtools index -@ {params.threads_per_task} -o {output.bai_NUMTs} {output.bam_NUMTs}
        bwa mem -t 4 -R "{params.rg}" {params.mt_ref_shift} {input.numt_fq1} {input.numt_fq2} | samtools sort -O bam -@ {params.threads_per_task} -o {output.bam_shifted_NUMTs}
        samtools index -@ {params.threads_per_task} -o {output.bai_shifted_NUMTs} {output.bam_shifted_NUMTs}
        """

rule mutect_calls_both:
    input:
        bam_chrM=rules.align_chrM_and_NUMTs.output.bam_chrM,
        bai_chrM=rules.align_chrM_and_NUMTs.output.bai_chrM,
        bam_shifted_chrM=rules.align_chrM_and_NUMTs.output.bam_shifted_chrM,
        bai_shifted_chrM=rules.align_chrM_and_NUMTs.output.bai_shifted_chrM,
        bam_NUMTs=rules.align_chrM_and_NUMTs.output.bam_NUMTs,
        bai_NUMTs=rules.align_chrM_and_NUMTs.output.bai_NUMTs,
        bam_shifted_NUMTs=rules.align_chrM_and_NUMTs.output.bam_shifted_NUMTs,
        bai_shifted_NUMTs=rules.align_chrM_and_NUMTs.output.bai_shifted_NUMTs
    output:
        vcf_chrM=ensure(temp(pj(chrM, 'variants', '{sample}.chrM_orig.vcf.gz')), non_empty=True),
        tbi_chrM=ensure(temp(pj(chrM, 'variants', '{sample}.chrM_orig.vcf.gz.tbi')), non_empty=True),
        stat_chrM=ensure(temp(pj(chrM, 'variants', '{sample}.chrM_orig.vcf.gz.stats')), non_empty=True),
        vcf_shift_chrM=ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted.vcf.gz')), non_empty=True),
        tbi_shift_chrM=ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted.vcf.gz.tbi')), non_empty=True),
        stat_shift_chrM=ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted.vcf.gz.stats')), non_empty=True),
        vcf_shift_back_chrM=ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted_backshifted.vcf.gz')), non_empty=True),
        tbi_shift_back_chrM=ensure(temp(pj(chrM,'variants','{sample}.chrM_shifted_backshifted.vcf.gz.tbi')), non_empty=True),
        vcf_NUMT=ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_orig.vcf.gz')), non_empty=True),
        tbi_NUMT=ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_orig.vcf.gz.tbi')), non_empty=True),
        stat_NUMT=ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_orig.vcf.gz.stats')), non_empty=True),
        vcf_shift_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted.vcf.gz')), non_empty=True),
        tbi_shift_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted.vcf.gz.tbi')), non_empty=True),
        stat_shift_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted.vcf.gz.stats')), non_empty=True),
        vcf_shift_back_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted_backshifted.vcf.gz')), non_empty=True),
        tbi_shift_back_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_shifted_backshifted.vcf.gz.tbi')), non_empty=True)
    conda: CONDA_VCF
    params:
        mt_ref=pj(ORIG_MT_fa),
        mt_ref_shift=pj(SHIFTED_MT_fa),
        chain=pj(MT_CHAIN),
        variants_dir_chrM=pj(chrM, 'variants'),
        variants_dir_NUMT=pj(chrM, 'variants', 'NUMTs')
    resources:
        n=2,
        mem_mb=1500
    shell:
        """
        mkdir -p {params.variants_dir_chrM}
        mkdir -p {params.variants_dir_NUMT}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" Mutect2 -R {params.mt_ref} -L chrM:4142-12425 --mitochondria-mode -I {input.bam_chrM} -O {output.vcf_chrM}
        tabix -f -p vcf {output.vcf_chrM}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" Mutect2 -R {params.mt_ref_shift} -L chrM:4142-12426 --mitochondria-mode -I {input.bam_shifted_chrM} -O {output.vcf_shift_chrM}
        tabix -f -p vcf {output.vcf_shift_chrM}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" LiftoverVcf -I {output.vcf_shift_chrM} -O {output.vcf_shift_back_chrM} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null
        tabix -f -p vcf {output.vcf_shift_back_chrM}

        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" Mutect2 -R {params.mt_ref} -L chrM:4142-12425 --mitochondria-mode -I {input.bam_NUMTs} -O {output.vcf_NUMT}
        tabix -f -p vcf {output.vcf_NUMT}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" Mutect2 -R {params.mt_ref_shift} -L chrM:4142-12426 --mitochondria-mode -I {input.bam_shifted_NUMTs} -O {output.vcf_shift_NUMT}
        tabix -f -p vcf {output.vcf_shift_NUMT}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" LiftoverVcf -I {output.vcf_shift_NUMT} -O {output.vcf_shift_back_NUMT} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null
        tabix -f -p vcf {output.vcf_shift_back_NUMT}
        """

rule merge_and_filter_both:
    input:
        o_vcf_chrM=rules.mutect_calls_both.output.vcf_chrM,
        o_tbi_chrM=rules.mutect_calls_both.output.tbi_chrM,
        sb_vcf_chrM=rules.mutect_calls_both.output.vcf_shift_back_chrM,
        sb_tbi_chrM=rules.mutect_calls_both.output.tbi_shift_back_chrM,
        orig_chrM=rules.mutect_calls_both.output.stat_chrM,
        shift_chrM=rules.mutect_calls_both.output.stat_shift_chrM,
        o_vcf_NUMT=rules.mutect_calls_both.output.vcf_NUMT,
        o_tbi_NUMT=rules.mutect_calls_both.output.tbi_NUMT,
        sb_vcf_NUMT=rules.mutect_calls_both.output.vcf_shift_back_NUMT,
        sb_tbi_NUMT=rules.mutect_calls_both.output.tbi_shift_back_NUMT,
        orig_NUMT=rules.mutect_calls_both.output.stat_NUMT,
        shift_NUMT=rules.mutect_calls_both.output.stat_shift_NUMT
    output:
        merged_vcf_chrM=ensure(temp(pj(chrM, 'variants', '{sample}.chrM_merged.vcf.gz')), non_empty=True),
        merged_tbi_chrM=ensure(temp(pj(chrM, 'variants', '{sample}.chrM_merged.vcf.gz.tbi')), non_empty=True),
        merged_stat_chrM=ensure(temp(pj(chrM,'variants','{sample}.chrM_merged.vcf.gz.stats')), non_empty=True),
        filtred_vcf_chrM=ensure(temp(pj(chrM,'variants','{sample}.chrM_filtred.vcf.gz')), non_empty=True),
        filtred_tbi_chrM=ensure(temp(pj(chrM,'variants','{sample}.chrM_filtred.vcf.gz.tbi')), non_empty=True),
        merged_vcf_NUMT=ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_merged.vcf.gz')), non_empty=True),
        merged_tbi_NUMT=ensure(temp(pj(chrM, 'variants', 'NUMTs', '{sample}.chrM_NUMT_merged.vcf.gz.tbi')), non_empty=True),
        merged_stat_NUMT=ensure(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMT_merged.vcf.gz.stats')),
        filtred_vcf_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMTs_filtred.vcf.gz')), non_empty=True),
        filtred_tbi_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','{sample}.chrM_NUMTs_filtred.vcf.gz.tbi')), non_empty=True)
    conda: CONDA_VCF
    params:
        mt_ref=pj(ORIG_MT_fa)
    resources:
        n=4,
        mem_mb=500
    shell:
        """
        gatk --java-options "-XX:ActiveProcessorCount=4 -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4" MergeMutectStats --stats {input.orig_chrM} --stats {input.shift_chrM} -O {output.merged_stat_chrM}
        gatk --java-options "-XX:ActiveProcessorCount=4 -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4" MergeVcfs -I {input.sb_vcf_chrM} -I {input.o_vcf_chrM} -O {output.merged_vcf_chrM}
        tabix -f -p vcf {output.merged_vcf_chrM}
        gatk --java-options "-XX:ActiveProcessorCount=4 -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4" FilterMutectCalls -OVI true -V {output.merged_vcf_chrM} -R {params.mt_ref} --mitochondria-mode True -O {output.filtred_vcf_chrM}
        tabix -f -p vcf {output.filtred_vcf_chrM}

        gatk --java-options "-XX:ActiveProcessorCount=4 -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4" MergeMutectStats --stats {input.orig_NUMT} --stats {input.shift_NUMT} -O {output.merged_stat_NUMT}
        gatk --java-options "-XX:ActiveProcessorCount=4 -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4" MergeVcfs -I {input.sb_vcf_NUMT} -I {input.o_vcf_NUMT} -O {output.merged_vcf_NUMT}
        tabix -f -p vcf {output.merged_vcf_NUMT}
        gatk --java-options "-XX:ActiveProcessorCount=4 -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4" FilterMutectCalls -OVI true -V {output.merged_vcf_NUMT} -R {params.mt_ref} --mitochondria-mode True -O {output.filtred_vcf_NUMT}
        tabix -f -p vcf {output.filtred_vcf_NUMT}
        """

rule mutect_bp_resolution_both:
    input:
        bam_chrM=rules.align_chrM_and_NUMTs.output.bam_chrM,
        bai_chrM=rules.align_chrM_and_NUMTs.output.bai_chrM,
        bam_shifted_chrM=rules.align_chrM_and_NUMTs.output.bam_shifted_chrM,
        bai_shifted_chrM=rules.align_chrM_and_NUMTs.output.bai_shifted_chrM,
        anno_chrM=rules.merge_and_filter_both.output.filtred_vcf_chrM,
        anno_tbi_chrM=rules.merge_and_filter_both.output.filtred_tbi_chrM,
        bam_NUMT=rules.align_chrM_and_NUMTs.output.bam_NUMTs,
        bai_NUMT=rules.align_chrM_and_NUMTs.output.bai_NUMTs,
        bam_shifted_NUMT=rules.align_chrM_and_NUMTs.output.bam_shifted_NUMTs,
        bai_shifted_NUMT=rules.align_chrM_and_NUMTs.output.bai_shifted_NUMTs,
        anno_NUMT=rules.merge_and_filter_both.output.filtred_vcf_NUMT,
        anno_tbi_NUMT=rules.merge_and_filter_both.output.filtred_tbi_NUMT
    output:
        vcf_BP_chrM=ensure(temp(pj(chrM, 'variants', 'gvcf', '{sample}.chrM_orig_BP.g.vcf.gz')), non_empty=True),
        tbi_BP_chrM=ensure(temp(pj(chrM, 'variants', 'gvcf', '{sample}.chrM_orig_BP.g.vcf.gz.tbi')), non_empty=True),
        stat_BP_chrM=ensure(temp(pj(chrM, 'variants', 'gvcf', '{sample}.chrM_orig_BP.g.vcf.gz.stats')), non_empty=True),
        vcf_shift_BP_chrM=ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_BP.g.vcf.gz')), non_empty=True),
        tbi_shift_BP_chrM=ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_BP.g.vcf.gz.tbi')), non_empty=True),
        stat_shift_BP_chrM=ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_BP.g.vcf.gz.stats')), non_empty=True),
        vcf_shift_back_BP_chrM=ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_backshifted_BP.g.vcf.gz')), non_empty=True),
        tbi_shift_back_BP_chrM=ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_shifted_backshifted_BP.g.vcf.gz.tbi')), non_empty=True),
        merged_vcf_BP_chrM=ensure(temp(pj(chrM,'variants', 'gvcf','{sample}.chrM_merged_BP.g.vcf.gz')), non_empty=True),
        merged_tbi_BP_chrM=ensure(temp(pj(chrM,'variants','gvcf','{sample}.chrM_merged_BP.g.vcf.gz.tbi')), non_empty=True),
        merged_vcf_BP_with_anno_chrM=ensure(temp(pj(chrM,'variants','gvcf','{sample}.chrM_merged_BP_annotated.g.vcf.gz')), non_empty=True),
        merged_vcf_BP_with_anno_tbi_chrM=ensure(temp(pj(chrM,'variants','gvcf','{sample}.chrM_merged_BP_annotated.g.vcf.gz.tbi')), non_empty=True),
        merged_vcf_BP_norm_chrM=ensure(temp(pj(chrM,'variants','gvcf','{sample}.chrM_merged_BP_norm.g.vcf.gz')), non_empty=True),
        vcf_BP_NUMT=ensure(temp(pj(chrM, 'variants', 'NUMTs', 'gVCF', '{sample}.chrM_NUMT_orig_BP_res.g.vcf.gz')), non_empty=True),
        tbi_BP_NUMT=ensure(temp(pj(chrM, 'variants', 'NUMTs', 'gVCF', '{sample}.chrM_NUMT_orig_BP_res.g.vcf.gz.tbi')), non_empty=True),
        stat_BP_NUMT=ensure(temp(pj(chrM, 'variants', 'NUMTs', 'gVCF', '{sample}.chrM_NUMT_orig_BP_res.g.vcf.gz.stats')), non_empty=True),
        vcf_shift_BP_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_shifted_BP_res.g.vcf.gz')), non_empty=True),
        tbi_shift_BP_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_shifted_BP_res.g.vcf.gz.tbi')), non_empty=True),
        stat_shift_BP_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_shifted_BP_res.g.vcf.gz.stats')), non_empty=True),
        vcf_shift_back_BP_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_shifted_backshifted_BP_res.g.vcf.gz')), non_empty=True),
        tbi_shift_back_BP_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_shifted_backshifted_BP_res.g.vcf.gz.tbi')), non_empty=True),
        merged_vcf_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_merged.g.vcf.gz')), non_empty=True),
        merged_tbi_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_merged.g.vcf.gz.tbi')), non_empty=True),
        merged_vcf_with_anno_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_merged_with_anno.g.vcf.gz')), non_empty=True),
        merged_vcf_with_anno_tbi_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_merged_with_anno.g.vcf.gz.tbi')), non_empty=True),
        merged_vcf_norm_NUMT=ensure(temp(pj(chrM,'variants','NUMTs','gVCF','{sample}.chrM_NUMT_merged_norm.g.vcf.gz')), non_empty=True)
    conda: CONDA_VCF
    params:
        mt_ref=pj(ORIG_MT_fa),
        mt_ref_shift=pj(SHIFTED_MT_fa),
        chain=pj(MT_CHAIN),
        variants_dir_chrM=pj(chrM, 'variants'),
        variants_dir_NUMT=pj(chrM, 'variants', 'NUMTs')
    resources:
        n=1,
        mem_mb=5000
    shell:
        """
        mkdir -p {params.variants_dir_chrM}
        mkdir -p {params.variants_dir_NUMT}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" Mutect2 -ERC BP_RESOLUTION -R {params.mt_ref} -L chrM:4142-12425 --mitochondria-mode -I {input.bam_chrM} -O {output.vcf_BP_chrM}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" Mutect2 -ERC BP_RESOLUTION -R {params.mt_ref_shift} -L chrM:4142-12426 --mitochondria-mode -I {input.bam_shifted_chrM} -O {output.vcf_shift_BP_chrM}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" LiftoverVcf -I {output.vcf_shift_BP_chrM} -O {output.vcf_shift_back_BP_chrM} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" MergeVcfs -I {output.vcf_shift_back_BP_chrM} -I {output.vcf_BP_chrM} -O {output.merged_vcf_BP_chrM}
        bcftools norm -d exact -o {output.merged_vcf_BP_norm_chrM} -O z {output.merged_vcf_BP_chrM}
        tabix {output.merged_vcf_BP_norm_chrM}
        bcftools annotate -a {input.anno_chrM} -c FILTER -O z -o {output.merged_vcf_BP_with_anno_chrM} {output.merged_vcf_BP_norm_chrM}
        tabix {output.merged_vcf_BP_with_anno_chrM}

        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" Mutect2 -ERC BP_RESOLUTION -R {params.mt_ref} -L chrM:4142-12425 --mitochondria-mode -I {input.bam_NUMT} -O {output.vcf_BP_NUMT}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" Mutect2 -ERC BP_RESOLUTION -R {params.mt_ref_shift} -L chrM:4142-12426 --mitochondria-mode -I {input.bam_shifted_NUMT} -O {output.vcf_shift_BP_NUMT}
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" LiftoverVcf -I {output.vcf_shift_BP_NUMT} -O {output.vcf_shift_back_BP_NUMT} -C {params.chain} -R {params.mt_ref} --REJECT /dev/null
        gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" MergeVcfs -I {output.vcf_shift_back_BP_NUMT} -I {output.vcf_BP_NUMT} -O {output.merged_vcf_NUMT}
        bcftools norm -d exact -o {output.merged_vcf_norm_NUMT} -O z {output.merged_vcf_NUMT}
        tabix {output.merged_vcf_norm_NUMT}
        bcftools annotate -a {input.anno_NUMT} -c FILTER -O z -o {output.merged_vcf_with_anno_NUMT} {output.merged_vcf_norm_NUMT}
        tabix {output.merged_vcf_with_anno_NUMT}
        """

rule estimate_mtdna_copy_number_wes:
    input:
        cov_file="stats/cov/{sample}.regions.bed.gz"
    output:
        cn_file=ensure(pj(chrM, "stats", "{sample}.mtDNA_CN.txt"), non_empty=True)
    resources:
        n=1,
        mem_mb=100
    shell:
        r"""
        zcat {input.cov_file} | \
        awk 'BEGIN {{ total_len=0; total_cov=0; mt_len=0; mt_cov=0; }}
             $1 ~ /^chr/ {{
                 len=$3-$2;
                 if ($1 == "chrM") {{
                     mt_len += len;
                     mt_cov += $4 * len;
                 }} else if ($1 ~ /^chr[0-9XY]+$/) {{
                     total_len += len;
                     total_cov += $4 * len;
                 }}
             }}
             END {{
                 if (total_len > 0 && mt_len > 0) {{
                     nuc_cov = total_cov / total_len;
                     mt_mean_cov = mt_cov / mt_len;
                     if (nuc_cov > 0) {{
                         mtdna_cn = (mt_mean_cov / nuc_cov) * 2;
                         print "{wildcards.sample}\t"mtdna_cn;
                     }} else {{
                         print "{wildcards.sample}\tNA";
                     }}
                 }} else {{
                     print "{wildcards.sample}\tNA";
                 }}
             }}' > {output.cn_file}
        """