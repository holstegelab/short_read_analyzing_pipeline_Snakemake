import pandas as pd
import read_stats
import os
import getpass
import utils
from shlex import quote
from common import *
onsuccess: shell("rm -fr logs/chrM/*")

wildcard_constraints:
    sample=PROCESSING_SAMPLE_PATTERN,


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
        copied=pj(chrM, "tar", "chrM_gvcfs.tar.copied"),
        checksum=pj(chrM, "tar", "chrM_gvcfs.tar.ADLER32")
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1",
        dcache_upload_slots=1,
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

rule chrm_extract_align_fused:
    """Extract chrM/NUMT reads and make four alignments in job scratch."""
    input:
        bam=pj(BAM, '{sample}.markdup.bam'),
        bai=pj(BAM, '{sample}.markdup.bam.bai')
    output:
        bam_chrM=temp(pj(chrM, '{sample}_chrM_orig.reads.bam')),
        bai_chrM=temp(pj(chrM, '{sample}_chrM_orig.reads.bai')),
        bam_shifted_chrM=temp(pj(chrM, '{sample}_chrM_shifted.reads.bam')),
        bai_shifted_chrM=temp(pj(chrM, '{sample}_chrM_shifted.reads.bai')),
        bam_NUMTs=temp(pj(chrM, 'NUMTs', '{sample}_NUMTs.realign.bam')),
        bai_NUMTs=temp(pj(chrM, 'NUMTs', '{sample}_NUMTs.realign.bai')),
        bam_shifted_NUMTs=temp(pj(chrM, 'NUMTs', '{sample}_NUMTs_shifted.reads.bam')),
        bai_shifted_NUMTs=temp(pj(chrM, 'NUMTs', '{sample}_NUMTs_shifted.reads.bai'))
    params:
        runner=srcdir('scripts/run_fused_chrm_extract_align.py'),
        numts=NUMTs,
        mt_ref=ORIG_MT_fa,
        mt_ref_shift=SHIFTED_MT_fa,
        threads_per_tool=2
    log:
        runner=pj(LOG, 'chrM', '{sample}.extract_align_fused.log'),
        io_profile=pj(LOG, 'chrM', '{sample}.extract_align_fused.io.json')
    conda: CONDA_MAIN
    priority: 21
    resources:
        time=get_time('chrm_extract_align_fused'),
        # Scheduling follows observed average CPU; samtools/bwa retain
        # two tool threads via params.threads_per_tool.
        n="2.0",
        mem_mb=2000,
        tmpdir=tmpdir,
        ssd_use="possible",
        ssd_gb=20
    shell:
        """
        python {params.runner:q} \
            --input-bam {input.bam:q} \
            --numts-bed {params.numts:q} \
            --original-reference {params.mt_ref:q} \
            --shifted-reference {params.mt_ref_shift:q} \
            --sample {wildcards.sample:q} \
            --output-bam-chrm {output.bam_chrM:q} \
            --output-bai-chrm {output.bai_chrM:q} \
            --output-bam-shifted-chrm {output.bam_shifted_chrM:q} \
            --output-bai-shifted-chrm {output.bai_shifted_chrM:q} \
            --output-bam-numts {output.bam_NUMTs:q} \
            --output-bai-numts {output.bai_NUMTs:q} \
            --output-bam-shifted-numts {output.bam_shifted_NUMTs:q} \
            --output-bai-shifted-numts {output.bai_shifted_NUMTs:q} \
            --metrics {log.io_profile:q} \
            --threads {params.threads_per_tool} \
            --memory-mb {resources.mem_mb} \
            --ssd-gb {resources.ssd_gb} \
            --shared-scratch-base {resources.tmpdir:q} \
            2> {log.runner:q}
        """


rule chrm_mutect_tail_fused:
    """Run Mutect, merge/filter, and BP resolution with scratch intermediates."""
    input:
        bam_chrM=pj(chrM, '{sample}_chrM_orig.reads.bam'),
        bai_chrM=pj(chrM, '{sample}_chrM_orig.reads.bai'),
        bam_shifted_chrM=pj(chrM, '{sample}_chrM_shifted.reads.bam'),
        bai_shifted_chrM=pj(chrM, '{sample}_chrM_shifted.reads.bai'),
        bam_NUMTs=pj(chrM, 'NUMTs', '{sample}_NUMTs.realign.bam'),
        bai_NUMTs=pj(chrM, 'NUMTs', '{sample}_NUMTs.realign.bai'),
        bam_shifted_NUMTs=pj(chrM, 'NUMTs', '{sample}_NUMTs_shifted.reads.bam'),
        bai_shifted_NUMTs=pj(chrM, 'NUMTs', '{sample}_NUMTs_shifted.reads.bai')
    output:
        chrM=ensure(pj(chrM, 'variants', 'gvcf', '{sample}.chrM_merged_BP_annotated.g.vcf.gz'), non_empty=True),
        chrM_tbi=ensure(pj(chrM, 'variants', 'gvcf', '{sample}.chrM_merged_BP_annotated.g.vcf.gz.tbi'), non_empty=True),
        numt=ensure(pj(chrM, 'variants', 'NUMTs', 'gVCF', '{sample}.chrM_NUMT_merged_with_anno.g.vcf.gz'), non_empty=True),
        numt_tbi=ensure(pj(chrM, 'variants', 'NUMTs', 'gVCF', '{sample}.chrM_NUMT_merged_with_anno.g.vcf.gz.tbi'), non_empty=True)
    params:
        runner=srcdir('scripts/run_fused_chrm_tail.py'),
        mt_ref=ORIG_MT_fa,
        mt_ref_shift=SHIFTED_MT_fa,
        chain=MT_CHAIN
    log:
        runner=pj(LOG, 'chrM', '{sample}.mutect_tail_fused.log'),
        io_profile=pj(LOG, 'chrM', '{sample}.mutect_tail_fused.io.json')
    conda: CONDA_VCF
    priority: 22
    resources:
        time=get_time('chrm_mutect_tail_fused'),
        # GATK's explicit ActiveProcessorCount settings remain unchanged.
        n="1.1",
        mem_mb=2500,
        tmpdir=tmpdir,
        ssd_use="possible",
        ssd_gb=lambda wildcards, input: ssd_gb_for_inputs(
            [input.bam_chrM, input.bam_shifted_chrM, input.bam_NUMTs, input.bam_shifted_NUMTs],
            factor=2.0,
            overhead_gb=8,
            minimum_gb=20,
        )
    shell:
        """
        python {params.runner:q} \
            --bam-chrm {input.bam_chrM:q} --bai-chrm {input.bai_chrM:q} \
            --bam-shifted-chrm {input.bam_shifted_chrM:q} --bai-shifted-chrm {input.bai_shifted_chrM:q} \
            --bam-numts {input.bam_NUMTs:q} --bai-numts {input.bai_NUMTs:q} \
            --bam-shifted-numts {input.bam_shifted_NUMTs:q} --bai-shifted-numts {input.bai_shifted_NUMTs:q} \
            --original-reference {params.mt_ref:q} \
            --shifted-reference {params.mt_ref_shift:q} \
            --chain {params.chain:q} \
            --sample {wildcards.sample:q} \
            --output-chrm-gvcf {output.chrM:q} --output-chrm-tbi {output.chrM_tbi:q} \
            --output-numt-gvcf {output.numt:q} --output-numt-tbi {output.numt_tbi:q} \
            --metrics {log.io_profile:q} \
            --memory-mb {resources.mem_mb} --ssd-gb {resources.ssd_gb} \
            --shared-scratch-base {resources.tmpdir:q} \
            2> {log.runner:q}
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
