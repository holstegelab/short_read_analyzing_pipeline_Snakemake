"""Aggregate an already completed cohort without exposing upstream producers.

This recovery entrypoint intentionally defines only leaf aggregation rules.
Missing per-sample inputs therefore fail validation inside the job instead of
causing Snakemake to reconstruct alignment, BAM, Kraken, or variant-calling
work.
"""

from pathlib import Path

from common import *


AGGREGATE_SCRIPT = str(Path(workflow.basedir) / "scripts" / "aggregate_existing_stats.py")
AGGREGATE_OUTPUTS = [
    f"{samplefile}.{suffix}"
    for samplefile in SAMPLE_FILES
    for suffix in (
        "oxo_quality.tab",
        "bam_quality.tab",
        "bam_rg_quality.tab",
        "sex_chrom.tab",
        "coverage.hdf5",
        "kraken.tab",
        "deepvariant_bcftools.tab",
        "phase_quality.tab",
    )
]


rule all:
    input:
        AGGREGATE_OUTPUTS


rule aggregate_existing_bam_quality:
    output:
        ensure("{samplefile}.bam_quality.tab", non_empty=True)
    resources:
        n="1.0",
        mem_mb=10000
    shell:
        "python {AGGREGATE_SCRIPT:q} --kind bam --samplefile {wildcards.samplefile:q} --output {output:q}"


rule aggregate_existing_bam_rg_quality:
    output:
        ensure("{samplefile}.bam_rg_quality.tab", non_empty=True)
    resources:
        n="1.0",
        mem_mb=1000
    shell:
        "python {AGGREGATE_SCRIPT:q} --kind bam-rg --samplefile {wildcards.samplefile:q} --output {output:q}"


rule aggregate_existing_oxo_quality:
    output:
        ensure("{samplefile}.oxo_quality.tab", non_empty=True)
    resources:
        n="1.0",
        mem_mb=1000
    shell:
        "python {AGGREGATE_SCRIPT:q} --kind oxo --samplefile {wildcards.samplefile:q} --output {output:q}"


rule aggregate_existing_sex:
    output:
        ensure("{samplefile}.sex_chrom.tab", non_empty=True)
    resources:
        n="1.0",
        mem_mb=1000
    shell:
        "python {AGGREGATE_SCRIPT:q} --kind sex --samplefile {wildcards.samplefile:q} --output {output:q}"


rule aggregate_existing_coverage:
    input:
        bam="{samplefile}.bam_quality.tab"
    output:
        ensure("{samplefile}.coverage.hdf5", non_empty=True)
    resources:
        n="1.0",
        mem_mb=14700
    shell:
        "python {AGGREGATE_SCRIPT:q} --kind coverage --samplefile {wildcards.samplefile:q} --output {output:q}"


rule aggregate_existing_kraken:
    output:
        ensure("{samplefile}.kraken.tab", non_empty=True)
    resources:
        n="1.0",
        mem_mb=500
    shell:
        "python {AGGREGATE_SCRIPT:q} --kind kraken --samplefile {wildcards.samplefile:q} --output {output:q}"


rule aggregate_existing_deepvariant:
    output:
        ensure("{samplefile}.deepvariant_bcftools.tab", non_empty=True)
    resources:
        n="1.0",
        mem_mb=4000
    shell:
        "python {AGGREGATE_SCRIPT:q} --kind deepvariant --samplefile {wildcards.samplefile:q} --output {output:q}"


rule aggregate_existing_phase:
    output:
        ensure("{samplefile}.phase_quality.tab", non_empty=True)
    resources:
        n="1.0",
        mem_mb=1000
    shell:
        "python {AGGREGATE_SCRIPT:q} --kind phase --samplefile {wildcards.samplefile:q} --output {output:q}"
