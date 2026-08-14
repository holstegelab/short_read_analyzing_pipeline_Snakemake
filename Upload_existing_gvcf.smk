"""Package and upload only already-existing DeepVariant gVCFs.

This recovery entrypoint deliberately exposes no alignment, BAM, Kraken, or
variant-calling producers. Missing source gVCFs therefore stop the DAG instead
of triggering sample reprocessing.
"""

import os
from pathlib import Path

from common import *


wildcard_constraints:
    samplefile=r"[\w\d_\-@]+",
    region=r"[\w\d]+",


PACKAGE_SCRIPT = str(Path(workflow.basedir) / "scripts" / "package_existing_gvcf.py")
EXPECTED_REMOTE_BASE = config.get("expected_gvcf_remote_base")
if not EXPECTED_REMOTE_BASE:
    raise ValueError(
        "Pass --config expected_gvcf_remote_base=dcache:<remote>:/exact/path"
    )

for _samplefile in SAMPLE_FILES:
    _actual_remote_base = remote_base_for_samplefile(_samplefile)
    if _actual_remote_base != EXPECTED_REMOTE_BASE:
        raise ValueError(
            f"Refusing upload: {_samplefile} resolves to {_actual_remote_base!r}, "
            f"not approved expected_gvcf_remote_base={EXPECTED_REMOTE_BASE!r}"
        )


WGS_SAMPLEFILES = [
    samplefile
    for samplefile in SAMPLE_FILES
    if any(
        "wgs" in sinfo["sample_type"]
        for sinfo in SAMPLEFILE_TO_SAMPLES[samplefile].values()
    )
]


def existing_gvcf_inputs_wgs(wildcards):
    parent = convert_to_level1(wildcards.region)
    paths = []
    for sample, sinfo in SAMPLEFILE_TO_SAMPLES[wildcards.samplefile].items():
        if "wgs" not in sinfo["sample_type"]:
            continue
        base = pj(
            DEEPVARIANT,
            "gVCF",
            parent,
            f"{sample}.{parent}.wg.vcf.gz",
        )
        paths.extend((base, base + ".tbi"))
    if not paths:
        raise ValueError(
            f"No WGS samples for {wildcards.samplefile} region {wildcards.region}"
        )
    return paths


def existing_gvcf_inputs_wes(wildcards):
    parent = convert_to_level1(wildcards.region)
    paths = []
    for sample, sinfo in SAMPLEFILE_TO_SAMPLES[wildcards.samplefile].items():
        if "wgs" in sinfo["sample_type"]:
            base = pj(
                DEEPVARIANT,
                "gVCF",
                "exome_extract",
                parent,
                f"{sample}.{parent}.wg.vcf.gz",
            )
        else:
            source_region = convert_to_level0(parent)
            base = pj(
                DEEPVARIANT,
                "gVCF",
                source_region,
                f"{sample}.{source_region}.wg.vcf.gz",
            )
        paths.extend((base, base + ".tbi"))
    return paths


def gvcf_remote_dir(samplefile, kind):
    remote_dir = os.path.join(
        remote_base_for_samplefile(samplefile),
        "gvcf",
        "deepvariant",
        "level2",
        kind,
    )
    expected = os.path.join(
        EXPECTED_REMOTE_BASE,
        "gvcf",
        "deepvariant",
        "level2",
        kind,
    )
    if remote_dir != expected:
        raise ValueError(
            f"Refusing upload to unexpected destination {remote_dir!r}; "
            f"expected {expected!r}"
        )
    return remote_dir


UPLOAD_DONE = [
    pj(GVCF_TAR, "deepvariant_gvcf_wes_uploads.done"),
]
if WGS_SAMPLEFILES:
    UPLOAD_DONE.append(pj(GVCF_TAR, "deepvariant_gvcf_wgs_uploads.done"))


rule all:
    input:
        UPLOAD_DONE


rule package_existing_deepvariant_level2_wgs:
    input:
        gvcfs=existing_gvcf_inputs_wgs
    output:
        tar=temp(
            pj(
                GVCF_TAR,
                "deepvariant_level2_wgs",
                "{samplefile}.{region}.dv.wgs.gvcf.tar",
            )
        )
    resources:
        n="1.5",
        mem_mb=4000,
        ssd_use="required",
        ssd_gb=lambda wc, input: ssd_gb_for_inputs(
            input.gvcfs,
            factor=0.5,
            overhead_gb=2,
            minimum_gb=4,
        )
    conda:
        CONDA_MAIN
    shell:
        "python {PACKAGE_SCRIPT:q} --kind wgs --samplefile {wildcards.samplefile:q} "
        "--region {wildcards.region:q} --output {output.tar:q}"


rule package_existing_deepvariant_level2_wes:
    input:
        gvcfs=existing_gvcf_inputs_wes
    output:
        tar=temp(
            pj(
                GVCF_TAR,
                "deepvariant_level2_wes",
                "{samplefile}.{region}.dv.wes.gvcf.tar",
            )
        )
    resources:
        n="1.5",
        mem_mb=4000,
        ssd_use="required",
        ssd_gb=lambda wc, input: ssd_gb_for_inputs(
            input.gvcfs,
            factor=0.5,
            overhead_gb=2,
            minimum_gb=4,
        )
    conda:
        CONDA_MAIN
    shell:
        "python {PACKAGE_SCRIPT:q} --kind wes --samplefile {wildcards.samplefile:q} "
        "--region {wildcards.region:q} --output {output.tar:q}"


rule upload_existing_deepvariant_level2_wgs:
    input:
        tar=pj(
            GVCF_TAR,
            "deepvariant_level2_wgs",
            "{samplefile}.{region}.dv.wgs.gvcf.tar",
        )
    output:
        copied=touch(
            pj(
                GVCF_TAR,
                "deepvariant_level2_wgs",
                "{samplefile}.{region}.dv.wgs.gvcf.tar.copied",
            )
        ),
        checksum=ensure(
            pj(
                GVCF_TAR,
                "deepvariant_level2_wgs",
                "{samplefile}.{region}.dv.wgs.gvcf.tar.ADLER32",
            ),
            non_empty=True,
        )
    params:
        ada_script=srcdir(ADA)
    resources:
        n="0.1",
        mem_mb=2000,
        dcache_upload_slots=1,
        dcache_use_add=config.get("dcache_use_add", 0),
        dcache_use_remove=config.get("dcache_use_remove", 0)
    run:
        remote_dir = gvcf_remote_dir(wildcards.samplefile, "wgs")
        copy_with_checksum(
            str(input.tar),
            remote_dir,
            os.path.basename(str(input.tar)),
            str(output.checksum),
            AGH_DCACHE_CONFIG,
            params.ada_script,
        )


rule upload_existing_deepvariant_level2_wes:
    input:
        tar=pj(
            GVCF_TAR,
            "deepvariant_level2_wes",
            "{samplefile}.{region}.dv.wes.gvcf.tar",
        )
    output:
        copied=touch(
            pj(
                GVCF_TAR,
                "deepvariant_level2_wes",
                "{samplefile}.{region}.dv.wes.gvcf.tar.copied",
            )
        ),
        checksum=ensure(
            pj(
                GVCF_TAR,
                "deepvariant_level2_wes",
                "{samplefile}.{region}.dv.wes.gvcf.tar.ADLER32",
            ),
            non_empty=True,
        )
    params:
        ada_script=srcdir(ADA)
    resources:
        n="0.1",
        mem_mb=2000,
        dcache_upload_slots=1,
        dcache_use_add=config.get("dcache_use_add", 0),
        dcache_use_remove=config.get("dcache_use_remove", 0)
    run:
        remote_dir = gvcf_remote_dir(wildcards.samplefile, "wes")
        copy_with_checksum(
            str(input.tar),
            remote_dir,
            os.path.basename(str(input.tar)),
            str(output.checksum),
            AGH_DCACHE_CONFIG,
            params.ada_script,
        )


rule finish_existing_deepvariant_wgs_uploads:
    input:
        expand(
            pj(
                GVCF_TAR,
                "deepvariant_level2_wgs",
                "{samplefile}.{region}.dv.wgs.gvcf.tar.copied",
            ),
            samplefile=WGS_SAMPLEFILES,
            region=level2_regions,
        )
    output:
        done=touch(pj(GVCF_TAR, "deepvariant_gvcf_wgs_uploads.done"))


rule finish_existing_deepvariant_wes_uploads:
    input:
        expand(
            pj(
                GVCF_TAR,
                "deepvariant_level2_wes",
                "{samplefile}.{region}.dv.wes.gvcf.tar.copied",
            ),
            samplefile=SAMPLE_FILES,
            region=level2_regions,
        )
    output:
        done=touch(pj(GVCF_TAR, "deepvariant_gvcf_wes_uploads.done"))
