from common import *

wildcard_constraints:
    sample=PROCESSING_SAMPLE_PATTERN,
    region = r"[\w\d]+",
    # readgroup="[\w\d_\-@]+"
onsuccess: shell("rm -fr logs/Deepvariant/*")
from common import *
from shlex import quote
import tarfile
import tempfile
import shutil
from pathlib import Path


module Tools:
    snakefile: 'Tools.smk'
    config: config
use rule * from Tools
mode = config.get("computing_mode", "WES")
DEEPVARIANT_APPTAINER = config.get("deepvariant_apptainer_output", f"{DEEPVARIANT}_apptainer")
DEEPVARIANT_NATIVE_PREFIX = config.get(
    "deepvariant_native_prefix",
    os.environ.get("DEEPVARIANT_NATIVE_PREFIX", DEEPVARIANT_NATIVE_RUNTIME),
)


DEEPVARIANT_LEASE_MODE = str(
    config.get('deepvariant_lease_mode', 'required')
).strip().lower()
if DEEPVARIANT_LEASE_MODE not in {'required', 'optional', 'disabled'}:
    raise ValueError(
        "deepvariant_lease_mode must be required, optional, or disabled"
    )


def get_deepvariant_native_runner(wildcards):
    if not DEEPVARIANT_NATIVE_PREFIX:
        raise ValueError(
            "DeepVariant requires --config deepvariant_native_prefix=/path/to/runtime "
            "or the DEEPVARIANT_NATIVE_PREFIX environment variable"
        )
    prefix = Path(DEEPVARIANT_NATIVE_PREFIX).expanduser()
    runner = prefix / "bin/run_deepvariant"
    ready = prefix / ".deepvariant-native.ready"
    if not ready.is_file() or not runner.is_file() or not os.access(runner, os.X_OK):
        raise FileNotFoundError(
            f"DeepVariant native runtime is incomplete at {prefix}. "
            "Run scripts/prepare_deepvariant_native.py --prefix PATH first."
        )
    return str(runner)


def level2_parent_level1(region):
    """Return the parent level1 region for a given level2 region."""
    return convert_to_level1(region)


def level2_interval(region, wgs):
    """Return the padded interval list for the level2 region, using WGS or WES bins."""
    return region_to_file(region, wgs=wgs, extension="bed")


def deepvariant_level2_inputs_wgs(wildcards):
    """Collect per-sample WGS gVCF paths for a given level2 region and samplefile."""
    region = wildcards.region
    samplefile = wildcards.samplefile
    parent_level1 = level2_parent_level1(region)
    files = []
    for sample, sinfo in SAMPLEFILE_TO_SAMPLES[samplefile].items():
        if 'wgs' not in sinfo['sample_type']:
            continue
        base_path = pj(DEEPVARIANT, "gVCF", parent_level1, f"{sample}.{parent_level1}.wg.vcf.gz")
        files.append(base_path)
        files.append(base_path + ".tbi")
    if not files:
        raise ValueError(f"No WGS DeepVariant gVCFs found for samplefile {samplefile} in region {region}")
    return files


def deepvariant_level2_inputs_wes(wildcards):
    """Collect per-sample exome gVCF paths for a given level2 region and samplefile."""
    region = wildcards.region
    samplefile = wildcards.samplefile
    parent_level1 = level2_parent_level1(region)
    files = []
    for sample, sinfo in SAMPLEFILE_TO_SAMPLES[samplefile].items():
        if 'wgs' in sinfo['sample_type']:
            base_path = pj(DEEPVARIANT, "gVCF", "exome_extract", parent_level1, f"{sample}.{parent_level1}.wg.vcf.gz")
        else:
            base0 = convert_to_level0(parent_level1)
            base_path = pj(DEEPVARIANT, "gVCF", base0, f"{sample}.{base0}.wg.vcf.gz")
        files.append(base_path)
        files.append(base_path + ".tbi")
    return files


def deepvariant_level2_wgs_samples(wildcards):
    return [
        sample
        for sample, sinfo in SAMPLEFILE_TO_SAMPLES[wildcards.samplefile].items()
        if 'wgs' in sinfo['sample_type']
    ]


def deepvariant_level2_samples(wildcards):
    return list(SAMPLEFILE_TO_SAMPLES[wildcards.samplefile].keys())


DEEPVARIANT_LEVEL2_WORKERS = int(config.get("deepvariant_level2_workers", 8))
if DEEPVARIANT_LEVEL2_WORKERS < 1:
    raise ValueError("deepvariant_level2_workers must be at least 1")


rule extract_and_tar_deepvariant_level2_wgs:
    input:
        gvcfs=deepvariant_level2_inputs_wgs
    output:
        tar=temp(pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar"))
    params:
        samples=deepvariant_level2_wgs_samples,
        region=lambda wc: wc.region,
        samplefile=lambda wc: wc.samplefile,
        interval=lambda wc: level2_interval(wc.region, wgs=True),
        dataset="wgs",
        lease_command=srcdir("scripts/zslurm_lease_client.py")
    threads: DEEPVARIANT_LEVEL2_WORKERS
    conda:
        CONDA_VCF
    resources:
        time = get_time('extract_and_tar_deepvariant_level2_wgs'),
        mem_mb=4000,
        n=lambda wildcards, threads: str(threads)
    script:
        "scripts/extract_and_tar_deepvariant_level2.py"


rule extract_and_tar_deepvariant_level2_wes:
    input:
        gvcfs=deepvariant_level2_inputs_wes
    output:
        tar=temp(pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar"))
    params:
        samples=deepvariant_level2_samples,
        region=lambda wc: wc.region,
        samplefile=lambda wc: wc.samplefile,
        interval=lambda wc: level2_interval(wc.region, wgs=False),
        dataset="wes",
        lease_command=srcdir("scripts/zslurm_lease_client.py")
    threads: DEEPVARIANT_LEVEL2_WORKERS
    conda:
        CONDA_VCF
    resources:
        time = get_time('extract_and_tar_deepvariant_level2_wes'),
        mem_mb=4000,
        n=lambda wildcards, threads: str(threads)
    script:
        "scripts/extract_and_tar_deepvariant_level2.py"


rule DeepVariant_all:
    input:
        expand("{dv}/{sample}.done",sample=sample_names, dv = DEEPVARIANT),


WGS_SAMPLEFILES = [
    samplefile
    for samplefile in SAMPLE_FILES
    if any('wgs' in sinfo['sample_type'] for sinfo in SAMPLEFILE_TO_SAMPLES[samplefile].values())
]


def deepvariant_region_inputs(wildcards):
    region = wildcards.region
    chunk = convert_to_level1(region)
    files = []
    for sample, sinfo in SAMPLEFILE_TO_SAMPLES[wildcards.samplefile].items():
        if 'wgs' not in sinfo['sample_type']:
            continue
        base = pj(DEEPVARIANT, "gVCF", chunk, f"{sample}.{chunk}.wg.vcf.gz")
        files.append(base)
        files.append(base + ".tbi")
    if not files:
        raise ValueError(f"No DeepVariant gVCFs found for samplefile {wildcards.samplefile} in region {region}")
    return files


rule copy_deepvariant_wgs_region_to_dcache:
    input:
        tar=pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar")
    output:
        copied=pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar.copied"),
        checksum=pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar.ADLER32")
    params:
        ada_script=srcdir(ADA)
    resources:
        time = get_time('copy_deepvariant_wgs_region_to_dcache'),
        mem_mb=2000,
        n="0.1",
        dcache_upload_slots=1,
        dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        remote_dir = os.path.join(remote_base_for_samplefile(wildcards.samplefile), "gvcf", "deepvariant", "level2", "wgs")
        remote_name = os.path.basename(input.tar)
        copy_with_checksum(str(input.tar), remote_dir, remote_name, str(output.checksum), AGH_DCACHE_CONFIG, params.ada_script)
        shell(f"touch {quote(str(output.copied))}")


rule copy_deepvariant_wes_region_to_dcache:
    input:
        tar=pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar")
    output:
        copied=pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar.copied"),
        checksum=pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar.ADLER32")
    params:
        ada_script=srcdir(ADA)
    resources:
        time = get_time('copy_deepvariant_wes_region_to_dcache'),
        mem_mb=2000,
        n="0.1",
        dcache_upload_slots=1,
        dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        remote_dir = os.path.join(remote_base_for_samplefile(wildcards.samplefile), "gvcf", "deepvariant", "level2", "wes")
        remote_name = os.path.basename(input.tar)
        copy_with_checksum(str(input.tar), remote_dir, remote_name, str(output.checksum), AGH_DCACHE_CONFIG, params.ada_script)
        shell(f"touch {quote(str(output.copied))}")


rule deepvariant_tar_wgs_all:
    input:
        expand(pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar.copied"), samplefile=WGS_SAMPLEFILES, region=level2_regions)
    output:
        done=touch(pj(GVCF_TAR, "deepvariant_gvcf_wgs_uploads.done"))


rule deepvariant_tar_wes_all:
    input:
        expand(pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar.copied"), samplefile=SAMPLE_FILES, region=level2_regions)
    output:
        done=touch(pj(GVCF_TAR, "deepvariant_gvcf_wes_uploads.done"))


def get_deepvariant_files(wildcards):#{{{
    sample = wildcards['sample']
    if 'wgs' in SAMPLEINFO[sample]['sample_type']:
        return [pj(DEEPVARIANT,  'gVCF', 'exome_extract', region, f'{sample}.{region}.wg.vcf.gz') for region in level1_regions]
    else:
        return [pj(DEEPVARIANT,  'gVCF', region, f'{sample}.{region}.wg.vcf.gz') for region in level0_regions]
#}}}


rule deepvariant_sample_done:
    input:
        get_deepvariant_files
    output:
        done=touch(pj(DEEPVARIANT, "{sample}.done"))
    resources:
        mem_mb = 100,
        n = "1.0"

def get_sequencing_mode(wildcards):#{{{
    return "WGS" if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else "WES"#}}}

def get_mem_mb_deepvariant(wildcards, attempt):#{{{
    # Size the normal attempt at the observed median process-tree PSS.  The
    # ZSlurm reservation is a packing estimate, not a hard memory limit.
    res = 8700
    return (attempt - 1) * 0.5 * res + res#}}}


def region_to_bed_file(wildcards):#{{{
    """Converts a region to a bed file location (see common.py and Tools.smk)"""
    sample = wildcards['sample']
    region = wildcards['region']
    return region_to_file(region, wgs='wgs' in SAMPLEINFO[sample]['sample_type'], extension='bed')#}}}

def region_to_bed_file_wgs(wildcards):#{{{
    region = wildcards['region']
    return region_to_file(region, wgs=True, extension='bed')#}}}

include: "Deepvariant_apptainer.smk"


rule deepvariant_phasing_fused:
    """Call, phase, merge, and extract one regional DeepVariant gVCF."""
    input:
        bed=region_to_bed_file,
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        validated_sex=pj(KMER,"{sample}.result.yaml")
    output:
        vcf=temp(pj(DEEPVARIANT, "VCF/{region}/{sample}.{region}.w.vcf.gz")),
        vcf_tbi=temp(pj(DEEPVARIANT, "VCF/{region}/{sample}.{region}.w.vcf.gz.tbi")),
        wstats=pj(STAT, "whatshap_dvphasing/{sample}.{region}.stats"),
        mwstats=pj(STAT, "whatshap_dvphasing/{sample}.{region}.merge_stats"),
        bcftools_stats=temp(pj(STAT, "deepvariant_bcftools/{sample}.{region}.bcftools_stats.txt")),
        bcftools_summary=ensure(pj(STAT, "deepvariant_bcftools/{sample}.{region}.summary.tsv"), non_empty=True),
        tmp_gvcf=temp(pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf")),
        gvcf=pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz"),
        gvcf_tbi=pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz.tbi"),
        gvcf_exome=ensure(pj(DEEPVARIANT, "gVCF/exome_extract/{region}/{sample}.{region}.wg.vcf.gz"), non_empty=True),
        gvcf_exome_tbi=ensure(pj(DEEPVARIANT, "gVCF/exome_extract/{region}/{sample}.{region}.wg.vcf.gz.tbi"), non_empty=True)
    log:
        runner=pj(LOG, "Deepvariant", "{sample}.{region}.deepvariant_phasing_fused.log"),
        io_profile=pj(LOG, "Deepvariant", "{sample}.{region}.deepvariant_phasing_fused.io.json")
    params:
        runner=srcdir('scripts/run_fused_deepvariant_phasing.py'),
        deepvariant_runner=get_deepvariant_native_runner,
        mode=get_sequencing_mode,
        haploid_contigs=lambda wc: 'chrX,chrX_KI270880v1_alt,chrX_KI270881v1_alt,chrX_KI270913v1_alt,chrY,chrY_KI270740v1_random' if wc.region.endswith('H') else 'chrNONE',
        ploidy=lambda wc: 1 if wc.region.endswith('H') else 2,
        skipsex=lambda wc, input: int(get_validated_sex_file(input) == 'female' and wc.region.startswith('Y')),
        interval_bed=lambda wc: region_to_file(region=wc.region, extension='bed', padding=True),
        merge_script=srcdir(MERGEPHASEDIRECT),
        stats_parser=srcdir('scripts/deepvariant_bcftools_stats_parser.py'),
        lease_mode=DEEPVARIANT_LEASE_MODE,
        lease_command=zslurm_lease_command(config)
    conda: CONDA_VCF
    resources:
        # Reserve observed average CPU; DeepVariant still runs 8 shards.
        n="7.5",
        nshards=8,
        mem_mb=get_mem_mb_deepvariant,
        attempt=lambda wildcards, attempt: attempt,
        time=get_time('deepvariant_phasing_fused'),
        tmpdir=tmpdir,
        ssd_use="possible",
        ssd_gb=16
    shell:
        """
        set +e
        printf '[deepvariant_phasing_fused] snakemake_attempt=%s\\n' \
            {resources.attempt:q} > {log.runner:q}
        rm -f -- {log.io_profile:q}
        python {params.runner:q} \
            --sample {wildcards.sample:q} \
            --region {wildcards.region:q} \
            --bed {input.bed:q} \
            --bam {input.bam:q} \
            --bai {input.bai:q} \
            --validated-sex {input.validated_sex:q} \
            --reference {REF:q} \
            --deepvariant-reference {REF_MALE:q} \
            --deepvariant-runner {params.deepvariant_runner:q} \
            --model-type {params.mode:q} \
            --haploid-contigs {params.haploid_contigs:q} \
            --ploidy {params.ploidy} \
            --skip-sex {params.skipsex} \
            --interval-bed {params.interval_bed:q} \
            --capture-auto-bed {INTERSECT_CAPTURE_KIT_AUTO_BED:q} \
            --capture-x-bed {INTERSECT_CAPTURE_KIT_X_BED:q} \
            --capture-y-bed {INTERSECT_CAPTURE_KIT_Y_BED:q} \
            --merge-script {params.merge_script:q} \
            --stats-parser {params.stats_parser:q} \
            --output-vcf {output.vcf:q} \
            --output-vcf-tbi {output.vcf_tbi:q} \
            --output-wstats {output.wstats:q} \
            --output-merge-stats {output.mwstats:q} \
            --output-bcftools-stats {output.bcftools_stats:q} \
            --output-bcftools-summary {output.bcftools_summary:q} \
            --output-tmp-gvcf {output.tmp_gvcf:q} \
            --output-gvcf {output.gvcf:q} \
            --output-gvcf-tbi {output.gvcf_tbi:q} \
            --output-exome-gvcf {output.gvcf_exome:q} \
            --output-exome-gvcf-tbi {output.gvcf_exome_tbi:q} \
            --metrics {log.io_profile:q} \
            --num-shards {resources.nshards} \
            --initial-cores {resources.n} \
            --initial-memory-mb {resources.mem_mb} \
            --low-cores 1 \
            --low-memory-mb 4000 \
            --attempt {resources.attempt} \
            --lease-mode {params.lease_mode:q} \
            --lease-command {params.lease_command:q} \
            --ssd-gb {resources.ssd_gb} \
            --shared-scratch-base {resources.tmpdir:q} \
            2>> {log.runner:q}
        status=$?
        set -e
        if [ "$status" -ne 0 ]; then
            cp -- {log.runner:q} "{log.runner}.attempt-{resources.attempt}.failed.log" || true
            if [ -s {log.io_profile:q} ]; then
                cp -- {log.io_profile:q} "{log.io_profile}.attempt-{resources.attempt}.failed.json" || true
            fi
        fi
        exit "$status"
        """
