import csv
import gzip
import json
import os
import re
import subprocess
import tempfile
import yaml

from common import *

current_dir = os.getcwd()
workflow_dir = str(srcdir("."))

gvcf_caller = config.get("caller", "Deepvariant")
glnexus_filtration = config.get("glnexus_filtration", "custom")
genotype_mode = config.get("genotype_mode", "WGS")
_CONTAINER_METADATA_CACHE = {}


def analysis_name():
    for key in ("analysis", "analysis_name", "analysis_id", "cohort_id"):
        value = config.get(key)
        if value:
            return str(value)
    raise ValueError(
        "Joint VCF transfer requires an analysis name. "
        "Provide one via '--config analysis=<name>' or one of "
        "'analysis_name', 'analysis_id', 'cohort_id'."
    )


def cohort_id():
    return str(config.get("cohort_id", analysis_name()))


def transfer_token_path():
    base_token = config.get("agh_processed", AGH_DCACHE_CONFIG)
    return config.get(
        "genbiodora_processed",
        pj(os.path.dirname(base_token), "genbiodora_processed.conf")
    )


def transfer_remote_name():
    remote_name = config.get("genbiodora_processed_remote")
    if remote_name:
        return str(remote_name)

    token_path = transfer_token_path()
    try:
        with open(token_path, "r") as handle:
            for line in handle:
                line = line.strip()
                if line.startswith("[") and line.endswith("]"):
                    return line[1:-1]
    except OSError:
        pass

    return "aghub_acdc_transfer_token"


if gvcf_caller == "HaplotypeCaller":
    glnexus_dirs = ["GLnexus_on_Haplotypecaller"]
elif gvcf_caller == "Deepvariant":
    glnexus_dirs = ["GLnexus_on_Deepvariant"]
elif gvcf_caller == "BOTH":
    glnexus_dirs = ["GLnexus_on_Haplotypecaller", "GLnexus_on_Deepvariant"]
else:
    raise ValueError(
        "invalid option provided to 'caller'; please choose either "
        "'HaplotypeCaller', 'Deepvariant' or 'BOTH'."
    )

if glnexus_filtration == "default":
    dir_appendix = "default"
elif glnexus_filtration == "custom":
    dir_appendix = "custom"
else:
    raise ValueError(
        "Invalid option provided to 'glnexus_filtration'; please choose "
        "either 'default' or 'custom'."
    )


def remote_annotated_dir(wildcards):
    return os.path.join(
        "projects",
        analysis_name(),
        f"{wildcards.genotype_mode}_{wildcards.types_of_gl}{dir_appendix}",
        "ANNOTATED"
    )


def remote_analysis_stat_dir():
    return os.path.join(
        "projects",
        analysis_name(),
        "stat"
    )


def remote_namespace_path(remote_dir, remote_name):
    return "/" + os.path.join(remote_dir, remote_name)


def remote_namespace_dir(remote_dir):
    return "/" + remote_dir


def read_vcf_samples(vcf_path):
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                return header[9:]
    raise ValueError(f"Could not find '#CHROM' header in {vcf_path}")


def compute_adler32_hex(path):
    adler_value = 1
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            adler_value = zlib.adler32(chunk, adler_value)
    return f"{adler_value & 0xffffffff:08x}"


def _git_output(args):
    try:
        return subprocess.check_output(
            ["git"] + args,
            cwd=current_dir,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return None


def pipeline_git_metadata():
    commit = _git_output(["rev-parse", "HEAD"])
    commit_short = _git_output(["rev-parse", "--short", "HEAD"])
    branch = _git_output(["branch", "--show-current"])
    status = _git_output(["status", "--short"])
    return {
        "commit": commit,
        "commit_short": commit_short,
        "branch": branch,
        "is_dirty": bool(status) if status is not None else None,
    }


def workflow_path(*segments):
    return pj(workflow_dir, *segments)


def configured_version(*keys):
    for key in keys:
        value = config.get(key)
        if value:
            return str(value)
    return None


def workflow_yaml_value(key, yaml_name="Snakefile.paths.yaml"):
    path = workflow_path(yaml_name)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r") as handle:
            data = yaml.safe_load(handle) or {}
    except Exception:
        return None
    value = data.get(key)
    return str(value) if value else None


def conda_dependency_version(env_name, package_name):
    path = workflow_path(env_name)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r") as handle:
            data = yaml.safe_load(handle) or {}
    except Exception:
        return None

    for dependency in data.get("dependencies", []):
        if not isinstance(dependency, str):
            continue
        match = re.match(rf"^{re.escape(package_name)}=([^=]+)$", dependency)
        if match:
            return match.group(1)
    return None


def gatk_path_from_config():
    return config.get("GATK") or workflow_yaml_value("GATK") or gatk


def gatk_version_from_path(path):
    match = re.search(r"gatk-package-([^/]+?)\.jar", str(path))
    if match:
        return match.group(1)
    match = re.search(r"gatk[_-](\d+(?:\.\d+)*)", str(path), flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return None


def container_metadata(snakefile_name, image_token, config_prefix):
    cache_key = (snakefile_name, image_token, config_prefix)
    if cache_key in _CONTAINER_METADATA_CACHE:
        return _CONTAINER_METADATA_CACHE[cache_key]

    image = config.get(f"{config_prefix}_container")
    snakefile_path = workflow_path(snakefile_name)
    if image is None and os.path.exists(snakefile_path):
        pattern = re.compile(r"container:\s*['\"]([^'\"]*" + re.escape(image_token) + r"[^'\"]*)['\"]")
        with open(snakefile_path, "r") as handle:
            for line in handle:
                match = pattern.search(line)
                if match:
                    image = match.group(1)
                    break

    version = config.get(f"{config_prefix}_version")
    if version is None and image and ":" in image.rsplit("/", 1)[-1]:
        version = image.rsplit(":", 1)[-1]

    result = {
        "version": str(version) if version else None,
        "container": image,
    }
    _CONTAINER_METADATA_CACHE[cache_key] = result
    return result


def tool_versions_metadata():
    gatk_path = gatk_path_from_config()
    return {
        "dragen_os": {
            "version": configured_version("dragen_os_version", "dragmap_version") or conda_dependency_version("envs/dragenos.yaml", "dragmap"),
            "command": config.get("dragmap", dragmap),
        },
        "glnexus": container_metadata("GLnexus.smk", "glnexus", "glnexus"),
        "deepvariant": container_metadata("Deepvariant.smk", "deepvariant", "deepvariant"),
        "haplotypecaller": {
            "version": configured_version("haplotypecaller_version", "gatk_version") or gatk_version_from_path(gatk_path),
            "command": config.get("gatk", gatk),
            "config_path": gatk_path,
            "tool": "GATK HaplotypeCaller",
        },
    }


def reference_build():
    return str(config.get("reference_build", config.get("genome_build", "GRCh38")))


def aligner_metadata():
    return {
        "name": "DRAGMAP",
        "command": dragmap,
    }


def caller_metadata(wildcards):
    return {
        "assay_type": wildcards.genotype_mode,
        "per_sample_caller": "DeepVariant" if wildcards.types_of_gl == "GLnexus_on_Deepvariant" else "HaplotypeCaller",
        "joint_caller": "GLnexus",
        "joint_filtration": glnexus_filtration,
    }


def samplefile_for_sample(sample):
    return os.path.basename(SAMPLEINFO[sample]["samplefile"])


def samplefile_stat_path(samplefile, suffix):
    return pj(get_samplefile_folder(samplefile), f"{samplefile}.{suffix}")


def capture_kit_for_sample(sample):
    sample_type = str(SAMPLEINFO[sample].get("sample_type", "")).lower()
    if "wgs" in sample_type:
        return "WGS"
    return SAMPLEINFO[sample].get("capture_kit")


def sample_metadata(sample):
    return {
        "samplefile": samplefile_for_sample(sample),
        "batch": SAMPLE_TO_BATCH.get(sample),
        "capture_kit": capture_kit_for_sample(sample),
    }


def _summary_key(value, missing="unknown"):
    return str(value) if value else missing


def increment_count(counts, value, missing="unknown"):
    key = _summary_key(value, missing=missing)
    counts[key] = counts.get(key, 0) + 1


def sample_summary(samples):
    batch_counts = {}
    capture_kit_counts = {}
    samplefile_data = {}

    for sample in samples:
        metadata = sample_metadata(sample)
        samplefile = metadata["samplefile"]
        if samplefile not in samplefile_data:
            samplefile_data[samplefile] = {
                "samplefile": samplefile,
                "n_samples": 0,
                "batches": {},
                "capture_kits": {},
            }

        samplefile_data[samplefile]["n_samples"] += 1
        increment_count(samplefile_data[samplefile]["batches"], metadata["batch"], missing="not_staged")
        increment_count(samplefile_data[samplefile]["capture_kits"], metadata["capture_kit"])
        increment_count(batch_counts, metadata["batch"], missing="not_staged")
        increment_count(capture_kit_counts, metadata["capture_kit"])

    return {
        "batch_counts": batch_counts,
        "capture_kit_counts": capture_kit_counts,
        "samplefiles": [samplefile_data[key] for key in sorted(samplefile_data)],
    }


def bed_file_for_region(wildcards):
    return region_to_file(
        wildcards.region,
        wgs=wildcards.genotype_mode == "WGS",
        extension="bed"
    )


def read_bed_region_span(bed_path):
    first_interval = None
    last_interval = None
    with open(bed_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"Invalid BED line in {bed_path}: {line}")
            interval = {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
            }
            if first_interval is None:
                first_interval = interval
            last_interval = interval
    if first_interval is None:
        raise ValueError(f"No BED coordinates found in {bed_path}")
    if first_interval["chrom"] == last_interval["chrom"]:
        return {
            "chrom": first_interval["chrom"],
            "start": first_interval["start"],
            "end": last_interval["end"],
            "bed_file": bed_path,
        }
    return {
        "start": first_interval,
        "end": last_interval,
        "bed_file": bed_path,
    }


def ingest_marker_inputs(wildcards):
    return expand(
        pj(
            current_dir,
            "{genotype_mode}_{types_of_gl}" + dir_appendix,
            "ANNOTATED",
            "{region}.annotated.irods.copied"
        ),
        genotype_mode=[wildcards.genotype_mode],
        types_of_gl=[wildcards.types_of_gl],
        region=level2_regions_diploid
    )


def stat_copy_inputs(wildcards):
    vcf_path = pj(
        current_dir,
        f"{wildcards.genotype_mode}_{wildcards.types_of_gl}{dir_appendix}",
        "ANNOTATED",
        f"{wildcards.region}.annotated.vcf.gz"
    )
    if not os.path.exists(vcf_path):
        return []

    samplefiles = sorted({samplefile_for_sample(sample) for sample in read_vcf_samples(vcf_path)})
    return expand(
        pj(STAT, "{samplefile}.analysis_stat.copied"),
        samplefile=samplefiles
    )


rule glnexus_to_irods_all:
    input:
        expand(
            pj(
                current_dir,
                "{genotype_mode}_{types_of_gl}" + dir_appendix,
                "ANNOTATED",
                ".ingest.created"
            ),
            genotype_mode=[genotype_mode],
            types_of_gl=glnexus_dirs,
        )


rule index_joint_vcf_for_irods:
    input:
        vcf=pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz")
    output:
        tbi=pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz.tbi")
    conda:
        CONDA_MAIN
    shell:
        """
        tabix -fp vcf {input.vcf}
        """


rule copy_analysis_stats_to_irods:
    input:
        bam=lambda wc: samplefile_stat_path(wc.samplefile, "bam_quality.tab"),
        bam_rg=lambda wc: samplefile_stat_path(wc.samplefile, "bam_rg_quality.tab"),
        oxo=lambda wc: samplefile_stat_path(wc.samplefile, "oxo_quality.tab"),
        sex=lambda wc: samplefile_stat_path(wc.samplefile, "sex_chrom.tab"),
        cov=lambda wc: samplefile_stat_path(wc.samplefile, "coverage.hdf5"),
        kraken=lambda wc: samplefile_stat_path(wc.samplefile, "kraken.tab")
    output:
        copied=temp(pj(STAT, "{samplefile}.analysis_stat.copied")),
        checksums=pj(STAT, "{samplefile}.analysis_stat.transfer_checksums.tsv")
    params:
        ada_script=srcdir(ADA),
        token=lambda wc: transfer_token_path(),
        remote_name=lambda wc: transfer_remote_name(),
        remote_dir=lambda wc: remote_analysis_stat_dir()
    resources:
        mem_mb=2000,
        n="0.1",
        dcache_use_add=config.get("dcache_use_add", 0),
        dcache_use_remove=config.get("dcache_use_remove", 0)
    run:
        remote_dir = params.remote_dir
        token = params.token
        remote_name = params.remote_name

        files = [
            (str(input.bam), os.path.basename(str(input.bam))),
            (str(input.bam_rg), os.path.basename(str(input.bam_rg))),
            (str(input.oxo), os.path.basename(str(input.oxo))),
            (str(input.sex), os.path.basename(str(input.sex))),
            (str(input.cov), os.path.basename(str(input.cov))),
            (str(input.kraken), os.path.basename(str(input.kraken))),
        ]

        checksum_entries = []
        for local_path, remote_file_name in files:
            checksum_tmp = f"{output.checksums}.{remote_file_name}.tmp"
            copy_with_checksum(
                local_path,
                remote_dir,
                remote_file_name,
                checksum_tmp,
                token,
                params.ada_script,
                remote_profile=remote_name
            )
            with open(checksum_tmp, "r") as handle:
                checksum_entries.append(f"{remote_file_name}\t{handle.readline().strip()}\n")

        with open(output.checksums, "w") as handle:
            handle.writelines(checksum_entries)

        shell(f"touch {quote(str(output.copied))}")


rule prepare_joint_vcf_irods_metadata:
    input:
        vcf=pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz")
    output:
        samples_tsv=pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.samples.tsv"),
        metadata_json=pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz.json")
    run:
        os.makedirs(os.path.dirname(str(output.samples_tsv)), exist_ok=True)

        samples = read_vcf_samples(str(input.vcf))
        local_vcf_checksum = compute_adler32_hex(str(input.vcf))
        remote_dir = remote_annotated_dir(wildcards)
        remote_samples_name = os.path.basename(str(output.samples_tsv))
        region_coordinates = read_bed_region_span(bed_file_for_region(wildcards))
        pipeline_version = pipeline_git_metadata()
        callers = caller_metadata(wildcards)
        sample_summary_data = sample_summary(samples)
        capture_kits = sorted(sample_summary_data["capture_kit_counts"])

        with open(output.samples_tsv, "w") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_id", "samplefile", "batch", "capture_kit"])
            for sample in samples:
                metadata = sample_metadata(sample)
                writer.writerow([
                    sample,
                    metadata["samplefile"],
                    metadata["batch"] or "",
                    metadata["capture_kit"] or "",
                ])

        metadata = {
            "cohort_id": cohort_id(),
            "n_of_samples": len(samples),
            "region": region_coordinates,
            "reference_build": reference_build(),
            "reference_fasta": REF,
            "samples_tsv_path": remote_namespace_path(remote_dir, remote_samples_name),
            "ADLER32": local_vcf_checksum,
            "pipeline_version": pipeline_version,
            "assay_type": callers["assay_type"],
            "capture_kit": capture_kits[0] if len(capture_kits) == 1 else capture_kits,
            "capture_kit_counts": sample_summary_data["capture_kit_counts"],
            "per_sample_caller": callers["per_sample_caller"],
            "joint_caller": callers["joint_caller"],
            "joint_filtration": callers["joint_filtration"],
            "aligner": aligner_metadata(),
            "tool_versions": tool_versions_metadata(),
            "batch_counts": sample_summary_data["batch_counts"],
            "samplefiles": sample_summary_data["samplefiles"],
            "stat_dir_path": remote_namespace_dir(remote_analysis_stat_dir()),
        }

        with open(output.metadata_json, "w") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
            handle.write("\n")


rule copy_joint_vcf_bundle_to_irods:
    input:
        vcf=pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz"),
        tbi=rules.index_joint_vcf_for_irods.output.tbi,
        stats_markers=stat_copy_inputs,
        samples_tsv=rules.prepare_joint_vcf_irods_metadata.output.samples_tsv,
        metadata_json=rules.prepare_joint_vcf_irods_metadata.output.metadata_json
    output:
        copied=temp(pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.irods.copied")),
        checksums=pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.transfer_checksums.tsv")
    params:
        ada_script=srcdir(ADA),
        token=lambda wc: transfer_token_path(),
        remote_name=lambda wc: transfer_remote_name(),
        remote_dir=remote_annotated_dir
    resources:
        mem_mb=2500,
        n="0.2",
        dcache_use_add=config.get("dcache_use_add", 0),
        dcache_use_remove=config.get("dcache_use_remove", 0)
    run:
        remote_dir = params.remote_dir
        token = params.token
        remote_name = params.remote_name

        files = [
            (str(input.vcf), os.path.basename(str(input.vcf))),
            (str(input.tbi), os.path.basename(str(input.tbi))),
            (str(input.samples_tsv), os.path.basename(str(input.samples_tsv))),
            (str(input.metadata_json), os.path.basename(str(input.metadata_json))),
        ]

        checksum_entries = []
        for local_path, remote_file_name in files:
            checksum_tmp = f"{output.checksums}.{remote_file_name}.tmp"
            copy_with_checksum(
                local_path,
                remote_dir,
                remote_file_name,
                checksum_tmp,
                token,
                params.ada_script,
                remote_profile=remote_name
            )
            with open(checksum_tmp, "r") as handle:
                checksum_entries.append(f"{remote_file_name}\t{handle.readline().strip()}\n")

        with open(output.checksums, "w") as handle:
            handle.writelines(checksum_entries)

        shell(f"touch {quote(str(output.copied))}")


rule create_remote_ingest_marker:
    input:
        ingest_marker_inputs
    output:
        marker=pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", ".ingest.created")
    params:
        ada_script=srcdir(ADA),
        token=lambda wc: transfer_token_path(),
        remote_name=lambda wc: transfer_remote_name(),
        remote_dir=remote_annotated_dir
    resources:
        mem_mb=1000,
        n="0.1",
        dcache_use_add=config.get("dcache_use_add", 0),
        dcache_use_remove=config.get("dcache_use_remove", 0)
    run:
        remote_dir = params.remote_dir
        token = params.token
        remote_name = params.remote_name
        os.makedirs(os.path.dirname(str(output.marker)), exist_ok=True)

        with tempfile.NamedTemporaryFile(prefix="irods_ingest_", suffix=".tmp", dir="/tmp", delete=False) as handle:
            ingest_tmp = handle.name
        with tempfile.NamedTemporaryFile(prefix="irods_ingest_checksum_", suffix=".tmp", dir="/tmp", delete=False) as handle:
            checksum_tmp = handle.name

        try:
            copy_with_checksum(
                ingest_tmp,
                remote_dir,
                ".ingest",
                checksum_tmp,
                token,
                params.ada_script,
                remote_profile=remote_name
            )
        finally:
            if os.path.exists(ingest_tmp):
                os.unlink(ingest_tmp)
            if os.path.exists(checksum_tmp):
                os.unlink(checksum_tmp)

        shell(f"touch {quote(str(output.marker))}")
