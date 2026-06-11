import gzip
import json
import os

from common import *

current_dir = os.getcwd()

gvcf_caller = config.get("caller", "BOTH")
glnexus_filtration = config.get("glnexus_filtration", "custom")
genotype_mode = config.get("genotype_mode", "WES")


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


def remote_namespace_path(remote_dir, remote_name):
    return "/" + os.path.join(remote_dir, remote_name)


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


rule glnexus_to_irods_all:
    input:
        expand(
            pj(
                current_dir,
                "{genotype_mode}_{types_of_gl}" + dir_appendix,
                "ANNOTATED",
                "{region}.annotated.irods.copied"
            ),
            genotype_mode=[genotype_mode],
            types_of_gl=glnexus_dirs,
            region=level2_regions_diploid
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

        with open(output.samples_tsv, "w") as handle:
            handle.write("sample_id\n")
            for sample in samples:
                handle.write(f"{sample}\n")

        metadata = {
            "cohort_id": cohort_id(),
            "n_of_samples": len(samples),
            "region": wildcards.region,
            "samples_tsv_path": remote_namespace_path(remote_dir, remote_samples_name),
            "adler32": local_vcf_checksum,
        }

        with open(output.metadata_json, "w") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
            handle.write("\n")


rule copy_joint_vcf_bundle_to_irods:
    input:
        vcf=pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz"),
        tbi=rules.index_joint_vcf_for_irods.output.tbi,
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
