from common import *

STORAGE_QUOTA_MODE = str(config.get('storage_quota_mode', 'auto')).strip().lower()
if STORAGE_QUOTA_MODE not in {'auto', 'required', 'disabled'}:
    raise ValueError("storage_quota_mode must be auto, required, or disabled")

wildcard_constraints:
    sample=r"[\w\d_\-@]+",


dcache_read_token = config.get(
    "token",
    config.get("dcache_read_token", DCACHE_READ_CONFIG),
)
dcache_read_prefix = config.get(
    "prefix",
    config.get("dcache_read_prefix", DCACHE_READ_PREFIX),
)
decryption_private_key = config.get(
    "path_to_decryption_private_key",
    DECRYPTION_PRIVATE_KEY,
)
decryption_passphrase_file = config.get(
    "path_to_decryption_passphrase_file", DECRYPTION_PASSPHRASE_FILE
)
decryption_passphrase_arg = (
    f"--password-file {decryption_passphrase_file}" if decryption_passphrase_file else ""
)


def dcache_target(*parts):
    return dcache_read_prefix.rstrip("/") + "/" + "/".join(str(part).strip("/") for part in parts)


def remote_cram_path(sample, suffix):
    return dcache_target(
        remote_base_for_sample(sample),
        "cram",
        f"{sample}.mapped_hg38.cram{suffix}",
    )


def remote_encrypted_cram(wildcards):
    return remote_cram_path(wildcards.sample, ".c4gh")


def remote_cram_index(wildcards):
    return remote_cram_path(wildcards.sample, ".crai")


rule Download_and_extract_bam_all:
    input:
        expand(pj(BAM, "{sample}.markdup.bam"), sample=sample_names),
        expand(pj(BAM, "{sample}.markdup.bam.bai"), sample=sample_names)
    default_target: True


rule download_encrypted_cram_from_tape:
    output:
        encrypted=temp(pj(CRAM, "{sample}.mapped_hg38.cram.c4gh")),
        crai_download=temp(pj(CRAM, "{sample}.mapped_hg38.cram.crai.downloaded"))
    params:
        remote_cram=remote_encrypted_cram,
        remote_crai=remote_cram_index,
        token=dcache_read_token,
        rclone=RCLONE,
        cram_dir=CRAM,
        quota_gate=srcdir('scripts/storage_quota.py'),
        quota_mode=STORAGE_QUOTA_MODE,
        lease_command=zslurm_lease_command(config)
    resources:
        mem_mb=2000,
        n="0.2",
        dcache_use_add=config.get("dcache_use_add", 0),
        dcache_use_remove=config.get("dcache_use_remove", 0)
    conda: CONDA_MAIN
    log:
        pj(LOG, "Download_and_extract_bam", "{sample}.download_cram.log")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.cram_dir}
        mkdir -p $(dirname {log})
        : > {log:q}
        python {params.quota_gate:q} --destination {output.encrypted:q} {output.crai_download:q} \
            --quota-mode {params.quota_mode:q} --cores {resources.n} --memory-mb {resources.mem_mb} \
            --lease-command {params.lease_command:q} --phase restore_download 2>> {log:q}
        {params.rclone:q} --config {params.token:q} -v copyto {params.remote_cram:q} {output.encrypted:q} >> {log:q} 2>&1
        {params.rclone:q} --config {params.token:q} -v copyto {params.remote_crai:q} {output.crai_download:q} >> {log:q} 2>&1
        """


rule decrypt_cram:
    input:
        encrypted=rules.download_encrypted_cram_from_tape.output.encrypted,
        crai_download=rules.download_encrypted_cram_from_tape.output.crai_download
    output:
        cram=temp(pj(CRAM, "{sample}.mapped_hg38.cram")),
        crai=temp(pj(CRAM, "{sample}.mapped_hg38.cram.crai"))
    params:
        script=srcdir("scripts/decrypt_c4gh.py"),
        private_key=decryption_private_key,
        passphrase_arg=decryption_passphrase_arg,
        quota_gate=srcdir('scripts/storage_quota.py'),
        quota_mode=STORAGE_QUOTA_MODE,
        lease_command=zslurm_lease_command(config)
    resources:
        mem_mb=1000,
        n="0.5"
    conda: CONDA_MAIN
    log:
        pj(LOG, "Download_and_extract_bam", "{sample}.decrypt_cram.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        : > {log:q}
        python {params.quota_gate:q} --destination {output.cram:q} {output.crai:q} \
            --quota-mode {params.quota_mode:q} --cores {resources.n} --memory-mb {resources.mem_mb} \
            --lease-command {params.lease_command:q} --phase restore_decrypt 2>> {log:q}
        python {params.script} --sk {params.private_key} {params.passphrase_arg} --replace {input.encrypted} {output.cram} >> {log:q} 2>&1
        cp {input.crai_download} {output.crai}
        """


rule cram_to_markdup_bam:
    input:
        cram=rules.decrypt_cram.output.cram,
        crai=rules.decrypt_cram.output.crai
    output:
        bam=ensure(pj(BAM, "{sample}.markdup.bam"), non_empty=True),
        bai=ensure(pj(BAM, "{sample}.markdup.bam.bai"), non_empty=True)
    params:
        ref=REF,
        bam_dir=BAM,
        quota_gate=srcdir('scripts/storage_quota.py'),
        quota_mode=STORAGE_QUOTA_MODE,
        lease_command=zslurm_lease_command(config)
    resources:
        n=8,
        mem_mb=14000
    conda: CONDA_MAIN
    log:
        pj(LOG, "Download_and_extract_bam", "{sample}.cram_to_bam.log")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.bam_dir}
        mkdir -p $(dirname {log})
        : > {log:q}
        python {params.quota_gate:q} --destination {output.bam:q} {output.bai:q} \
            --quota-mode {params.quota_mode:q} --cores {resources.n} --memory-mb {resources.mem_mb} \
            --lease-command {params.lease_command:q} --phase restore_bam 2>> {log:q}
        samtools view --threads {resources.n} --reference {params.ref} -O BAM -o {output.bam} {input.cram} 2>> {log:q}
        samtools index -@ {resources.n} {output.bam} {output.bai} 2>> {log}
        """
