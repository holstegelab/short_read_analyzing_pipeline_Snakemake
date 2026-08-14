from common import *
wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    # readgroup="[\w\d_\-@]+"

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
#use rule * from Aligner

rule Encrypt_all:
    input: 
        expand("{cram}/{sample}.mapped_hg38.cram.copied",sample=sample_names, cram = CRAM)

sk = pj(RESOURCES,".c4gh/master_key_for_encryption")
pk1 = config.get("path_to_public_key_1",  pj(RESOURCES, ".c4gh/recipient1.pub"))
pk2 = config.get("path_to_public_key_2", pj(RESOURCES, ".c4gh/recipient2.pub"))

PKs = [pk1, pk2]


agh_dcache = config.get('agh_processed', pj(RESOURCES,".agh/agh_processed.conf"))
CRAM_DELIVERY = "delivery"
CRAM_DELIVERY_PACKAGE = pj(CRAM_DELIVERY, "projectmine_cram_decryption_package.zip")
CRAM_REFERENCE_DICT = os.path.splitext(REF)[0] + ".dict"


def cram_delivery_marker(wildcards):
    """Require one shared reference/tool upload for samples with a dCache target."""
    target = remote_base_for_sample(wildcards.sample)
    if parse_dcache_uri(target) is None:
        return []
    samplefile = os.path.basename(SAMPLEINFO[wildcards.sample]["samplefile"])
    return pj(CRAM_DELIVERY, f"{samplefile}.dcache_assets.copied")


rule build_cram_delivery_package:
    input:
        readme=srcdir("delivery/cram_decryption_package/README.md"),
        decrypt=srcdir("delivery/cram_decryption_package/decrypt_cram.py"),
        environment=srcdir("delivery/cram_decryption_package/environment.yml"),
        key_instructions=srcdir("delivery/cram_decryption_package/KEY_PACKAGE_INSTRUCTIONS.md"),
        short_readme=srcdir("delivery/cram_decryption_package/KORTE_README_CRAM_DECRYPTIE.md"),
        bam_revert=srcdir("scripts/bam_revert.py")
    output:
        package=CRAM_DELIVERY_PACKAGE
    params:
        builder=srcdir("scripts/build_cram_delivery_package.py")
    conda: CONDA_MAIN
    resources:
        time=600,
        mem_mb=200,
        n="0.1"
    shell:
        """
        python {params.builder} \
          --readme {input.readme} \
          --decrypt-script {input.decrypt} \
          --environment {input.environment} \
          --key-instructions {input.key_instructions} \
          --short-readme {input.short_readme} \
          --bam-revert {input.bam_revert} \
          --output {output.package}
        """


rule copy_cram_delivery_assets_to_dcache:
    input:
        package=rules.build_cram_delivery_package.output.package,
        bam_revert=srcdir("scripts/bam_revert.py"),
        fasta=REF,
        fai=REF + ".fai",
        dictionary=CRAM_REFERENCE_DICT
    output:
        copied=touch(pj(CRAM_DELIVERY, "{samplefile}.dcache_assets.copied")),
        fasta_sum=pj(CRAM_DELIVERY, "{samplefile}.reference.fa.ADLER32"),
        fai_sum=pj(CRAM_DELIVERY, "{samplefile}.reference.fa.fai.ADLER32"),
        dict_sum=pj(CRAM_DELIVERY, "{samplefile}.reference.dict.ADLER32"),
        bam_revert_sum=pj(CRAM_DELIVERY, "{samplefile}.bam_revert.py.ADLER32"),
        package_sum=pj(CRAM_DELIVERY, "{samplefile}.decryption_package.zip.ADLER32")
    resources:
        time=21600,
        mem_mb=500,
        n="0.2",
        dcache_upload_slots=1
    params:
        ada_script=srcdir(ADA)
    run:
        target = remote_base_for_samplefile(wildcards.samplefile)
        if parse_dcache_uri(target) is None:
            raise ValueError(
                f"CRAM delivery assets require a dCache target, got {target!r}"
            )

        target_reference = os.path.join(target, "reference")
        target_tools = os.path.join(target, "tools")
        transfers = [
            (input.fasta, target_reference, os.path.basename(input.fasta), output.fasta_sum),
            (input.fai, target_reference, os.path.basename(input.fai), output.fai_sum),
            (
                input.dictionary,
                target_reference,
                os.path.basename(input.dictionary),
                output.dict_sum,
            ),
            (
                input.bam_revert,
                target_tools,
                os.path.basename(input.bam_revert),
                output.bam_revert_sum,
            ),
            (
                input.package,
                target_tools,
                os.path.basename(input.package),
                output.package_sum,
            ),
        ]
        for local_path, remote_dir, remote_name, checksum_path in transfers:
            copy_with_checksum(
                str(local_path),
                remote_dir,
                remote_name,
                str(checksum_path),
                agh_dcache,
                params.ada_script,
            )


rule Encrypt_crams:
    input: pj(CRAM,"{sample}.mapped_hg38.cram")
    output: enCRAM=temp(pj(CRAM,"{sample}.mapped_hg38.cram.c4gh"))
    params:
        private_key = sk,
        public_key = expand("--recipient_pk {PKs}", PKs = PKs)
    conda: CONDA_MAIN         
    resources:
        time = get_time('Encrypt_crams'),
        mem_mb=200,
        n="0.3"
    shell:
        """
        python -m crypt4gh encrypt --sk {params.private_key}  {params.public_key} < {input} > {output}
        """

rule copy_to_dcache:
    input:
        cram=rules.Encrypt_crams.output.enCRAM,
        crai=pj(CRAM,"{sample}.mapped_hg38.cram.crai"),
        delivery=cram_delivery_marker
    resources:
        time = get_time('copy_to_dcache'),
        mem_mb=1500,
        n="0.2",
        dcache_upload_slots=1,
        # cram is uploaded here -> give back its share of the start_sample reservation
        # (see active_release_upload in common.py). The remainder frees at finished.
        active_use_remove=active_release_upload
    params:
        ada_script = srcdir(ADA) #temporarily as ada on Snellius is out of date
    output:
        copied = (pj(CRAM,"{sample}.mapped_hg38.cram.copied")),
        sum = (pj(CRAM,"{sample}.mapped_hg38.cram.ADLER32"))
    run:

        sample = SAMPLEINFO[wildcards['sample']]
        target = sample['target']
        samplefile = os.path.basename(sample['samplefile'])

        if target is None:
            target = os.path.join(sample['study'], samplefile)

        if target.endswith('/'):
            target = target[:-1]

        target_cram = os.path.join(target, "cram")

        input_cram = os.path.basename(input['cram'])
        input_crai = os.path.basename(input['crai'])
        copy_with_checksum(
            str(input.cram),
            target_cram,
            input_cram,
            str(output.sum),
            agh_dcache,
            params.ada_script,
        )

        crai_checksum = str(output.sum) + ".crai.tmp"
        try:
            copy_with_checksum(
                str(input.crai),
                target_cram,
                input_crai,
                crai_checksum,
                agh_dcache,
                params.ada_script,
            )
        finally:
            try:
                os.unlink(crai_checksum)
            except FileNotFoundError:
                pass

        shell("touch {output.copied}")
