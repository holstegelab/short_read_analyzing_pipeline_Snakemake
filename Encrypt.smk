from common import *
import zlib
wildcard_constraints:
    sample="[\w\d_\-@]+",
    # readgroup="[\w\d_\-@]+"

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
#use rule * from Aligner

rule Encrypt_all:
    input: 
        expand("{cram}/{sample}.mapped_hg38.cram.copied",sample=sample_names, cram = CRAM)
    default_target: True

sk = pj(RESOURCES,".c4gh/master_key_for_encryption")
pk1 = config.get("path_to_public_key_1",  pj(RESOURCES, ".c4gh/recipient1.pub"))
pk2 = config.get("path_to_public_key_2", pj(RESOURCES, ".c4gh/recipient2.pub"))

PKs = [pk1, pk2]


agh_dcache = config.get('agh_processed', pj(RESOURCES,".agh/agh_processed.conf"))


rule Encrypt_crams:
    input: pj(CRAM,"{sample}.mapped_hg38.cram")
    output: enCRAM=temp(pj(CRAM,"{sample}.mapped_hg38.cram.c4gh"))
    params:
        private_key = sk,
        public_key = expand("--recipient_pk {PKs}", PKs = PKs)
    conda: CONDA_MAIN         
    resources:
        mem_mb=200,
        n="0.3"
    shell:
        """
        crypt4gh encrypt --sk {params.private_key}  {params.public_key} < {input} > {output}
        """

rule copy_to_dcache:
    input:
        cram=rules.Encrypt_crams.output.enCRAM,
        crai=pj(CRAM,"{sample}.mapped_hg38.cram.crai")
    resources:
        mem_mb=200,
        n="0.1"
    output:
        copied = pj(CRAM,"{sample}.mapped_hg38.cram.copied"),
        sum = pj(CRAM,"{sample}.mapped_hg38.cram.ADLER32")
    run:

        sample = SAMPLEINFO[wildcards['sample']]
        target = sample['target']

        if target is None:
            target = sample['study'] + '/' + sample['samplefile'] + '/' + input.cram.split('/')[-1]
        if target.endswith('/'):
            target = target[:-1]

        shell("rclone --config {agh_dcache} copy {input.cram} agh_processed:{target}/")
        shell("rclone --config {agh_dcache} copy {input.crai} agh_processed:{target}/")
        shell("{ada} --tokenfile {agh_dcache} --api https://dcacheview.grid.surfsara.nl:22880/api/v1 --checksum {target}/$(basename {input.cram}) | awk '{{print$2}}' | awk -F '=' '{{print$2}}' > {output.sum}")
        with open(input.cram, 'rb') as f:
            data = f.read()
            ADLER32_local = zlib.adler32(data)
        with open(output.sum, 'r') as sum_file:
            ADLER32_remote = sum_file.readline().rstrip('\n')
        if f'{ADLER32_local:x}' != ADLER32_remote:
            shell("rclone --config {agh_dcache} delete agh_processed:{target}/$(basename {input.cram})")
            shell("rclone --config {agh_dcache} delete agh_processed:{target}/$(basename {input.crai})")
            shell("rm {input.cram}")
            raise ValueError(f"Checksums do not match for {input.cram}. Local: {ADLER32_local:x}, Remote: {ADLER32_remote}")
        else:
            shell("touch {output.copied}")



