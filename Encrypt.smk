from common import *
import zlib
import time
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
        mem_mb=2000,
        n="0.1"
    output:
        copied = pj(CRAM,"{sample}.mapped_hg38.cram.copied"),
        sum = pj(CRAM,"{sample}.mapped_hg38.cram.ADLER32")
    run:

        sample = SAMPLEINFO[wildcards['sample']]
        target = sample['target']

        if target is None:
            target = sample['study'] + '/' + sample['samplefile'] 

        if target.endswith('/'):
            target = target[:-1]

        input_cram = os.path.basename(input['cram'])
        input_crai = os.path.basename(input['crai'])

        with open(input.cram, 'rb') as f:
            data = f.read()
            ADLER32_local = zlib.adler32(data)
       
        ADLER32_remote = ''
        retry_counter = 0
        while f'{ADLER32_local:08x}' != ADLER32_remote and retry_counter <= 3:
            try:
                shell("rclone --config {agh_dcache} copy {input.cram} agh_processed:{target}/{input_cram}")
                shell("rclone --config {agh_dcache} copy {input.crai} agh_processed:{target}/{input_crai}")
                shell("{ada} --tokenfile {agh_dcache} --api https://dcacheview.grid.surfsara.nl:22880/api/v1 --checksum {target}/{input_cram} | awk '{{print$2}}' | awk -F '=' '{{print$2}}' > {output.sum}")
            except:
                shell("{ada} --tokenfile {agh_dcache} --api https://dcacheview.grid.surfsara.nl:22880/api/v1 --checksum {target}/{input_cram} | awk '{{print$2}}' | awk -F '=' '{{print$2}}' > {output.sum}")
            with open(output.sum, 'r') as sum_file:
                ADLER32_remote = sum_file.readline().rstrip('\n')
            retry_counter += 1
            if f'{ADLER32_local:08x}' != ADLER32_remote:
                shell("rclone --config {agh_dcache} delete agh_processed:{target}/{input_cram}")
                shell("rclone --config {agh_dcache} deletefile agh_processed:{target}/{input_crai}")
                time.sleep(60)
                

        if f'{ADLER32_local:08x}' != ADLER32_remote:
            raise ValueError(f"Checksums do not match for {input.cram} after 3 retries. Local: {ADLER32_local:08x}, Remote: {ADLER32_remote}")
        else:
            shell("touch {output.copied}")



