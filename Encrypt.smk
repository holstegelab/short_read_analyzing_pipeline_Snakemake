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
        python -m crypt4gh encrypt --sk {params.private_key}  {params.public_key} < {input} > {output}
        """

rule copy_to_dcache:
    input:
        cram=rules.Encrypt_crams.output.enCRAM,
        crai=pj(CRAM,"{sample}.mapped_hg38.cram.crai")
    resources:
        mem_mb=1500,
        n="0.2"
    params:
        ada_script = srcdir(ADA) #temporarily as ada on Snellius is out of date
    output:
        copied = temp(pj(CRAM,"{sample}.mapped_hg38.cram.copied")),
        sum = temp(pj(CRAM,"{sample}.mapped_hg38.cram.ADLER32"))
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

        ADLER32_local = 1
        with open(input.cram, 'rb') as f:
            for chunk in iter(lambda: f.read(16 * 1024 * 1024), b''):
                ADLER32_local = zlib.adler32(chunk, ADLER32_local)
        ADLER32_local &= 0xffffffff

        import sys
        sys.stderr.write('PATHS: ' + target + ' ' +  sample['study'] + ' ' +  sample['samplefile'] + ' ' +  input_cram)
        sys.stderr.flush()
       
        ADLER32_remote = ''
        retry_counter = 0
        sys.stderr.write(f'\nCH {ADLER32_local:08x}\n')
        sys.stderr.flush()

        while f'{ADLER32_local:08x}' != ADLER32_remote and retry_counter <= 3:
            shell("rclone --config {agh_dcache} mkdir -v agh_processed:{target_cram}")            
            
            shell("rclone --config {agh_dcache} -v copyto {input.cram} agh_processed:{target_cram}/{input_cram}")
            
            shell("rclone --config {agh_dcache} -v copyto {input.crai} agh_processed:{target_cram}/{input_crai}")

            shell("{params.ada_script} --tokenfile {agh_dcache} --api https://dcacheview.grid.surfsara.nl:22880/api/v1 --checksum {target_cram}/{input_cram} | awk '{{print $2}}' | awk -F '=' '{{print $2}}' > {output.sum}")
              
            with open(output.sum, 'r') as sum_file:
                ADLER32_remote = sum_file.readline().rstrip('\n')
            retry_counter += 1

            sys.stderr.write(f'\nCOMPARE #{ADLER32_local:08x}# != #{ADLER32_remote}#\n')
            sys.stderr.write(f'{target}\n')
            sys.stderr.write(f'{input_cram}\n')
            sys.stderr.write(f'{input_crai}\n')
            sys.stderr.flush()

            if f'{ADLER32_local:08x}' != ADLER32_remote:
                shell("rclone --config {agh_dcache} -v deletefile  agh_processed:{target_cram}/{input_cram}")
                shell("rclone --config {agh_dcache} -v deletefile agh_processed:{target_cram}/{input_crai}")
                time.sleep(60)
                

        if f'{ADLER32_local:08x}' != ADLER32_remote:
            raise ValueError(f"Checksums do not match for {input.cram} after 3 retries. Local: {ADLER32_local:08x}, Remote: {ADLER32_remote}")
        else:
            shell("touch {output.copied}")



