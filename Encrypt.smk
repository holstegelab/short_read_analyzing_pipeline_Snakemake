from common import *
wildcard_constraints:
    sample="[\w\d_\-@]+",
    # readgroup="[\w\d_\-@]+"

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

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
    input: rules.mCRAM.output.cram
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
        crai=rules.mCRAM.output.crai
    resources:
        mem_mb=200,
        n="0.1"
    output:
        pj(CRAM,"{sample}.mapped_hg38.cram.copied")
    run:
        sample = SAMPLEINFO[wildcards['sample']]
        target = sample['target']

        if target is None:
            target = sample['study'] + '/' + sample['samplefile']
        if target.endswith('/'):
            target = target[:-1]

        shell("rclone --config {agh_dcache} copy {input.cram} agh_processed:{target}/")    
        shell("rclone --config {agh_dcache} copy {input.crai} agh_processed:{target}/")    
        shell("touch {output}")

