import pandas as pd
import read_stats
import os
import getpass

configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")
gatk = config['gatk']
samtools = config['samtools']
bcftools = config['bcftools']
dragmap = config['dragmap']
verifybamid2 = config['verifybamid2']


tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())
tmpdir_alternative = os.path.join(config['tmpdir'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

wildcard_constraints:
    sample="[\w\d_\-@]+",
    # readgroup="[\w\d_\-@]+"


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

rule Encrypt_all:
    input: 
        expand("{cram}/{sample}.mapped_hg38.cram.copied",sample=sample_names, cram = config['CRAM'])
    default_target: True

sk = config.get("private_key", os.path.join(config['RES'],".c4gh/master_key_for_encryption"))
pk1 = config.get("public_key", os.path.join(config['RES'], ".c4gh/recipient1.pub"))
pk2 = config.get("public_key_2", os.path.join(config['RES'],".c4gh/recipient2.pub"))

PKs = [pk1, pk2]


agh_dcache = config.get('agh_processed', os.path.join(config['RES'],".agh/agh_processed.conf"))


rule Encrypt_crams:
    input: rules.mCRAM.output.cram
    output: enCRAM=temp(os.path.join(config['CRAM'],"{sample}.mapped_hg38.cram.c4gh"))
    params:
            private_key = sk,
            public_key = expand("--recipient_pk {PKs}", PKs = PKs)
    conda: "envs/preprocess.yaml"            
    shell:
        """
        crypt4gh encrypt --sk {params.private_key}  {params.public_key} < {input} > {output}
        """

rule copy_to_dcache:
    input:
        cram=rules.Encrypt_crams.output.enCRAM,
        crai=rules.mCRAM.output.crai
    output:
        os.path.join(config['CRAM'],"{sample}.mapped_hg38.cram.copied")
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

