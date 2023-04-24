import pandas as pd
import read_stats
import os
import getpass

configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")

ref = os.path.join(config['RES'],config['ref'])
tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

wildcard_constraints:
    sample="[\w\d_\-@]+",
    extension='sam|bam|cram',
    filetype = 'fq|fastq',
    sex = 'male|female',
    batchnr='[\d]+',
    sex_ref = 'ref_male|ref_female'
    # readgroup="[\w\d_\-@]+"

from read_samples import *
import utils
from common import *
sex_ref = ['ref_male', 'ref_female']
sex_ref_hash = ['ref_male_hash', 'ref_female_hash']
sex_ref_str = ['ref_male_str', 'ref_female_str']
RES = config['RES']
rule Reference_preparation_all:
    input: expand('{RES}/{sex_ref_hash}', RES = RES, sex_ref_hash = sex_ref_hash),
            expand('{RES}/{sex_ref_str}', RES = RES, sex_ref_str = sex_ref_str),
    default_target: True


rule create_hash:
    input:
            fasta = os.path.join(config['RES'], config['{sex_ref}']),
            dir = os.path.join(config['RES'], config['{sex_ref}_dir']),
            bed = os.path.join(config['RES'], config['{sex_ref}_bed'])
    output: os.path.join(config['RES'], config['{sex_ref}_hash'])
    resources: mem_mb = 60000,
                n = 32
    shell: 
            """
            dragen-os --build-hash-table true --ht-reference {input.fasta} --output-directory {input.dir} --ht-mask-bed {input.bed} 
            """

rule ComposeSTRTableFile:
    input: fasta = os.path.join(config['RES'], config['{sex_ref}']),
    output: str_file = os.path.join(config['RES'], config['{sex_ref}_str'])
    resources: mem_mb = 12000,
                n = 2
    shell:
        """
        gatk ComposeSTRTableFile --java-options "-Xmx{resources.mem_mb}M" -R {input} -O {output}
        """


