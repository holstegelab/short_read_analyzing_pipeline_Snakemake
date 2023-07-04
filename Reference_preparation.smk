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
    sex_ref = 'ref_male|ref_female',
    sex_ref_hash = ''
    # readgroup="[\w\d_\-@]+"

from read_samples import *
import utils
from common import *

sex_ref = ['ref_male', 'ref_female']
sex_ref_hash = ['ref_male_hash', 'ref_female_hash']
sex_ref_str = ['ref_male_str', 'ref_female_str']
RES = config['RES']

rule Reference_preparation_all:
    input:
        expand('{RES}{sex_ref_hash}', RES=RES, sex_ref_hash=config['ref_male_hash']),
        expand('{RES}{sex_ref_str}', RES=RES, sex_ref_str=config['ref_male_str']),
        expand('{RES}{sex_ref_hash}', RES=RES, sex_ref_hash=config['ref_female_hash']),
        expand('{RES}{sex_ref_str}', RES=RES, sex_ref_str=config['ref_female_str']),
        os.path.join(config['RES'],config['SHIFTED_MT_fai']),
        os.path.join(config['RES'],config['ORIG_MT_fai']),
        os.path.join(config['RES'],config['SHIFTED_MT_dict']),
        os.path.join(config['RES'],config['ORIG_MT_dict'])
    default_target: True


rule create_fai:
    input: fasta=ancient(os.path.join(config['RES'], config['ref_male']))
    output: fai = os.path.join(config['RES'], config['ref_male_fai'])
    conda: "envs/preprocess.yaml"
    shell: "samtools faidx {input}"

rule create_fai_1:
    input: fasta=ancient(os.path.join(config['RES'], config['ref_female']))
    output: fai = os.path.join(config['RES'], config['ref_female_fai'])
    conda: "envs/preprocess.yaml"
    shell: "samtools faidx {input}"

rule create_fai_chrM_shifted:
    input: mt_ref_shift = ancient(os.path.join(config['RES'], config['SHIFTED_MT_fa']))
    output: fai = os.path.join(config['RES'], config['SHIFTED_MT_fai'])
    conda: "envs/preprocess.yaml"
    shell: "samtools faidx {input}"

rule create_dict_for_chrM_orig_reference:
    input: mt_ref_shift = ancient(os.path.join(config['RES'], config['ORIG_MT_fa']))
    output: fai = os.path.join(config['RES'], config['ORIG_MT_fai'])
    conda: "envs/preprocess.yaml"
    shell: "samtools faidx {input}"

rule create_dict:
    input: fasta=ancient(os.path.join(config['RES'], config['ref_male'])),
    output: dict = os.path.join(config['RES'], config['ref_male_dict'])
    conda: "envs/gatk.yaml"
    shell: "gatk CreateSequenceDictionary -R {input}"

rule create_dict_1:
    input: fasta= ancient(os.path.join(config['RES'], config['ref_female'])),
    output: dict =os.path.join(config['RES'], config['ref_female_dict'])
    conda: "envs/gatk.yaml"
    shell: "gatk CreateSequenceDictionary -R {input}"

rule create_dict_for_chrM_shifted_reference:
    input: mt_ref_shift = ancient(os.path.join(config['RES'], config['SHIFTED_MT_fa']))
    output: dict = os.path.join(config['RES'], config['SHIFTED_MT_dict'])
    conda: "envs/gatk.yaml"
    shell:  """
            gatk CreateSequenceDictionary -R {input}
            """

rule create_dict_for_chrM_reference:
    input: mt_ref_shift = ancient(os.path.join(config['RES'], config['ORIG_MT_fa']))
    output: dict = os.path.join(config['RES'], config['ORIG_MT_dict'])
    conda: "envs/gatk.yaml"
    shell:  """
            gatk CreateSequenceDictionary -R {input}
            """

rule create_hash:
    input:
        fasta=ancient(os.path.join(config['RES'], config['ref_male'])),
        dir=ancient(os.path.join(config['RES'], config['ref_male_dir'])),
        bed=ancient(os.path.join(config['RES'], config['ref_male_bed'])),
        dict = ancient(os.path.join(config['RES'], config['ref_male_dict'] ))
    conda: "envs/preprocess.yaml"
    output: hash = os.path.join(config['RES'], config[sex_ref_hash[0]])
    params:
        dragmap = os.path.join(config['RES'],config['SOFTWARE'],'dragen-os')
    resources:
        mem_mb=50000,
        n=64
    shell:  """
            {params.dragmap} --build-hash-table true --ht-reference {input.fasta} --output-directory {input.dir} --ht-mask-bed {input.bed} 
            """

rule ComposeSTRTableFile:
    input:
        fasta=ancient(os.path.join(config['RES'], config['ref_male'])),
        dict= ancient(os.path.join(config['RES'],config['ref_male_dict']))
    output:
        str_file=os.path.join(config['RES'], config['ref_male_str'])
    conda: "envs/preprocess.yaml"
    resources:
        mem_mb=12000,
        n=2
    shell:"""
        gatk ComposeSTRTableFile --java-options "-Xmx{resources.mem_mb}M" -R {input.fasta} -O {output.str_file}
        """

rule create_hash_1:
    input:
        fasta=ancient(os.path.join(config['RES'], config['ref_female'])),
        dir=ancient(os.path.join(config['RES'], config['ref_female_dir'])),
        bed=ancient(os.path.join(config['RES'], config['ref_female_bed'])),
        dict= ancient(os.path.join(config['RES'],config['ref_female_dict']))
    conda: "envs/preprocess.yaml"
    output: hash = os.path.join(config['RES'], config['ref_female_hash'])
    params:
        dragmap = os.path.join(config['RES'],config['SOFTWARE'],'dragen-os')
    resources:
        mem_mb=50000,
        n=64
    shell:  """
            {params.dragmap} --build-hash-table true --ht-reference {input.fasta} --output-directory {input.dir} --ht-mask-bed {input.bed} 
            """

rule ComposeSTRTableFile_1:
    input:
        fasta=ancient(os.path.join(config['RES'], config['ref_female'])),
        dict= ancient(os.path.join(config['RES'],config['ref_female_dict']))
    output:
        str_file = os.path.join(config['RES'], config['ref_female_str'])
    conda: "envs/preprocess.yaml"
    resources:
        mem_mb=12000,
        n=2
    shell:"""
        gatk ComposeSTRTableFile --java-options "-Xmx{resources.mem_mb}M" -R {input.fasta} -O {output.str_file}
        """


