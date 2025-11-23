import pandas as pd
import read_stats
import os
import getpass

wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    extension=r'sam|bam|cram',
    filetype = r'fq|fastq',
    sex = r'male|female',
    batchnr=r'[\d]+',
    sex_ref = r'ref_male|ref_female',
    sex_ref_hash = ''
    # readgroup="[\w\d_\-@]+"

from read_samples import *
import utils
from common import *
import constants

rule Reference_preparation_all:
    input:
        expand('{RES}/{sex_ref_hash}', RES=RESOURCES, sex_ref_hash=REF_MALE_HASH),
        expand('{RES}/{sex_ref_str}', RES=RESOURCES, sex_ref_str=REF_MALE_STR),
        expand('{RES}/{sex_ref_hash}', RES=RESOURCES, sex_ref_hash=REF_FEMALE_HASH),
        expand('{RES}/{sex_ref_str}', RES=RESOURCES, sex_ref_str=REF_FEMALE_STR),
        SHIFTED_MT_fai,
        ORIG_MT_fai,
        SHIFTED_MT_dict,
        ORIG_MT_dict,
        ORIG_MT_fa + ".amb",
        ORIG_MT_fa + ".ann",
        ORIG_MT_fa + ".bwt",
        ORIG_MT_fa + ".pac",
        ORIG_MT_fa + ".sa",
        SHIFTED_MT_fa + ".amb",
        SHIFTED_MT_fa + ".ann",
        SHIFTED_MT_fa + ".bwt",
        SHIFTED_MT_fa + ".pac",
        SHIFTED_MT_fa + ".sa",
    default_target: True


rule create_fai:
    input: fasta=ancient(REF_MALE),
    output: fai = REF_MALE_FAI,
    conda: "envs/preprocess.yaml"
    shell: "samtools faidx {input}"

rule create_fai_1:
    input: fasta=ancient(REF_FEMALE),
    output: fai = REF_FEMALE_FAI
    conda: "envs/preprocess.yaml"
    shell: "samtools faidx {input}"

rule create_fai_chrM_shifted:
    input: mt_ref_shift = ancient(SHIFTED_MT_fa)
    output: fai = SHIFTED_MT_fai
    conda: "envs/preprocess.yaml"
    shell: "samtools faidx {input}"

rule create_dict_for_chrM_orig_reference:
    input: mt_ref_shift = ancient(ORIG_MT_fa)
    output: fai = ORIG_MT_fai
    conda: "envs/preprocess.yaml"
    shell: "samtools faidx {input}"

rule create_dict:
    input: fasta=ancient(REF_MALE),
    output: dict = REF_MALE_DICT
    conda: CONDA_VCF
    shell: "gatk CreateSequenceDictionary -R {input}"

rule create_dict_1:
    input: fasta = ancient(REF_FEMALE),
    output: dict = REF_FEMALE_DICT
    conda: CONDA_VCF
    shell: "gatk CreateSequenceDictionary -R {input}"

rule create_dict_for_chrM_shifted_reference:
    input: mt_ref_shift = ancient(SHIFTED_MT_fa)
    output: dict = SHIFTED_MT_dict
    conda: CONDA_VCF
    shell:  """
            gatk CreateSequenceDictionary -R {input}
            """

rule create_dict_for_chrM_reference:
    input: mt_ref_shift = ancient(ORIG_MT_fa)
    output: dict = ORIG_MT_dict
    conda: CONDA_VCF
    shell:  """
            gatk CreateSequenceDictionary -R {input}
            """

rule bwa_index_chrM_orig:
    input:
        fa=ancient(ORIG_MT_fa)
    output:
        amb=ORIG_MT_fa + ".amb",
        ann=ORIG_MT_fa + ".ann",
        bwt=ORIG_MT_fa + ".bwt",
        pac=ORIG_MT_fa + ".pac",
        sa=ORIG_MT_fa + ".sa"
    conda: CONDA_MAIN
    shell:
        """
        bwa index {input.fa}
        """

rule bwa_index_chrM_shifted:
    input:
        fa=ancient(SHIFTED_MT_fa)
    output:
        amb=SHIFTED_MT_fa + ".amb",
        ann=SHIFTED_MT_fa + ".ann",
        bwt=SHIFTED_MT_fa + ".bwt",
        pac=SHIFTED_MT_fa + ".pac",
        sa=SHIFTED_MT_fa + ".sa"
    conda: CONDA_MAIN
    shell:
        """
        bwa index {input.fa}
        """

rule create_hash:
    input:
        fasta=ancient(REF_MALE),
        dir=ancient(REF_MALE_DIR),
        bed=ancient(REF_MALE_BED),
        dict = ancient(REF_MALE_DICT)
    conda: "envs/preprocess.yaml"
    output: hash = REF_MALE_HASH
    params:
        dragmap = dragmap
    resources:
        mem_mb=50000,
        n=64
    shell:  """
            {params.dragmap} --build-hash-table true --ht-reference {input.fasta} --output-directory {input.dir} --ht-mask-bed {input.bed} 
            """

rule ComposeSTRTableFile:
    input:
        fasta=ancient(REF_MALE),
        dict= ancient(REF_MALE_DICT)
    output:
        str_file=REF_MALE_STR
    conda: "envs/vcf_handling.yaml"
    resources:
        mem_mb=12000,
        n=2
    shell:"""
        gatk ComposeSTRTableFile --java-options "-Xmx{resources.mem_mb}M" -R {input.fasta} -O {output.str_file}
        """

rule create_hash_1:
    input:
        fasta=ancient(REF_FEMALE),
        dir=ancient(REF_FEMALE_DIR),
        bed=ancient(REF_FEMALE_BED),
        dict= ancient(REF_FEMALE_DICT)
    conda: "envs/preprocess.yaml"
    output: hash = REF_FEMALE_HASH
    params:
        dragmap = dragmap
    resources:
        mem_mb=50000,
        n=64
    shell:  """
            {params.dragmap} --build-hash-table true --ht-reference {input.fasta} --output-directory {input.dir} --ht-mask-bed {input.bed} 
            """

rule ComposeSTRTableFile_1:
    input:
        fasta=ancient(REF_FEMALE),
        dict= ancient(REF_FEMALE_DICT)
    output:
        str_file = REF_FEMALE_STR
    conda: "envs/vcf_handling.yaml"
    resources:
        mem_mb=12000,
        n=2
    shell:"""
        gatk ComposeSTRTableFile --java-options "-Xmx{resources.mem_mb}M" -R {input.fasta} -O {output.str_file}
        """


