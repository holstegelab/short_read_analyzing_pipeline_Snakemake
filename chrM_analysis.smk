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

ref = config['RES'] + config['ref']

tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())
tmpdir_alternative = os.path.join(config['tmpdir'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

wildcard_constraints:
    sample="[\w\d_\-@]+",


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

module Stat:
    snakefile: 'Stat.smk'
    config: config
use rule verifybamid from Stat

mode = config.get("computing_mode", "WES")
cur_dir = os.getcwd()

def generate_gvcf(wildcards):
    """Generate gvcf file name."""
    res = []
    for sample in sample_names:
        sinfo = SAMPLEINFO[sample]
        if sinfo['sex'] == 'F':
            for chrom in main_chrs_ploidy_female:
                res.append(os.path.join(cur_dir, config['gVCF'], chrom, sample + '.' + chrom + '.female.g.vcf.gz'))
        else:
            for chrom in main_chrs_ploidy_male:
                res.append(os.path.join(cur_dir, config['gVCF'], chrom, sample + '.' + chrom + '.male.g.vcf.gz'))
    return res

rule gVCF_all:
    input:
        generate_gvcf,
        rules.Aligner_all.input
    default_target: True