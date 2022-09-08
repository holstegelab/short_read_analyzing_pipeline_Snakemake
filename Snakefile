import pandas as pd
import read_stats
import os
configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")


gatk = config['miniconda'] + config['gatk']
samtools = config['miniconda'] + config['samtools']
bcftools = config['miniconda'] + config['bcftools']
dragmap = config['miniconda'] + config['dragmap']
cutadapt = config['miniconda'] + config['cutadapt']
verifybamid2 = config['miniconda'] + config['verifybamid2']
ref = config['RES'] + config['ref']

wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+"

# main chromosomes from GRCh38 splitted into 99 bins
bins = config['RES'] + config['bin_file_ref']
import csv
chr_p = []
with open(bins) as file:
    tsv_file = csv.reader(file, delimiter="\t")
    for line in tsv_file:
        if line[1] == str('0'):
            out = line[0] + ':' + str('1') + '-' + line[2]
        else:
            out = line[0] + ':' + line[1] + '-' + line[2]
        chr_p.append(out)
# print(chrs)

from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO = load_samplefiles('.', config)


# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
module gVCF:
    snakefile: 'gVCF.smk'
    config: config
module DBImport:
    snakefile: 'DBImport.smk'
    config: config
module Genotype:
    snakefile: 'Genotype.smk'
    config: config
module VQSR:
    snakefile: 'VQSR.smk'
    config: config
module Stat:
    snakefile: 'Stat.smk'
    config: config
module SV_delly:
    snakefile: 'SV_delly.smk'
    config: config
module CNV_with_cnvkit_Module:
    snakefile: 'CNV_with_cnvkit_Module.smk'
    config: config
use rule * from Aligner
use rule * from gVCF
use rule * from SV_delly
use rule * from CNV_with_cnvkit_Module
use rule * from DBImport
use rule * from Genotype
use rule * from VQSR
use rule * from Stat


rule all:
    input:
        rules.Aligner_all.input,
        rules.gVCF_all.input,
        rules.Genotype_all.input,
        rules.VQSR_all.input,
        rules.Stat_all.input,
        rules.SV_delly_all.input,
        rules.CNV_with_cnvkit_Module_all.input
    default_target: True

