import pandas as pd
import read_stats
import os
configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")

gatk = config['gatk']
samtools = config['samtools']
bcftools = config['bcftools']
dragmap = config['dragmap']
verifybamid2 = config['verifybamid2']
ref = config['RES'] + config['ref']

wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+"

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
module Combine_gVCF:
    snakefile: 'Combine_gVCF.smk'
    config: config
use rule * from Aligner
use rule * from gVCF
use rule * from SV_delly
use rule * from CNV_with_cnvkit_Module
# use rule * from DBImport
use rule * from Genotype
use rule * from Stat
use rule * from VQSR

VQSR = config.get("VQSR", "NO")
if VQSR == "RUN_VQSR":
    VQSR_rule = rules.VQSR_all.input,
elif VQSR == "NO" or VQSR == "NO_VQSR" or VQSR == "NO_RUN":
    VQSR_rule = []
else:
    raise ValueError(
        "invalid option provided to 'VQSR'; please choose either 'RUN_VQSR' or 'NO_VQSR'."
    )

SV = config.get("SV", "RUN_SV")
if SV == "RUN_SV":
    SV_rule = rules.SV_delly_all.input
else:
    SV_rule = []

CNV = config.get("CNV", "RUN_CNV")
if CNV == "RUN_CNV":
    CNV_rule = rules.CNV_with_cnvkit_Module_all.input
else:
    CNV_rule = []

gVCF_combine_method = config.get("Combine_gVCF_method", "COMBINE_GVCF")
if gVCF_combine_method == "DBIMPORT":
    rule_all_combine = rules.DBImport_all.input
    use rule * from DBImport
elif gVCF_combine_method == "COMBINE_GVCF":
    rule_all_combine = rules.Combine_gVCF_all.input
    use rule * from Combine_gVCF
else:
    raise ValueError(
        "invalid option provided to 'Combine_gVCF_method'; please choose either 'COMBINE_GVCF' or 'DBIMPORT'."
    )

rule all:
    input:
        rules.Aligner_all.input,
        rules.gVCF_all.input,
        rules.Genotype_all.input,
        rule_all_combine,
        VQSR_rule,
        rules.Stat_all.input,
        SV_rule,
        CNV_rule
    default_target: True





