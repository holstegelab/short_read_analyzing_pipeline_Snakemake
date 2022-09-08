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


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO = load_samplefiles('.', config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

rule SV_delly_all:
    input:
            config['DELLY'] + 'Filtred_SV.vcf',
            rules.Aligner_all.input,
    # input: expand("{delly}/Filtred_SV.vcf", delly = config['DELLY'])
    default_target: True

rule delly_call:
    input:
        bam = rules.markdup.output.mdbams
    output:
        call = config['DELLY'] + '/first_call/{sample}.bcf'
    conda: 'preprocess'
    shell: "delly call -g {ref} -o {output} {input}"

rule delly_merge:
    input: expand('{delly}/first_call/{sample}.bcf', delly = config['DELLY'], sample = sample_names)
    output: config['DELLY'] + 'Merged_sites.bcf'
    conda: 'preprocess'
    shell: "delly merge -o {output} {input}"

rule delly_genotype:
    input:
        bam = rules.markdup.output.mdbams,
        sites = rules.delly_merge.output
    output: calls = config['DELLY'] + '/geno_call/{sample}_geno.bcf'
    conda: 'preprocess'
    shell: "delly call -g {ref} -v {input.sites} -o {output.calls} {input.bam}"

rule bcf_merge:
    input: expand('{delly}/geno_call/{sample}_geno.bcf', delly = config['DELLY'], sample = sample_names)
    output: bcf_merge  = config['DELLY'] + '/Merged_Genotyped.bcf'
    conda: 'preprocess'
    shell: "bcftools merge -m id -O b -o {output.bcf_merge} {input}"

rule filter:
    input: rules.bcf_merge.output.bcf_merge
    output: bcf_filter = os.path.join(config['DELLY'], 'Filtred_SV.bcf')
    conda: 'preprocess'
    shell: "delly filter -f germline -o {output} {input}"

rule vcf_prod:
    input: rules.filter.output
    output: os.path.join(config['DELLY'], 'Filtred_SV.vcf')
    conda: 'preprocess'
    shell: "bcftools convert -O z -o {output} {input}"


