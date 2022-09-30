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
    readgroup="[\w\d_\-@]+",

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
use rule * from gVCF

rule Combine_gVCF_all:
    input: expand("{gvcf}/MERGED/cohort_{chr}.g.vcf.gz", gvcf = config['gVCF'], chr = main_chrs)
    default_target: True

rule reblock_gvcf:
    input: gvcf = rules.HaplotypeCaller.output.gvcf
    output: gvcf_reblock = config['gVCF'] + "/reblock/{chr}/{sample}.{chr}.g.vcf.gz"
    log: Reblock=config['LOG'] + "/{sample}_{chr}_reblock.log"
    benchmark:
        config['BENCH'] + "/{sample}_{chr}_reblock.txt"
    conda: "envs/preprocess.yaml"
    params:
        dbsnp=config['RES'] + config['dbsnp'],
    shell:
        "{gatk} ReblockGVCF --keep-all-alts -D {params.dbsnp} -R {ref} -V {input.gvcf} -O {output.gvcf_reblock} -G StandardAnnotation -G AS_StandardAnnotation 2> {log}"

rule combinegvcfs:
    input: expand("{gvcf}/reblock/{chr}/{sample}.{chr}.g.vcf.gz", gvcf = config['gVCF'], sample = sample_names, allow_missing=True)
    output: chr_gvcfs = config['gVCF'] + "/MERGED/cohort_{chr}.g.vcf.gz"
    log: combine =config['LOG'] + "/{chr}_combine.log"
    benchmark:
        config['BENCH'] + "/{chr}_reblock.txt"
    conda: "envs/preprocess.yaml"
    params: inputs = expand("--variant {gvcf}/reblock/{chr}/{sample}.{chr}.g.vcf.gz", gvcf = config['gVCF'], sample = sample_names, allow_missing=True)
    shell:
        "{gatk} CombineGVCFs -G StandardAnnotation -G AS_StandardAnnotation {params.inputs} -O {output} -R {ref} 2> {log}"