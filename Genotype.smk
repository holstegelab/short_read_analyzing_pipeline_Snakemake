configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")
gatk = config['miniconda'] + config['gatk']
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
module gVCF:
    snakefile: 'gVCF.smk'
    config: config
module DBImport:
    snakefile: 'DBImport.smk'
    config: config
use rule * from DBImport

rule Genotype_all:
    input:
        rules.gVCF_all.input,
        expand("{vcf}/ALL_chrs.vcf.gz", vcf=config['VCF']),
    default_target: True


# genotype
# multiple samplefiles
rule GenotypeDBI:
    input:
        rules.GenomicDBImport.output.dbi
    output:
        raw_vcfDBI=config['VCF'] + "/Merged_raw_DBI_{chr}.vcf.gz"
    log: config['LOG'] + '/' + "GenotypeDBI.{chr}.log"
    benchmark: config['BENCH'] + "/GenotypeDBI.{chr}.txt"
    params:
        dbsnp = config['RES'] + config['dbsnp']
    conda: "preprocess"
    priority: 40
    shell:
            "{gatk} GenotypeGVCFs -R {ref} -V gendb://{input} -O {output} -D {params.dbsnp} --intervals {wildcards.chr} 2> {log}"


rule Mergechrs:
    input:
        expand(config['VCF'] + "/Merged_raw_DBI_{chr}.vcf.gz", chr = chr)
    params:
        vcfs = expand("-I {dir}/Merged_raw_DBI_{chr}.vcf.gz", dir = config['VCF'], chr = chr)
    conda: "preprocess"
    log: config['LOG'] + '/' + "Mergechrs.log"
    benchmark: config['BENCH'] + "/Mergechrs.txt"
    output:
        vcf = config['VCF'] + "/ALL_chrs.vcf.gz"
    priority: 45
    shell:
        "{gatk} GatherVcfs {params.vcfs} -O {output} -R {ref} 2> {log} && {gatk} IndexFeatureFile -I {output} "