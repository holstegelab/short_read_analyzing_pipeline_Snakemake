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

module gVCF:
    snakefile: 'gVCF.smk'
    config: config
use rule * from gVCF

rule DBImport_all:
    input:
        expand("done_{chr}", chr = main_chrs),
    default_target: True

#Genomics DBImport instead CombineGVCFs
rule GenomicDBImport:
    input:
        expand("{gvcfs}/{chr}/{sample}.{chr}.g.vcf.gz", gvcfs = config['gVCF'], sample = sample_names, chr = main_chrs)
    log: config['LOG'] + '/' + "GenomicDBImport.{chr}.log"
    benchmark: config['BENCH'] + "/GenomicDBImport.{chr}.txt"
    output:
        dbi=directory("genomicsdb_{chr}"),
        gvcf_list = temp("{chr}_gvcfs.list"),
        ready = touch(temp('done_{chr}'))
    threads: config['GenomicDBImport']['n']
    params: batches = '75'
    priority: 30
    conda: "preprocess"
    shell:
        "ls gvcf/{wildcards.chr}/*.g.vcf.gz > {output.gvcf_list} && {gatk} GenomicsDBImport --reader-threads {threads} \
        -V {wildcards.chr}_gvcfs.list --intervals {wildcards.chr}  -R {ref} --genomicsdb-workspace-path {output.dbi} --batch-size {params.batches} \
         --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader 2> {log}"

# rule update_DBImport:
#     input:
#         dir = rules.GenomicDBImport.output.dbi,
#         sample_file = '?'
#     output: touch(temp('done_update.{chr}'))
#     log: config['LOG'] + '/' + "GenomicDBImport.{chr}.log"
#     benchmark: config['BENCH'] + "/GenomicDBImport_update.{chr}.txt"
#     params: batches = '75'
#     priority: 30
#     conda: "preprocess"
#     shell: "{gatk} GenomicsDBImport --genomicsdb-update-workspace-path {input.dir} -V {input.sample_file} --intervals {chr} -R {ref} \ "
#             " --batch-size {params.batches} --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader 2> {log}"