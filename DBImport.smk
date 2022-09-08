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

bins = config['RES'] + config['bin_file_ref']

rule DBImport_all:
    input:
        expand(["done_p{chr_p}.{chr}.txt"], zip, chr = main_chrs_db, chr_p = chr_p),
        rules.gVCF_all.input,
        # expand("{chr}_gvcfs.list", chr = main_chrs)
    default_target: True

rule make_glist:
    input: g= expand("{gvcfs}/{chr}/{sample}.{chr}.g.vcf.gz",gvcfs=config['gVCF'],sample=sample_names,chr=main_chrs),
    output: glist = temp("{chr}_gvcfs.list"),
    shell: "ls gvcf/{wildcards.chr}/*.g.vcf.gz > {output.glist}"


#Genomics DBImport instead CombineGVCFs
rule GenomicDBImport:
    input:
        g = expand("{gvcfs}/{chr}/{sample}.{chr}.g.vcf.gz", gvcfs = config['gVCF'], sample = sample_names, chr = main_chrs),
        glist = rules.make_glist.output.glist
    log: config['LOG'] + '/' + "GenomicDBImport.{chr_p}.{chr}.log"
    benchmark: config['BENCH'] + "/GenomicDBImport.{chr_p}.{chr}.txt"
    output:
        dbi=directory("genomicsdb_{chr}.p{chr_p}"),
        # gvcf_list = temp("{chr}_gvcfs.list"),
        ready = touch(temp('done_p{chr_p}.{chr}.txt'))
    threads: config['GenomicDBImport']['n']
    params:
        batches = '75',
        #gvcfs = lambda wildcards expand(" -V {gvcfs}/{chr}/{sample}.{chr}.g.vcf.gz", gvcfs = config['gVCF'], sample = sample_names, chr = wildcards.chr),
        intervals = config['RES'] + config['bin_file_ref'] + '/{chr}/hg38_mainchr_bins{chr_p}.bed.interval_list'
    priority: 30
    # conda: "preprocess"
    shell:
            # "ls gvcf/{wildcards.chr}/*.g.vcf.gz > {output.gvcf_list} && "
            # "{gatk} GenomicsDBImport --reader-threads {threads} -V {input.glist} \
            "{gatk} GenomicsDBImport --reader-threads {threads} -V {input.glist} \
                --intervals {params.intervals}  -R {ref} --genomicsdb-workspace-path {output.dbi} --batch-size {params.batches} \
             --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader 2> {log}"


# rule update_DBImport:
#     input:
#         dir = '/path/to/existing/DBI_dir'
#         sample_file = 'path/to/samplefiles/or/glists' (? produce glist inside rule: "ls gvcf/{wildcards.chr}/*.g.vcf.gz > {output.gvcf_list} &&" where gvcf path to folder with gvcfs and gvcf_list = "{chr}_gvcfs.list")
#     output: ready = touch(temp('done_p{chr_p}.{chr}.txt'))
#     log: config['LOG'] + '/' + "GenomicDBImport_update.{chr_p}.{chr}.log"
#     benchmark: config['BENCH'] + "/GenomicDBImport_update.{chr_p}.{chr}.txt"
#     params:
#         batches='75',
#         intervals=config['RES'] + config['bin_file_ref'] + '/{chr}/hg38_mainchr_bins{chr_p}.bed.interval_list'
#     priority: 30
#     conda: "preprocess"
#     shell: "{gatk} GenomicsDBImport --genomicsdb-update-workspace-path {input.dir} -V {?} --intervals {params.intervals} -R {ref} \ "
#             " --batch-size {params.batches} --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader 2> {log}"

