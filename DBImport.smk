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
verifybamid2 = config['verifybamid2']
cur_dir = os.getcwd()

ref = config['RES'] + config['ref']

wildcard_constraints:
    sample="[\w\d_\-@]+",
    mode = "WES|WGS"
    # readgroup="[\w\d_\-@]+",

from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO = load_samplefiles('.', config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()
tmpdir_alternative = os.path.join(config['tmpdir'],getpass.getuser())
tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())

module Aligner:
    snakefile: 'Aligner.smk'
    config: config

module gVCF:
    snakefile: 'gVCF.smk'
    config: config
use rule * from gVCF

bins = config['RES'] + config['bin_file_ref']
mode = config.get("computing_mode", "WES")
rule DBImport_all:
    input:
        expand(["labels/done_p{chr_p}.{chr}.{mode}.txt"], zip, chr = main_chrs_db, chr_p = chr_p, mode = [mode]*853),
        rules.gVCF_all.input,
        # expand("{chr}_gvcfs.list", chr = main_chrs)
    default_target: True

DBImethod = config.get("DBI_method", "new")
DBIpath = config.get("DBIpath", "genomicsdb_")
if DBImethod == "new":
    # if want to
    DBI_method_params = "--genomicsdb-workspace-path "
    path_to_dbi = "genomicsdb_"
    labels = []
elif DBImethod == "update" and len(DBIpath) != 1:
    DBI_method_params = "--genomicsdb-update-workspace-path "
    path_to_dbi = DBIpath
    labels = expand(["labels/done_backup_{samplefile}_{mode}_{chr}.p{chr_p}"],zip,chr=main_chrs_db,chr_p=chr_p,mode=[mode] * 853,samplefile=SAMPLE_FILES * 853)
elif DBImethod == "update" and len(DBIpath) == 1:
    raise ValueError(
        "If you want to update existing DB please provide path to this DB in format 'DBIpath=/path/to/directory_with_DB-s/genomicsdb_'"
        "Don't provide {chr}.p{chr_p} part of path!"
    )
else:
    raise ValueError(
        "invalid option provided to 'DBImethod'; please choose either 'new' or 'update'."
    )



rule backup_gdbi:
    input: gdbi = path_to_dbi + '{chr}.p{chr_p}_{mode}'
    output: label = touch(temp('labels/done_backup_{samplefile}_{mode}_{chr}.p{chr_p}'))
    params: tar = "{samplefile}_{mode}_gdbi_{chr}.p{chr_p}.tar.gz"
    shell: """
            mkdir -p BACKUPS/previous &&
            find . -maxdepth 2 -name '*_gdbi_{chr}.p{chr_p}.tar.gz' -type f -print0 | xargs -0r mv -t BACKUPS/previous/ && 
            tar -czv -f BACKUPS/{params.tar} {input}
            """


def get_mem_mb_GenomicDBI(wildcrads, attempt):
    return attempt*(config['GenomicDBImport']['mem'])


rule GenomicDBImport:
    input:
        g=expand("{cd}/{gvcfs}/reblock/{chr}/{sample}.{chr}.{mode}.g.vcf.gz",gvcfs=config['gVCF'],sample=sample_names, mode = [mode], cd = cur_dir, allow_missing=True),
        intervals=os.path.join(config['RES'],config['kit_folder'],'BINS','interval_list','{chr}_{chr_p}.interval_list'),
        labels = labels
    log: config['LOG'] + "/GenomicDBImport.{chr_p}.{chr}.{mode}.log"
    benchmark: config['BENCH'] + "/{chr}_{chr_p}.{mode}_GenomicDBImport.txt"
    conda: "envs/preprocess.yaml"
    output:
        ready=touch(temp('labels/done_p{chr_p}.{chr}.{mode}.txt'))
    threads: config['GenomicDBImport']['n']
    params:
        inputs=expand(" -V {gvcfs}/reblock/{chr}/{sample}.{chr}.{mode}.g.vcf.gz",gvcfs=config['gVCF'],sample=sample_names,mode = [mode],allow_missing=True),
        dbi=os.path.join(path_to_dbi + "{chr}.p{chr_p}_{mode}"),
        method=DBI_method_params,
        batches='75',
    priority: 30
    resources: mem_mb = get_mem_mb_GenomicDBI,
            tmpdir= config['tmpdir']
    shell:
        """{gatk} GenomicsDBImport --java-options "-Xmx{resources.mem_mb}M"  --reader-threads {threads} {params.inputs} \
            --intervals {input.intervals}  -R {ref} {params.method} {params.dbi}/ --batch-size {params.batches} --tmp-dir {resources.tmpdir} \
         --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader 2> {log}"""
