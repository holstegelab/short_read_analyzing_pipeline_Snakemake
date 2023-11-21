import pandas as pd
import read_stats
import os
import getpass
import read_samples
from common import *
import utils
current_dir = os.getcwd()

wildcard_constraints:
    sample="[\w\d_\-@]+",
    mode = "WES|WGS"
    # readgroup="[\w\d_\-@]+",

module Aligner:
    snakefile: 'Aligner.smk'
    config: config

module gVCF:
    snakefile: 'gVCF.smk'
    config: config
use rule * from gVCF

sample_types = config.get("sample_types","WES")
if sample_types == 'WES':
    parts = level2_regions
else:
    parts = level3_regions

regions = []
for part in parts:
    regions.append(convert_to_level0(part))
rule DBImport_all:
    input:
        expand(['labels/done_{region}.p{part}.txt'],zip,region = regions, part = parts),
        rules.gVCF_all.input,
        # expand("{chr}_gvcfs.list", chr = main_chrs)
    default_target: True

DBImethod = config.get("DBI_method", "new")
DBIpath = config.get("DBIpath", "genomicsdb_")

if DBImethod == "new":
    # if want to
    DBI_method_params = "--genomicsdb-workspace-path "
    path_to_dbi = DBIpath
    labels = []
elif DBImethod == "update" and len(DBIpath) != 1:
    DBI_method_params = "--genomicsdb-update-workspace-path "
    path_to_dbi = DBIpath
    number_of_splits = len(regions)
    labels = expand(["labels/done_backup_{samplefile}_{region}.p{part}"],zip, region = regions, part = parts,samplefile=SAMPLE_FILES * number_of_splits)
elif DBImethod == "update" and len(DBIpath) == 1:
    raise ValueError(
        "If you want to update existing DB please provide path to this DB in format 'DBIpath=/path/to/directory_with_DB-s/genomicsdb_'"
        "Don't provide {chr}.p{chr_p} part of path!"
    )
else:
    raise ValueError(
        "invalid option provided to 'DBImethod'; please choose either 'new' or 'update'."
    )
def region_to_IL_file(wildcards):#{{{
    """Converts a region to a bed file location (see common.py and Tools.smk)"""
    # sample = wildcards['sample']
    region = wildcards['parts']
    # return region_to_file(region, wgs='wgs' in SAMPLEINFO[sample]['sample_type'], extension='bed')#}}}
    if sample_types == 'WES':
        return region_to_file(region,extension='interval_list')
    elif sample_types == 'WGS':
        return region_to_file(region, wgs=True, extension='interval_list')#}}}




rule backup_gdbi:
    input: gdbi = path_to_dbi + '{region}.p{part}'
    output: label = touch(temp('labels/done_backup_{samplefile}_{region}.p{part}'))
    params: tar = "{samplefile}_gdbi_{region}.p{part}.tar.gz"
    shell: """
            mkdir -p BACKUPS/previous &&
            find . -maxdepth 2 -name '*_gdbi_{region}.p{part}.tar.gz' -type f -print0 | xargs -0r mv -t BACKUPS/previous/ && 
            tar -czv -f BACKUPS/{params.tar} {input}
            """


def get_mem_mb_GenomicDBI(wildcrads, attempt):
    return attempt*3500


rule GenomicDBImport:
    input:
        g=expand("{cd}/{GVCF}/{region}/{sample}.{region}.wg.vcf.gz",cd = current_dir, GVCF = GVCF, sample=sample_names,allow_missing=True),
        labels = labels
    log: pj(LOG, "GenomicDBImport.{region}.p{part}.log")
    benchmark: pj(BENCH, "{region}.p{part}_GenomicDBImport.txt")
    conda: CONDA_VCF
    output:
        ready=touch(temp('labels/done_{region}.p{part}.txt'))
    threads: 2
    params:
        inputs=expand(" -V {cd}/{GVCF}/{region}/{sample}.{region}.wg.vcf.gz",cd = current_dir, GVCF = GVCF, sample=sample_names,allow_missing=True),
        dbi=os.path.join(path_to_dbi + "{region}.p{part}"),
        method=DBI_method_params,
        batches='75',
        intervals = region_to_IL_file,
        ref = REF_MALE
    priority: 30
    resources: mem_mb = get_mem_mb_GenomicDBI,
            tmpdir= config['tmpdir']
    shell:
        """{gatk} GenomicsDBImport --java-options "-Xmx{resources.mem_mb}M"  --reader-threads {threads} {params.inputs} \
            --intervals {params.intervals}  -R {params.ref} {params.method} {params.dbi}/ --batch-size {params.batches} --tmp-dir {resources.tmpdir} \
         --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader 2> {log}"""
