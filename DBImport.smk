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

bins = config['RES'] + config['bin_file_ref']

rule DBImport_all:
    input:
        expand(["done_p{chr_p}.{chr}.txt"], zip, chr = main_chrs_db, chr_p = chr_p),
        rules.gVCF_all.input,
        # expand("{chr}_gvcfs.list", chr = main_chrs)
    default_target: True

rule make_glist:
    input: g= expand("{gvcfs}/{chr}/{sample}.{chr}.g.vcf.gz",gvcfs=config['gVCF'],sample=sample_names,chr=main_chrs),
    output: glist = ("{chr}_gvcfs.list"),
    shell: "ls gvcf/{wildcards.chr}/*.g.vcf.gz > {output.glist}"


DBImethod = config.get("DBI_method", "new")
DBIpath = config.get("DBIpath", "")
if DBImethod == "new":
    DBI_method_params = "--genomicsdb-workspace-path "
    path_to_dbi = "genomicsdb_"
elif DBImethod == "update" and len(DBIpath) != 0:
    DBI_method_params = "--genomicsdb-update-workspace-path "
    path_to_dbi = DBIpath
elif DBImethod == "update" and len(DBIpath) == 0:
    raise ValueError(
        "If you want to update existing DB please provide path to this DB in format 'DBIpath=/path/to/directory_with_DB-s/genomicsdb_'"
        "Don't provide {chr}.p{chr_p} part of path!"
    )
else:
    raise ValueError(
        "invalid option provided to 'DBImethod'; please choose either 'new' or 'update'."
    )

#Genomics DBImport instead CombineGVCFs
rule GenomicDBImport:
    input:
        g = expand("{gvcfs}/{chr}/{sample}.{chr}.g.vcf.gz", gvcfs = config['gVCF'], sample = sample_names, chr = main_chrs),
        glist = rules.make_glist.output.glist
    log: config['LOG'] + '/' + "GenomicDBImport.{chr_p}.{chr}.log"
    benchmark: config['BENCH'] + "/GenomicDBImport.{chr_p}.{chr}.txt"
    conda: "envs/preprocess.yaml"
    output:
        #dbi = directory(os.path.join(path_to_dbi + "{chr}.p{chr_p}")),
        # dbi=directory("genomicsdb_{chr}.p{chr_p}"),
        # gvcf_list = temp("{chr}_gvcfs.list"),
        ready = touch(temp('done_p{chr_p}.{chr}.txt'))
    threads: config['GenomicDBImport']['n']
    params:
        dbi = os.path.join(path_to_dbi + "{chr}.p{chr_p}"),
        method = DBI_method_params,
        batches = '75',
        #gvcfs = lambda wildcards expand(" -V {gvcfs}/{chr}/{sample}.{chr}.g.vcf.gz", gvcfs = config['gVCF'], sample = sample_names, chr = wildcards.chr),
        intervals = config['RES'] + config['bin_file_ref'] + '/{chr}/hg38_mainchr_bins{chr_p}.bed.interval_list'
    priority: 30
    # conda: "envs/preprocess.yaml"
    shell:
            # "ls gvcf/{wildcards.chr}/*.g.vcf.gz > {output.gvcf_list} && "
            # "{gatk} GenomicsDBImport --reader-threads {threads} -V {input.glist} \
            "{gatk} GenomicsDBImport --reader-threads {threads} -V {input.glist} \
                --intervals {params.intervals}  -R {ref} {params.method} {params.dbi}/ --batch-size {params.batches} \
             --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader 2> {log}"

if DBImethod == "update":
    rule SelectVariants_For_Genotype:
        input:
            gdbi = rules.GenomicDBImport.params.dbi,
            ready = rules.GenomicDBImport.output.ready
        output: gvcf= config['gVCF'] + "/SELECTED/{chr}_p{chr_p}.g.vcf"
        log: config['LOG'] + '/' + "SelectVariants_For_Genotype.{chr_p}.{chr}.log"
        benchmark: config['BENCH'] + "/SelectVariants_For_Genotype.{chr_p}.{chr}.txt"
        params: SN = expand("-sn {samples}", samples = sample_names)
        conda: "envs/preprocess.yaml"
        shell:
            "{gatk} SelectVariants -R {ref} -V gendb:///{input.gdbi} {params.SN} -O {output}"

