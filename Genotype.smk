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
module DBImport:
    snakefile: 'DBImport.smk'
    config: config
module Combine_gVCF:
    snakefile: 'Combine_gVCF.smk'
    config: config


bins = config['RES'] + config['bin_file_ref']
gVCF_combine_method = config.get("Combine_gVCF_method", "COMBINE_GVCF")
if gVCF_combine_method == "DBIMPORT":
    use rule * from DBImport
    rule_all_combine = rules.DBImport_all.input
elif gVCF_combine_method == "COMBINE_GVCF":
    use rule * from Combine_gVCF
    rule_all_combine = rules.Combine_gVCF_all.input
else:
    raise ValueError(
        "invalid option provided to 'Combine_gVCF_method'; please choose either 'COMBINE_GVCF' or 'DBIMPORT'."
    )

rule Genotype_all:
    input:
        rule_all_combine,
        # expand(config['VCF'] + "/Merged_raw_DBI_{chr_p}.vcf.gz", chr_p = chr_p),
        expand("{vcf}/ALL_chrs.vcf.gz", vcf=config['VCF']),
        expand("{vcf}/Merged_norm.vcf", vcf=config['VCF_Final']),
    default_target: True

if gVCF_combine_method == "DBIMPORT":
    #use rule * from DBImport
    DBImethod = config.get("DBI_method", "new")
    merged_input = expand(["{dir}/Merged_raw_DBI_{chr}.p{chr_p}.vcf.gz"],zip,chr=main_chrs_db,chr_p=chr_p, dir = [config['VCF']]*99)
    vcfs_to_merge = expand(["-I {dir}/Merged_raw_DBI_{chr}.p{chr_p}.vcf.gz"],zip,chr=main_chrs_db,chr_p=chr_p,dir=[config['VCF']] * 99)
    if DBImethod == "new":
        rule GenotypeDBI:
            input:
                test = rules.GenomicDBImport.output.ready
            output:
                raw_vcfDBI=config['VCF'] + "/Merged_raw_DBI_{chr}.p{chr_p}.vcf.gz"
            log: config['LOG'] + '/' + "GenotypeDBI_{chr}.p{chr_p}.log"
            benchmark: config['BENCH'] + "/GenotypeDBI_{chr}.p{chr_p}.txt"
            params:
                dir = rules.GenomicDBImport.params.dbi,
                dbsnp = config['RES'] + config['dbsnp'],
                intervals= config['RES'] + config['bin_file_ref'] + '/{chr}/hg38_mainchr_bins{chr_p}.bed.interval_list'
            conda: "envs/preprocess.yaml"
            priority: 40
            shell:
                    "{gatk} GenotypeGVCFs -R {ref} -V gendb://{params.dir} -O {output} -D {params.dbsnp} --intervals {params.intervals} 2> {log}"
    elif DBImethod == "update":
        rule GenotypeDBI:
            input:
                gvcf = rules.SelectVariants_For_Genotype.output.gvcf,
                test = rules.GenomicDBImport.output.ready
            output:
                raw_vcfDBI=config['VCF'] + "/Merged_raw_DBI_{chr}.p{chr_p}.vcf.gz"
            log: config['LOG'] + '/' + "GenotypeDBI_{chr}.p{chr_p}.log"
            benchmark: config['BENCH'] + "/GenotypeDBI_{chr}.p{chr_p}.txt"
            params:
                dbsnp = config['RES'] + config['dbsnp'],
                intervals= config['RES'] + config['bin_file_ref'] + '/{chr}/hg38_mainchr_bins{chr_p}.bed.interval_list'
            conda: "envs/preprocess.yaml"
            priority: 40
            shell:
                    "{gatk} GenotypeGVCFs -R {ref} -V {input.gvcf} -O {output} -D {params.dbsnp} --intervals {params.intervals} 2> {log}"
    else:
        raise ValueError(
            "invalid option provided to 'DBImethod'; please choose either 'new' or 'update'."
        )
elif gVCF_combine_method == "COMBINE_GVCF":
    #use rule * from Combine_gVCF
    merged_input = expand("{dir}/MERGED/cohort_{chr}.g.vcf.gz", dir = config['VCF'], chr = main_chrs)
    vcfs_to_merge = expand("-I {dir}/MERGED/cohort_{chr}.g.vcf.gz", dir = config['VCF'], chr = main_chrs)
    rule GenotypeDBI:
        input: gvcf = rules.combinegvcfs.output.chr_gvcfs
        output: raw_vcf = config['VCF'] + "/MERGED/cohort_{chr}.g.vcf.gz"
        log: config['LOG'] + '/' + "GenotypeDBI_{chr}.log"
        benchmark: config['BENCH'] + "/GenotypeDBI_{chr}.txt"
        params:
            dbsnp=config['RES'] + config['dbsnp'],
        conda: "envs/preprocess.yaml"
        priority: 40
        shell:
            "{gatk} GenotypeGVCFs -R {ref} -V {input.gvcf} -O {output} -D {params.dbsnp} --intervals {wildcards.chr} 2> {log}"


rule Mergechrs:
    input:
        merged_input
        # expand(["{dir}/Merged_raw_DBI_{chr}.p{chr_p}.vcf.gz"],zip,chr=main_chrs_db,chr_p=chr_p, dir = [config['VCF']]*99)
    params: vcfs = vcfs_to_merge
    conda: "envs/preprocess.yaml"
    log: config['LOG'] + '/' + "Mergechrs.log"
    benchmark: config['BENCH'] + "/Mergechrs.txt"
    output:
        vcf = config['VCF'] + "/ALL_chrs.vcf.gz"
    priority: 45
    shell:
        "{gatk} GatherVcfs {params.vcfs} -O {output} -R {ref} 2> {log}"

rule index_mergechrs:
    input: rules.Mergechrs.output.vcf
    conda: "envs/preprocess.yaml"
    log: config['LOG'] + '/' + "Mergechrsed_index.log"
    benchmark: config['BENCH'] + "/Mergechrsed_index.txt"
    output:
        vcfidx = config['VCF'] + "/ALL_chrs.vcf.gz.idx"
    priority: 45
    shell:
        "{gatk} IndexFeatureFile -I {input} -O {output} 2> {log}"


rule norm:
    input:
        vcf = rules.Mergechrs.output.vcf,
        idx = rules.index_mergechrs.output.vcfidx
    output:
        normVCF=config['VCF_Final'] + "/Merged_norm.vcf",
    log: config['LOG'] + '/' + "normalization.log"
    benchmark: config['BENCH'] + "/normalization.txt"
    priority: 80
    conda: "envs/preprocess.yaml"
    shell:
        "({bcftools} norm -f {ref} {input} -m -both -O v | {bcftools} norm -d exact -f {ref} > {output.normVCF}) 2> {log}"

rule norm_idx:
    input:
        normVCF = rules.norm.output.normVCF
    output:
        idx=config['VCF_Final'] + "/Merged_norm.vcf.idx"
    log: config['LOG'] + '/' + "normalization.log"
    benchmark: config['BENCH'] + "/idx_normalizated.txt"
    priority: 80
    conda: "envs/preprocess.yaml"
    shell:
        "{gatk} IndexFeatureFile -I {input.normVCF} -O {output.idx} 2> {log}"


