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



from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

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
gVCF_combine_method = config.get("Combine_gVCF_method", "DBIMPORT")
mode = config.get("computing_mode", "WES")
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
        # expand(["{vcf}/Merged_raw_DBI_{chr}.p{chr_p}.{mode}.vcf.gz"],zip,chr=main_chrs_db,chr_p=chr_p, vcf = [config['VCF']]*853, mode = [mode]*853),
        # expand("{vcf}/ALL_chrs.{mode}.vcf.gz", vcf=config['VCF'], mode = mode),
        expand("{stat}/BASIC.{chr}.{mode}.variant_calling_detail_metrics", stat = config['STAT'], mode = mode, chr = chr),
        expand("{vcf}/PER_chr/{chr}_{mode}_merged.vcf.gz",  vcf=config['VCF'], mode = mode, chr = chr),
        expand("{vcf}/PER_chr/{chr}_{mode}_merged.vcf.gz.tbi", vcf=config['VCF'], mode = mode, chr = chr),
        expand("{vcf}/{chr}/Merged_norm.{chr}.{mode}.vcf.gz", vcf = config['VCF_Final'], mode = mode, chr = chr),
        expand("{vcf}/{chr}/Merged_norm.{chr}.{mode}.vcf.gz.tbi", vcf = config['VCF_Final'], mode = mode, chr = chr)
    default_target: True


def get_mem_mb_genotype(wildcrads, attempt):
    return attempt*int(config['GenotypeDBI']['mem'])

def get_parts_capture_kit(wildcards):
    chr = wildcards.chr
    parts = wildcards.chr_p
    if mode == 'WES':
        capture_kit_parts_path = os.path.join(config['RES'], config['kit_folder'], config['MERGED_CAPTURE_KIT'], 'interval_list', config['MERGED_CAPTURE_KIT'] + '_' + chr + '_' + parts + '.interval_list')
    else:
        capture_kit_parts_path = os.path.join(config['RES'], config['kit_folder'], 'BINS', 'interval_list', chr + '_' + parts + '.interval_list')
    return capture_kit_parts_path

if gVCF_combine_method == "DBIMPORT":
    #use rule * from DBImport
    DBImethod = config.get("DBI_method", "new")
    rule GenotypeDBI:
        input:
            test = rules.GenomicDBImport.output.ready
        output:
            raw_vcfDBI=config['VCF'] + "/Merged_raw_DBI_{chr}.p{chr_p}.{mode}.vcf.gz"
        log: config['LOG'] + '/' + "GenotypeDBI_{chr}.p{chr_p}.{mode}.log"
        benchmark: config['BENCH'] + "/GenotypeDBI_{chr}.p{chr_p}.{mode}.txt"
        params:
            dir = rules.GenomicDBImport.params.dbi,
            dbsnp = config['RES'] + config['dbsnp'],
            intervals = get_parts_capture_kit
        conda: "envs/preprocess.yaml"
        resources: mem_mb = get_mem_mb_genotype
        priority: 40
        shell:
                """{gatk} GenotypeGVCFs --java-options "-Xmx{resources.mem_mb}M" -R {ref} -V gendb://{params.dir} -O {output} -D {params.dbsnp} --intervals {params.intervals} 2> {log}"""


elif gVCF_combine_method == "COMBINE_GVCF":
    #use rule * from Combine_gVCF
    rule GenotypeDBI:
        input: gvcf = rules.combinegvcfs.output.chr_gvcfs
        output: raw_vcf = config['VCF'] + "/Merged_raw_DBI_{chr}.p{chr_p}.{mode}.vcf.gz"
        log: config['LOG'] + '/' + "GenotypeDBI_{chr}.{chr_p}.{mode}.log"
        benchmark: config['BENCH'] + "/GenotypeDBI_{chr}.{chr_p}.{mode}.txt"
        params:
            dbsnp=config['RES'] + config['dbsnp'],
        conda: "envs/preprocess.yaml"
        resources: mem_mb = get_mem_mb_genotype
        priority: 40
        shell:
            "{gatk} GenotypeGVCFs -R {ref} -V {input.gvcf} -O {output} -D {params.dbsnp} --intervals {wildcards.chr} 2> {log}"




rule merge_per_chr:
    input:
        vcfs = lambda wildcards: expand("{dir}/Merged_raw_DBI_{chr}.p{{chr_p}}.{mode}.vcf.gz".format(dir=config['VCF'], chr=wildcards.chr, mode=mode), chr_p=valid_chr_p[wildcards.chr])
    output: per_chr_vcfs = os.path.join(config['VCF'], "PER_chr", "{chr}_{mode}_merged.vcf.gz")
    params: inputs = lambda wildcards: expand("-I {dir}/Merged_raw_DBI_{chr}.p{{chr_p}}.{mode}.vcf.gz".format(dir=config['VCF'], chr=wildcards.chr, mode=mode), chr_p=valid_chr_p[wildcards.chr])
    log: config['LOG'] + "/merge_per_chr_{chr}.{mode}.log"
    benchmark: config['BENCH'] + "/merge_per_chr_{chr}.{mode}.txt"
    priority: 45
    conda: "envs/preprocess.yaml"
    resources: mem_mb = config['merge_per_chr']['mem']
    shell:
        """{gatk} GatherVcfs {params.inputs} -O {output} -R {ref} 2> {log}"""


rule index_per_chr:
    input: rules.merge_per_chr.output.per_chr_vcfs
    conda: "envs/preprocess.yaml"
    log: config['LOG'] + "/Mergechrsed_index.{chr}.{mode}.log"
    benchmark: config['BENCH'] + "/Mergechrsed_index.{chr}.{mode}.txt"
    output:
        vcfidx = os.path.join(config['VCF'], "PER_chr", "{chr}_{mode}_merged.vcf.gz.tbi")
    priority: 46
    resources: mem_mb=config['index_per_chr']['mem']
    shell:
        "{gatk} IndexFeatureFile -I {input} -O {output} 2> {log}"


rule norm_per_chr:
    input:
        vcf = rules.merge_per_chr.output.per_chr_vcfs,
        idx = rules.index_per_chr.output.vcfidx
    output:
        normVCF=config['VCF_Final'] + "/{chr}/Merged_norm.{chr}.{mode}.vcf.gz",
    log: config['LOG'] + "/normalization_per_chr.{chr}.{mode}.log"
    benchmark: config['BENCH'] + "/normalization_per_chr.{chr}.{mode}.txt"
    priority: 50
    conda: "envs/preprocess.yaml"
    threads: config['norm_per_chr']['n']
    shell:
        "({bcftools} norm --threads {threads} -f {ref} {input} -m -both -O v | {bcftools} norm --threads {threads} --check-ref ws -d exact -f {ref} -O z > {output.normVCF}) 2> {log}"

rule norm_idx_per_chr:
    input:
        normVCF = rules.norm_per_chr.output.normVCF
    output:
        idx=config['VCF_Final'] + "/{chr}/Merged_norm.{chr}.{mode}.vcf.gz.tbi"
    log: config['LOG'] + "/idx_normalizated.{chr}.{mode}.log"
    benchmark: config['BENCH'] + "/idx_normalizated.{chr}.{mode}.txt"
    priority: 55
    conda: "envs/preprocess.yaml"
    shell:
        "{gatk} IndexFeatureFile -I {input.normVCF} -O {output.idx} 2> {log}"


# basic stats
# include hom-het ratio, titv ratio, etc.
rule basic_stats_per_chr:
    input:
        vcf = rules.norm_per_chr.output.normVCF,
        tbi = rules.norm_idx_per_chr.output.idx
    output:
        os.path.join(config['STAT'], "BASIC.{chr}.{mode}.variant_calling_detail_metrics"),
        os.path.join(config['STAT'], "BASIC.{chr}.{mode}.variant_calling_summary_metrics")
    priority: 90
    log: os.path.join(config['LOG'], "VCF_stats.{chr}.{mode}.log")
    benchmark: os.path.join(config['BENCH'], "VCF_stats.{chr}.{mode}.txt")
    params: dbsnp = os.path.join(config['RES'], config['dbsnp'])
    conda: "envs/preprocess.yaml"
    threads: config['basic_stats_per_chr']['n']
    shell:
        "{gatk} CollectVariantCallingMetrics \
        -R {ref} -I {input.vcf} -O stats/BASIC.{wildcards.chr}.{wildcards.mode} \
        --DBSNP {params.dbsnp} --THREAD_COUNT {threads} 2> {log}"



# rule Mergechrs:
#     input:
        # merged_input = expand(["{dir}/Merged_raw_DBI_{chr}.p{chr_p}.{mode}.vcf.gz"],zip,chr=main_chrs_db,chr_p=chr_p, dir = [config['VCF']]*853, mode = [mode]*853)
#     params: vcfs_to_merge = expand(["-I {dir}/Merged_raw_DBI_{chr}.p{chr_p}.{mode}.vcf.gz"],zip,chr=main_chrs_db,chr_p=chr_p,dir=[config['VCF']] * 853, mode = [mode]*853)
#     conda: "envs/preprocess.yaml"
#     log: config['LOG'] + "/Mergechrs.{mode}.log"
#     benchmark: config['BENCH'] + "/Mergechrs.{mode}.txt"
#     output:
#         vcf = config['VCF'] + "/ALL_chrs.{mode}.vcf.gz"
#     priority: 45
#     shell:
#         "{gatk} GatherVcfs {params.vcfs} -O {output} -R {ref} 2> {log}"
#
# rule index_mergechrs:
#     input: rules.Mergechrs.output.vcf
#     conda: "envs/preprocess.yaml"
#     log: config['LOG'] + '/' + "Mergechrsed_index.{mode}.log"
#     benchmark: config['BENCH'] + "/Mergechrsed_index.{mode}.txt"
#     output:
#         vcfidx = config['VCF'] + "/ALL_chrs.{mode}.vcf.gz.tbi"
#     priority: 46
#     shell:
#         "{gatk} IndexFeatureFile -I {input} -O {output} 2> {log}"
#
# rule norm:
#     input:
#         vcf = rules.Mergechrs.output.vcf,
#         idx = rules.index_mergechrs.output.vcfidx
#     output:
#         normVCF=config['VCF_Final'] + "/Merged_norm.{mode}.vcf",
#     log: config['LOG'] + '/' + "normalization.{mode}.log"
#     benchmark: config['BENCH'] + "/normalization.{mode}.txt"
#     priority: 50
#     conda: "envs/preprocess.yaml"
#     shell:
#         "({bcftools} norm -f {ref} {input} -m -both -O v | {bcftools} norm -d exact -f {ref} > {output.normVCF}) 2> {log}"
#
# rule norm_idx:
#     input:
#         normVCF = rules.norm.output.normVCF
#     output:
#         idx=config['VCF_Final'] + "/Merged_norm.{mode}.vcf.idx"
#     log: config['LOG'] + '/' + "normalization.{mode}.log"
#     benchmark: config['BENCH'] + "/idx_normalizated.{mode}.txt"
#     priority: 55
#     conda: "envs/preprocess.yaml"
#     shell:
#         "{gatk} IndexFeatureFile -I {input.normVCF} -O {output.idx} 2> {log}"
#
#
