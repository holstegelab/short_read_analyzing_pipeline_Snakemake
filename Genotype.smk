import os
from itertools import repeat
wildcard_constraints:
    sample="[\w\d_\-@]+",

from common import *

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module DBImport:
    snakefile: 'DBImport.smk'
    config: config
module Combine_gVCF:
    snakefile: 'Combine_gVCF.smk'
    config: config


gVCF_combine_method = config.get("Combine_gVCF_method", "DBIMPORT")
mode = config.get("computing_mode", "WES")


#parts = level2_regions
#work around problem in genomicdbimport with many contigs
level2_range_special = [('A', 2,x,2) for x in range(0,99)] + \
               [('X', 1,x,2) for x in range(0,5)] + \
               [('X', 1, x, 1) for x in range(0, 5)] + \
               [('Y', 1, x,2) for x in range(0,2)] + \
               [('Y', 1, x, 1) for x in range(0, 2)] 
parts = get_regions(level2_range_special)


regions = []
for part in parts:
    regions.append(convert_to_level0(part))

print(regions, parts)

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
        [f"{VCF}/merged_{region}.p{part}.{mode}.vcf.gz" for region, part in zip(regions,parts)]
        #expand("{stat}/BASIC.{chr}.{mode}.variant_calling_detail_metrics", stat = config['STAT'], mode = mode, chr = main_chr),
        #expand("{vcf}/PER_chr/{chr}_{mode}_merged.vcf.gz",  vcf=config['VCF'], mode = mode, chr = main_chr),
        #expand("{vcf}/PER_chr/{chr}_{mode}_merged.vcf.gz.tbi", vcf=config['VCF'], mode = mode, chr = main_chr),
        #expand("{vcf}/{chr}/Merged_norm.{chr}.{mode}.vcf.gz", vcf = config['VCF_Final'], mode = mode, chr = main_chr),
        #expand("{vcf}/{chr}/Merged_norm.{chr}.{mode}.vcf.gz.tbi", vcf = config['VCF_Final'], mode = mode, chr = main_chr)
    default_target: True


def get_mem_mb_genotype(wildcrads, attempt):
    return attempt*int(20000)

def region_to_IL_file(wildcards):#{{{
    """Converts a region to a interval_list file location (see common.py and Tools.smk)"""
    region = wildcards['part']
    # WGS files have fewer regions so DBI works faster and could use multiple cores
    return region_to_file(region, wgs=mode=='WGS', extension='interval_list')#}}}


#if gVCF_combine_method == "DBIMPORT":
rule GenotypeDBI:
    input:
        ready=rules.GenomicDBImport.output.ready
    output:
        raw_vcfDBI=pj(VCF, "merged_{region}.p{part}.{mode}.vcf.gz")
    params:
        dir = rules.GenomicDBImport.params.dbi,
        intervals = region_to_IL_file,
        ploidy = lambda wildcards: 1 if 'H' in wildcards['part'] else 2
    conda: CONDA_VCF
    resources: 
        mem_mb_java = get_mem_mb_genotype
        mem_mb=40000
    priority: 40
    shell:"""
        {gatk} GenotypeGVCFs --java-options "-Xmx{resources.mem_mb}M" -R {REF} -V gendb://{params.dir} -O {output} -D {DBSNP} --intervals {params.intervals} -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation -A StrandBiasBySample -A AssemblyComplexity --annotate-with-num-discovered-alleles --genomicsdb-shared-posixfs-optimizations  --ploidy {params.ploidy}
        """


#elif gVCF_combine_method == "COMBINE_GVCF":
#    rule GenotypeDBI:
#        input: 
#            gvcf = rules.combinegvcfs.output.chr_gvcfs
#        output: 
#            raw_vcf = pj(VCF, "/merged_{region}.p{part}.{mode}.vcf.gz"
#        conda: CONDA_VCF
#        resources: mem_mb = get_mem_mb_genotype
#        priority: 40
#        shell:"""
#            {gatk} GenotypeGVCFs -R {ref} -V {input.gvcf} -O {output} -D {DBSNP} --intervals {wildcards.chr}
#            """




rule merge_per_chr:
    input:
        vcfs = lambda wildcards: expand("{dir}/Merged_raw_DBI_{chr}.p{{chr_p}}.{mode}.vcf.gz".format(dir=VCF, chr=wildcards.chr, mode=mode), chr_p=valid_chr_p[wildcards.chr])
    output: 
        per_chr_vcfs = os.path.join(VCF, "PER_chr", "{chr}_{mode}_merged.vcf.gz")
    params: 
        inputs = lambda wildcards: expand("-I {dir}/Merged_raw_DBI_{chr}.p{{chr_p}}.{mode}.vcf.gz".format(dir=VCF, chr=wildcards.chr, mode=mode), chr_p=valid_chr_p[wildcards.chr])
    priority: 45
    conda: "envs/vcf_handling.yaml"
    resources: 
        mem_mb = 400
    shell: """
            {gatk} GatherVcfs {params.inputs} -O {output} -R {REF}
            """


rule index_per_chr:
    input: rules.merge_per_chr.output.per_chr_vcfs
    output:
        vcfidx = os.path.join(VCF, "PER_chr", "{chr}_{mode}_merged.vcf.gz.tbi")
    priority: 46
    conda: "envs/preprocess.yaml"
    resources: 
        mem_mb=450
    shell:
        "{gatk} IndexFeatureFile -I {input} -O {output}"


rule norm_per_chr:
    input:
        vcf = rules.merge_per_chr.output.per_chr_vcfs,
        idx = rules.index_per_chr.output.vcfidx
    output:
        normVCF=pj(VCF, "{chr}/Merged_norm.{chr}.{mode}.vcf.gz"),
    priority: 50
    conda: "envs/preprocess.yaml"
    threads: 150
    shell:
        "({bcftools} norm --threads {threads} -f {REF} {input} -m -both -O v | {bcftools} norm --threads {threads} --check-ref ws -d exact -f {REF} -O z > {output.normVCF})"

rule norm_idx_per_chr:
    input:
        normVCF = rules.norm_per_chr.output.normVCF
    output:
        idx=pj(VCF,"{chr}/Merged_norm.{chr}.{mode}.vcf.gz.tbi")
    priority: 55
    conda: CONDA_VCF
    shell:
        "{gatk} IndexFeatureFile -I {input.normVCF} -O {output.idx}"


# basic stats
# include hom-het ratio, titv ratio, etc.
rule basic_stats_per_chr:
    input:
        vcf = rules.norm_per_chr.output.normVCF,
        tbi = rules.norm_idx_per_chr.output.idx
    output:
        os.path.join(STAT, "BASIC.{chr}.{mode}.variant_calling_detail_metrics"),
        os.path.join(STAT, "BASIC.{chr}.{mode}.variant_calling_summary_metrics")
    priority: 90
    conda: CONDA_VCF
    resources:
        n=2,
        mem_mb=4500
    shell:"""{gatk} CollectVariantCallingMetrics \
        -R {REF} -I {input.vcf} -O stats/BASIC.{wildcards.chr}.{wildcards.mode} \
        --DBSNP {DBSNP} --THREAD_COUNT {resources.n}"""



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
