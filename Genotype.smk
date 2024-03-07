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


gVCF_combine_method = config.get("gvcf_combine_method", "DBIMPORT")
genotype_mode = config.get("genotype_mode", "WES") #or WGS
genotype_alg = config.get("genotype_alg", "GenotypeGVCFs") #or GnarlyGenotyper
genotype_level = int(config.get("genotype_level", 3))
DBIpath = config.get("DBIpath", "genomicsdb_")


print(f"Genotype algorithm: {genotype_alg}")
print(f"Genotype level: {genotype_level}")
print(f"Genotype caller: {genotype_alg}")

if genotype_level == 2:
    parts = get_regions(level2_range_so)
elif genotype_level == 3:
    parts = get_regions(level3_range_so)
elif genotype_level == 4:
    parts = get_regions(level4_range_so)
else:
    raise RuntimeError(f'Unknown level {gneotype_level}')


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
        [f"{genotype_alg}/{VCF}/merged_{region}.{genotype_mode}.vcf.gz" for region in parts],
        [f"{genotype_alg}/{VCF}/merged_{region}.{genotype_mode}.vcf.gz.tbi" for region in parts],
        [f"{genotype_alg}/{VCF}/merged_{region}.{genotype_mode}.rescaled.vcf.gz" for region in parts],
        # [f"{genotype_alg}/{VCF}/ANNOTATED/merged_{region}.{genotype_mode}_annotated.vcf.gz" for region in parts]
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
    region = wildcards['region']
    # WGS files have fewer regions so DBI works faster and could use multiple cores
    return region_to_file(region, wgs=genotype_mode=='WGS', classic=True, padding=False, extension='interval_list')#}}}


#if gVCF_combine_method == "DBIMPORT":
rule GenotypeDBI:
    input:
        ready='labels/done_p{region}.txt',
        dir=DBIpath + "p{region}",
        intervals = region_to_IL_file
    output:
        raw_vcfDBI=ensure(expand(pj("{genotype_alg}",VCF, "merged_{region}.{genotype_mode}.vcf.gz"), genotype_alg = genotype_alg, allow_missing=True), non_empty=True),
        tbi = ensure(expand(pj("{genotype_alg}",VCF, "merged_{region}.{genotype_mode}.vcf.gz.tbi"), genotype_alg = genotype_alg, allow_missing=True), non_empty=True)
    params:
        ploidy = lambda wildcards: 1 if 'H' in wildcards['region'] else 2,
        annotations = lambda wildcards: '' if genotype_alg == 'GnarlyGenotyper' else "-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation -A StrandBiasBySample -A AssemblyComplexity --keep-combined-raw-annotations",
        genotype_alg = genotype_alg
    conda: CONDA_VCF
    resources: 
        mem_mb_java = get_mem_mb_genotype,
        mem_mb=lambda wildcards: 5000 if genotype_alg == 'GnarlyGenotyper' else 10000 #need to make this sample nr. and region level dependent. 
    priority: 40
    shell:"""
        {gatk} {params.genotype_alg} --java-options "-Xmx{resources.mem_mb}M" -R {REF} -V gendb://{input.dir} -O {output.raw_vcfDBI} -D {DBSNP} --intervals {input.intervals} {params.annotations} --annotate-with-num-discovered-alleles --genomicsdb-shared-posixfs-optimizations  --ploidy {params.ploidy} --only-output-calls-starting-in-intervals
        """

rule posterior_phasing:
    input:
       rules.GenotypeDBI.output.raw_vcfDBI
    output:
        rescaled=pj(VCF, "merged_{region}.{genotype_mode}.rescaled.vcf.gz"),
        tbi = ensure(pj(VCF, "merged_{region}.{genotype_mode}.rescaled.vcf.gz.tbi"), non_empty=True)
        # processed=pj(VCF, "merged_{region}.{genotype_mode}.rescaled_phased.vcf.gz")
    params:
        process_phasing='',
        posterior_rescale=''
    conda: CONDA_VCF
    resources:
        mem_mb=1000,
        n=1
    priority: 40
    shell:"""
            gatk CalculateGenotypePosteriors -V {input} -O {output.rescaled}  
        """

