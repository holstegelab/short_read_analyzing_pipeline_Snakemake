import os
from itertools import repeat
wildcard_constraints:
    sample="[\w\d_\-@]+",
    genotype_level="[2|3|4]",
    region="[0-9XYAFOH]+",
    genotype_alg="[GenotypeGVCFs|GnarlyGenotyper]",
    genotype_mode="[WES|WGS]"

from common import *

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()



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
    raise RuntimeError(f'Unknown level {genotype_level}')


rule Genotype_all:
    input:
        # rule_all_combine,
        # expand(["{vcf}/Merged_raw_DBI_{chr}.p{chr_p}.{mode}.vcf.gz"],zip,chr=main_chrs_db,chr_p=chr_p, vcf = [config['VCF']]*853, mode = [mode]*853),
        [f"{genotype_alg}/{VCF}/merged_{region}.{genotype_mode}.vcf.gz" for region in parts],
        [f"{genotype_alg}/{VCF}/merged_{region}.{genotype_mode}.vcf.gz.tbi" for region in parts],
        [f"{genotype_alg}/{VCF}/rescaled/{region}_{genotype_mode}.rescaled.vcf.gz" for region in parts],
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
       vcf = expand(pj("{genotype_alg}",VCF, "merged_{region}.{genotype_mode}.vcf.gz"), genotype_alg = genotype_alg, genotype_mode = genotype_mode, allow_missing=True),
    output:
        rescaled=expand(pj("{genotype_alg}", VCF, "rescaled", "{region}_{genotype_mode}.rescaled.vcf.gz"), genotype_alg = genotype_alg, genotype_mode = genotype_mode, allow_missing=True),
        tbi = ensure(expand(pj("{genotype_alg}", VCF, "rescaled", "{region}_{genotype_mode}.rescaled.vcf.gz.tbi"), genotype_alg = genotype_alg, genotype_mode = genotype_mode, allow_missing = True), non_empty=True)
    conda: CONDA_VCF
    resources:
        mem_mb=5000,
        n=1
    priority: 40
    shell:"""
            mkdir -p {genotype_alg}/vcfs/rescaled
            gatk CalculateGenotypePosteriors --java-options "-Xmx{resources.mem_mb}M" -V {input} -O {output.rescaled}  
        """

