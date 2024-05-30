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
        expand("{current_dir}/{genotype_alg}/{VCF}/ANNOTATED/{region}.annotated.vcf.gz", VCF = VCF, current_dir = current_dir, genotype_alg = genotype_alg, region = parts),
        # rule_all_combine,
        # expand(["{vcf}/Merged_raw_DBI_{chr}.p{chr_p}.{mode}.vcf.gz"],zip,chr=main_chrs_db,chr_p=chr_p, vcf = [config['VCF']]*853, mode = [mode]*853),
        # [f"{genotype_alg}/{VCF}/merged_{region}.{genotype_mode}.vcf.gz" for region in parts],
        # [f"{genotype_alg}/{VCF}/merged_{region}.{genotype_mode}.vcf.gz.tbi" for region in parts],
        # [f"{genotype_alg}/{VCF}/rescaled/{region}_{genotype_mode}.rescaled.vcf.gz" for region in parts],

        # [f"{current_dir}/{genotype_alg}/VCF/ANNOTATED/{region}.annotated.vcf.gz" for region in parts],
        # [f"{current_dir}/{genotype_alg}/{VCF}/ANNOTATED/{region}.annotated.vcf.gz.tbi" for region in parts],
    default_target: True


def get_mem_mb_genotype(wildcrads, attempt):
    return attempt*int(30000)

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
        mem_mb=lambda wildcards: 15000 if genotype_alg == 'GnarlyGenotyper' else 25000 #need to make this sample nr. and region level dependent.
    priority: 40
    shell:"""
        {gatk} {params.genotype_alg} --java-options "-Xmx{resources.mem_mb}M" -R {REF} -V gendb://{input.dir} -O {output.raw_vcfDBI} -D {DBSNP} --intervals {input.intervals} {params.annotations} --annotate-with-num-discovered-alleles --genomicsdb-shared-posixfs-optimizations  --ploidy {params.ploidy} --only-output-calls-starting-in-intervals
        """

rule posterior_phasing:
    input:
       vcf = pj("{genotype_alg}",VCF, "merged_{region}.{genotype_mode}.vcf.gz")
    output:
        rescaled=pj("{genotype_alg}", VCF, "rescaled", "{region}_{genotype_mode}.rescaled.vcf.gz"),
        tbi = ensure(pj("{genotype_alg}", VCF, "rescaled", "{region}_{genotype_mode}.rescaled.vcf.gz.tbi"), non_empty=True)
    conda: CONDA_VCF
    resources:
        mem_mb=15000,
        n=2
    priority: 40
    shell:"""
            mkdir -p {genotype_alg}/vcfs/rescaled
            gatk CalculateGenotypePosteriors --java-options "-Xmx{resources.mem_mb}M" -V {input} -O {output.rescaled}  
        """



rule extract_positions:
    input: vcf = pj("{genotype_alg}", VCF, "rescaled", "{region}_{genotype_mode}.rescaled.vcf.gz"),
    output: vcf = temp(pj(current_dir, "{genotype_alg}", VCF,  "ANNOTATED_temp" , "{region}_pos_only.vcf"))
    conda: CONDA_MAIN
    shell: "bcftools view --drop-genotypes -O v -o {output.vcf} {input.vcf}"

rule annotate_revel:
    input: vcf = temp(pj(current_dir, "{genotype_alg}", VCF,  "ANNOTATED_temp" , "{region}_pos_only.vcf"))
    output: vcf_annotated = pj(current_dir, "{genotype_alg}", VCF,  "ANNOTATED_temp" , "{region}.annotated_pos_only.vcf")
    conda: CONDA_MAIN
    resources: n = "4",
            mem_mb = 6000
    # log: pj(current_dir,"logs","glnexus","annotate_revel_{region}.{genotype_mode}.{types_of_gl}.log")
    params: temp_dir = pj(current_dir, "{genotype_alg}", VCF,  "ANNOTATED_temp")
    shell:
        """
        mkdir -p {params.temp_dir} &&

        bcftools annotate -a {REVEL} -h {REVEL_header} -c CHROM,POS,REF,ALT,REVEL {input.vcf} -O v -o {output.vcf_annotated} --threads {resources.n} 2> {log}
        """

rule annotate_gene:
    input: temp_vcf = pj(current_dir, "{genotype_alg}", VCF,  "ANNOTATED_temp" , "{region}.annotated_pos_only.vcf")
    output: vcf_annotated=temp(pj(current_dir, "{genotype_alg}", VCF,  "ANNOTATED_temp" ,"{region}.annotated.hg38_multianno.vcf")),
    conda: CONDA_ANNOVAR
    params:
        out=pj(current_dir, "{genotype_alg}", VCF,  "ANNOTATED_temp" ,"{region}.annotated")
    # log: pj(current_dir,"logs","glnexus","annotate_gene_{region}.{genotype_mode}.{types_of_gl}.log")
    resources: n = "2",
                mem_mb = 5000
    shell:
        """
        perl {annovar} {input.temp_vcf} {annovar_db} -out {params.out} -protocol ensGene,refGene -operation g,g -vcfinput -buildver hg38 -thread {resources.n} 2> {log}
        """

rule bring_anno_to_samples:
    input: vcf_annotated = pj(current_dir, "{genotype_alg}", VCF,  "ANNOTATED_temp" ,"{region}.annotated.hg38_multianno.vcf"),
            samples_vcf =  pj("{genotype_alg}", VCF, "rescaled", "{region}_{genotype_mode}.rescaled.vcf.gz"),
    output: vcf_anno_samples = pj(current_dir, "{genotype_alg}", VCF,  "ANNOTATED" ,"{region}.annotated.vcf.gz"),
            # tbi = ensure(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz.tbi"), non_empty=True)
    conda: CONDA_MAIN
    resources:
        n = "2",
        mem_mb = 6000
    shell:
        """
        bgzip {input.vcf_annotated}
        tabix -p vcf {input.vcf_annotated}.gz
        bcftools annotate -a {input.vcf_annotated}.gz -c INFO -O z -o {output.vcf_anno_samples} {input.samples_vcf}
        tabix -p vcf {output.vcf_anno_samples}  
        """

