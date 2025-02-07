import os
from itertools import repeat
wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    genotype_level=r"[2|3|4]",
    region=r"[0-9XYAFOH]+",
    genotype_alg=r"(GenotypeGVCFs|GnarlyGenotyper)",
    genotype_mode=r"(WES|WGS)"

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
        [f"{genotype_alg}/{VCF}/merged_{region}.{genotype_mode}.vcf.gz" for region in parts],
        expand(['{genotype_alg}/{VCF}/ANNOTATED/{region}_{genotype_mode}.annotated.vcf.gz.tbi'], genotype_mode = genotype_mode, VCF = VCF, genotype_alg = genotype_alg, region = parts),
     # pj("{genotype_alg}", VCF, "ANNOTATED", "{region}_{genotype_mode}.annotated.vcf.gz")
        # rule_all_combine,
        # expand(["{vcf}/Merged_raw_DBI_{chr}.p{chr_p}.{mode}.vcf.gz"],zip,chr=main_chrs_db,chr_p=chr_p, vcf = [config['VCF']]*853, mode = [mode]*853),

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
        raw_vcfDBI=ensure(pj("{genotype_alg}",VCF, "merged_{region}.{genotype_mode}.vcf.gz"), non_empty=True),
        tbi = ensure(pj("{genotype_alg}",VCF, "merged_{region}.{genotype_mode}.vcf.gz.tbi"), non_empty=True),
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
    run:
        if wildcards['region'].endswith("H"):
            shell("mkdir -p {genotype_alg}/vcfs/rescaled")
            shell("ln -r {input.vcf} {output.rescaled}")
            shell("ln -r {input.vcf}.tbi {output.tbi}")
        else:
            shell("mkdir -p {genotype_alg}/vcfs/rescaled")
            shell("gatk CalculateGenotypePosteriors --java-options \"-Xmx{resources.mem_mb}M\" -V {input.vcf} -O {output.rescaled}")
    #
    # shell:
    #     """
    #         mkdir -p {genotype_alg}/vcfs/rescaled
    #         gatk CalculateGenotypePosteriors --java-options "-Xmx{resources.mem_mb}M" -V {input} -O {output.rescaled}
    #     """
    #

rule split_multiallelic:
    input: vcf=pj("{genotype_alg}",VCF,"rescaled","{region}_{genotype_mode}.rescaled.vcf.gz"),
    output: temp_tsv = temp(pj("{genotype_alg}",VCF,"rescaled","{region}_{genotype_mode}.rescaled.split.tsv.gz")),
            vcf=pj("{genotype_alg}",VCF,"rescaled","{region}_{genotype_mode}.rescaled.split.vcf.gz"),
    conda: CONDA_MAIN
    resources:
        mem_mb=24000,
        n=12
    shell:
        """
        {bcftools_patched} view -U -c1:nonmajor {input.vcf}  --threads {resources.n} -Ou | {bcftools_patched} query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %AS_FS]\t%AC\t%AN\n' | awk '{{l2=split($4,a,","); if(l2 > 1 && $7 > 0) {{split($6,ac,","); ref=1.0; for (i in ac) ref-= (ac[i]/$7); failcount=0; for (i in ac) failcount+=(ref<(ac[i]/$7)); OFS="\t"; print $1, $2, $3, $4, $4, $5, l2, $6, ref, failcount}}}}' | bgzip --threads {resources.n}  > {output.temp_tsv}
        tabix -s1 -b2 -e2 {output.temp_tsv}
        {bcftools_patched} view -U -c1:nonmajor {input.vcf} --threads {resources.n} -Ou | {bcftools_patched} annotate -a {output.temp_tsv} -h {multiallelic_hdr} -c CHROM,POS,REF,ALT,MA_ALT,MA_FILTER,NMA_ALT,MA_ALT_AC,MA_REF_AF,MA_REFLOW_COUNT --threads {resources.n} -Ov | {bcftools_patched} norm -f {REF} -m- -c w --thread {resources.n} -o {output.vcf} --output-type z
        """


rule extract_positions:
    input: vcf = pj("{genotype_alg}", VCF, "rescaled", "{region}_{genotype_mode}.rescaled.split.vcf.gz"),
    output: vcf = temp(pj("{genotype_alg}", VCF,  "ANNOTATED_temp" , "{region}_{genotype_mode}_pos_only.vcf"))
    conda: CONDA_MAIN
    shell: "bcftools view --drop-genotypes -O v -o {output.vcf} {input.vcf}"

rule annotate_revel:
    input: vcf = pj("{genotype_alg}", VCF,  "ANNOTATED_temp" , "{region}_{genotype_mode}_pos_only.vcf")
    output: vcf_annotated = pj("{genotype_alg}", VCF,  "ANNOTATED_temp" , "{region}_{genotype_mode}.annotated_pos_only.vcf")
    conda: CONDA_MAIN
    resources: n = "4",
            mem_mb = 6000
    # log: pj(current_dir,"logs","glnexus","annotate_revel_{region}.{genotype_mode}.{types_of_gl}.log")
    params: temp_dir = pj("{genotype_alg}", VCF,  "ANNOTATED_temp")
    shell:
        """
        mkdir -p {params.temp_dir} &&

        bcftools annotate -a {REVEL} -h {REVEL_header} -c CHROM,POS,REF,ALT,REVEL {input.vcf} -O v -o {output.vcf_annotated} --threads {resources.n}
        """

rule annotate_gene:
    input: temp_vcf = pj("{genotype_alg}", VCF,  "ANNOTATED_temp" , "{region}_{genotype_mode}.annotated_pos_only.vcf")
    output: vcf_annotated=pj("{genotype_alg}", VCF,  "ANNOTATED_temp" ,"{region}_{genotype_mode}.annotated.hg38_multianno.vcf"),
    conda: CONDA_ANNOVAR
    params:
        out=pj("{genotype_alg}", VCF,  "ANNOTATED_temp" ,"{region}_{genotype_mode}.annotated")
    # log: pj(current_dir,"logs","glnexus","annotate_gene_{region}.{genotype_mode}.{types_of_gl}.log")
    resources: n = "2",
                mem_mb = 5000
    shell:
        """
        perl {annovar} {input.temp_vcf} {annovar_db} -out {params.out} -protocol ensGene,refGene -operation g,g -vcfinput -buildver hg38 -thread {resources.n} 
        """

rule bring_anno_to_samples:
    input: vcf_annotated = pj("{genotype_alg}", VCF,  "ANNOTATED_temp" ,"{region}_{genotype_mode}.annotated.hg38_multianno.vcf"),
            samples_vcf =  pj("{genotype_alg}", VCF, "rescaled", "{region}_{genotype_mode}.rescaled.vcf.gz"),
    output:
            vcf_anno_samples = "{genotype_alg}/{VCF}/ANNOTATED/{region}_{genotype_mode}.annotated.vcf.gz",
            # vcf_anno_samples = pj("{genotype_alg}", VCF,  "ANNOTATED" ,"{region}_{genotype_mode}.annotated.vcf.gz"),
            # tbi = ensure(pj("{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz.tbi"), non_empty=True)
    conda: CONDA_MAIN
    resources:
        n = "2",
        mem_mb = 6000
    shell:
        """
        bgzip {input.vcf_annotated}
        tabix -p vcf {input.vcf_annotated}.gz
        bcftools annotate -a {input.vcf_annotated}.gz -c INFO -O z -o {output.vcf_anno_samples} {input.samples_vcf}
        """
        # tabix -p vcf {output.vcf_anno_samples}


rule checking:
    input: vcf = "{genotype_alg}/{VCF}/ANNOTATED/{region}_{genotype_mode}.annotated.vcf.gz"
    output: vcf = "{genotype_alg}/{VCF}/ANNOTATED/{region}_{genotype_mode}.annotated.vcf.gz.tbi"
    conda: CONDA_MAIN
    shell:
        """
        
        tabix -p vcf {input.vcf}
        rm -rf {genotype_alg}/{VCF}/ANNOTATED_temp/{region}_{genotype_mode}*
        """