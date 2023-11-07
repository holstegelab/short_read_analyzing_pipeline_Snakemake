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
    chr = "[\w\d]+",
    # readgroup="[\w\d_\-@]+"

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner
module Deepvariant:
    snakefile: 'Deepvariant.smk'
    config: config
module gVCF:
    snakefile: 'gVCF.smk'
    config: config
# use rule * from Deepvariant

gvcf_caller = config.get("caller", "HaplotypeCaller")
glnexus_filtration = config.get("glnexus_filtration", "default")



if gvcf_caller == "HaplotypeCaller":
    use rule * from gVCF
    rule_gvcf_all_input = rules.gVCF_all.input
    gvcf_input = expand("{cd}/{GVCF}/reblock/{region}/{sample}.{region}.g.vcf.gz",cd = current_dir, GVCF = GVCF, sample=sample_names,allow_missing=True),
    glnexus_dir = ["GLnexus_on_Haplotypecaller"]

elif gvcf_caller == "Deepvariant":
    use rule * from Deepvariant
    rule_gvcf_all_input = rules.DeepVariant_all.input
    gvcf_input = expand("{cd}/{DEEPVARIANT}/gVCF/{region}/{sample}.{region}.wg.vcf.gz", cd = current_dir, DEEPVARIANT = DEEPVARIANT, sample = sample_names, allow_missing=True)
    glnexus_dir = ["GLnexus_on_Deepvariant"]

elif gvcf_caller == "BOTH":
    use rule * from gVCF
    use rule * from Deepvariant
    rule_gvcf_all_input = [rules.gVCF_all.input, rules.DeepVariant_all.input]
    gvcf_input = expand("{cd}/{GVCF}/reblock/{region}/{sample}.{region}.g.vcf.gz",cd=current_dir,GVCF=GVCF,sample=sample_names,allow_missing=True),
    glnexus_dir = ["GLnexus_on_Haplotypecaller", "GLnexus_on_Deepvariant"]

else:
    raise ValueError(
        "invalid option provided to 'caller'; please choose either 'HaplotypeCaller' or 'Deepvariant'."
    )

if glnexus_filtration == 'default':
    dir_appendix = "default"
elif glnexus_filtration == 'custom':
    dir_appendix = "custom"

def conf_filter(wildcards):
    if glnexus_filtration == 'default':
        if gvcf_caller == "Deepvariant":
            if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
                conf_filters = "DeepVariantWGS"
            else:
                conf_filters = "DeepVariantWES"
        elif gvcf_caller == "HaplotypeCaller":
            conf_filters = "gatk"
        # dir_appendix = "default"
    elif glnexus_filtration == 'custom':
        conf_filters = "/gpfs/home1/gozhegov/short_read_analyzing_pipeline_Snakemake/Glnexus_preset.yml"
        # dir_appendix = "custom"
    return conf_filters

# def get_chrom_capture_kit_bed(wildcards):
#     chr = wildcards.chr
#     if mode == 'WES':
#         capture_kit_chr_path_bed = pj(config['RES'], config['kit_folder'], config['MERGED_CAPTURE_KIT'] + '_hg38', config['MERGED_CAPTURE_KIT'] + '_hg38_' + chr + '.interval_list.bed')
#     elif mode == 'WGS':
#         capture_kit_chr_path_bed = chr
#     return capture_kit_chr_path_bed

def region_to_bed_file(wildcards):#{{{
    """Converts a region to a bed file location (see common.py and Tools.smk)"""
    sample = wildcards['sample']
    region = wildcards['region']
    return region_to_file(region, wgs='wgs' in SAMPLEINFO[sample]['sample_type'], extension='bed')#}}}


rule GLnexus_all:
    input:
        expand("{cur_dir}/{types_of_gl}{appendix}/{region}.vcf.gz", cur_dir = current_dir, region = level1_regions, types_of_gl = glnexus_dir, appendix = dir_appendix),
        expand("{cur_dir}/{types_of_gl}{appendix}/{region}.vcf.gz.tbi", cur_dir = current_dir, region = level1_regions, types_of_gl = glnexus_dir, appendix = dir_appendix),
        rule_gvcf_all_input
    default_target: True

rule glnexus:
    input: gvcf_input
    output: vcf = pj(current_dir, glnexus_dir[0] + dir_appendix, "{region}.vcf.gz")
    container: "docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    params: bed = region_to_bed_file,
            mem_gb = 7,
            scratch_dir =  temp(current_dir + '/' + tmpdir + "/{region}_glnexus.DB"),
            conf_filters = conf_filter
    log: pj(current_dir,LOG,"{region}.glnexus.log")
    benchmark:
        pj(current_dir,BENCH,"{region}.glnexus.txt")
    threads: 4
    resources: mem_mb = 7000
    shell:
        """
        glnexus_cli  --dir {params.scratch_dir} --bed {params.bed} --threads {threads} --mem-gbytes {params.mem_gb} --config {params.conf_filters}  {input} | bcftools view - | bgzip -@ {threads} -c > {output} 2> {log}
        """
rule index_deep:
    input: rules.glnexus.output.vcf
    output: tbi = pj(current_dir, glnexus_dir[0] + dir_appendix, "{region}.vcf.gz.tbi")
    conda: "envs/preprocess.yaml"
    shell: "gatk IndexFeatureFile -I {input}"

if gvcf_caller == "BOTH":
    use rule glnexus as glnexus_2 with:
        input: gvcf_input = expand("{cd}/{DEEPVARIANT}/gVCF/{region}/{sample}.{region}.wg.vcf.gz", cd = current_dir, DEEPVARIANT = DEEPVARIANT, sample = sample_names, allow_missing=True)
        output: vcf=pj(current_dir,glnexus_dir[1] + dir_appendix,"{region}.vcf.gz")
        params: scratch_dir =  temp(current_dir + '/' + tmpdir + "/{region}_glnexus_2.DB"),
                bed= region_to_bed_file,
                mem_gb= 7,
                conf_filters= conf_filter
        benchmark:
            pj(current_dir,BENCH,"{region}.glnexus_2.txt")
        log: pj(current_dir,LOG,"{region}.glnexus_2.log")
    use rule index_deep as index_deep_2 with:
        input: rules.glnexus_2.output.vcf
        output: tbi = pj(current_dir, glnexus_dir[1] + dir_appendix, "{region}.vcf.gz.tbi")

# rule norma_gln:
#     input: rules.glnexus.output.vcf
#     output: normVCF = pj(current_dir, "glnexus_norm", "{chr}", "{chr}.vcf.gz")
#     log: pj(current_dir,config['LOG'],"{chr}.norma_gln.log")
#     conda: "envs/preprocess.yaml"
#     benchmark:
#         pj(current_dir,config['BENCH'],"{chr}.norma_gln.txt")
#     shell: "bcftools norm -f {ref} {input} -m -both -O v | bcftools norm --check-ref ws -d exact -f {ref} -O z > {output.normVCF} 2> {log}"
