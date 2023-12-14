import pandas as pd
import read_stats
import os
import getpass
import read_samples
from common import *
import utils
current_dir = os.getcwd()
onsuccess: shell("rm -fr logs/GLnexus/*")
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

gvcf_caller = config.get("caller", "HaplotypeCaller")
glnexus_filtration = config.get("glnexus_filtration", "custom")
sample_types = config.get("sample_types","WES")


if gvcf_caller == "HaplotypeCaller":
    use rule * from gVCF
    rule_gvcf_all_input = rules.gVCF_all.input
    gvcf_input = expand("{cd}/{GVCF}/reblock/{region}/{sample}.{region}.wg.vcf.gz",cd = current_dir, GVCF = GVCF, sample=sample_names,allow_missing=True),
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
    gvcf_input = expand("{cd}/{GVCF}/reblock/{region}/{sample}.{region}.wg.vcf.gz",cd=current_dir,GVCF=GVCF,sample=sample_names,allow_missing=True),
    glnexus_dir = ["GLnexus_on_Haplotypecaller", "GLnexus_on_Deepvariant"]

else:
    raise ValueError(
        "invalid option provided to 'caller'; please choose either 'HaplotypeCaller' or 'Deepvariant'."
    )

if glnexus_filtration == 'default':
    dir_appendix = "default"
elif glnexus_filtration == 'custom':
    dir_appendix = "custom"
else:
    raise ValueError(
        "Invalid option provided to 'glnexus_filtration';"
        "please choose either 'default' or 'custom'"
        "custom preset located at /gpfs/work3/0/qtholstg/hg38_res_v2/software/Glnexus_preset.yml"
    )

def conf_filter(wildcards):
    if glnexus_filtration == 'default':
        if gvcf_caller == "Deepvariant":
            if sample_types == 'WES':
                conf_filters = "DeepVariantWES"
            else:
                conf_filters = "DeepVariantWGS"
        elif gvcf_caller == "HaplotypeCaller":
            conf_filters = "gatk"
    elif glnexus_filtration == 'custom':
        conf_filters = "/gpfs/work3/0/qtholstg/hg38_res_v2/software/Glnexus_preset.yml"
    return conf_filters

def region_to_bed_file(wildcards):#{{{
    """Converts a region to a bed file location (see common.py and Tools.smk)"""
    region = wildcards['parts']
    if sample_types == 'WES':
        return region_to_file(region,extension='bed')
    elif sample_types == 'WGS':
        return region_to_file(region, wgs=True, extension='bed')#}}}

if sample_types == 'WES':
    parts = level2_regions_diploid
else:
    parts = level3_regions_diploid


rule GLnexus_all:
    input:
        expand("{cur_dir}/{types_of_gl}{appendix}/{region}/{parts}.vcf.gz", cur_dir = current_dir, region = "F", parts = parts, types_of_gl = glnexus_dir, appendix = dir_appendix),
        expand("{cur_dir}/{types_of_gl}{appendix}/{region}/{parts}.vcf.gz.tbi", cur_dir = current_dir, region = "F", parts = parts, types_of_gl = glnexus_dir, appendix = dir_appendix),
        rule_gvcf_all_input
    default_target: True

rule glnexus:
    input: gvcf_input
    output: vcf = pj(current_dir, glnexus_dir[0] + dir_appendix, "{region}", "{parts}.vcf.gz")
    container: "docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    params: bed = region_to_bed_file,
            mem_gb = 7,
            scratch_dir =  temp(current_dir + '/' + tmpdir + "/{region}_{parts}_glnexus.DB"),
            conf_filters = conf_filter
    log: pj(current_dir,LOG,"GLnexus","{region}.{parts}.glnexus.log")
    threads: 4
    resources: mem_mb = 7000
    shell:
        """
        rm -rf {params.scratch_dir} &&
        glnexus_cli  --dir {params.scratch_dir} --bed {params.bed} --threads {threads} --mem-gbytes {params.mem_gb} --config {params.conf_filters}  {input} 2> {log} | bcftools view -  2>> {log}| bgzip -@ {threads} -c > {output} 2>> {log}
        """
rule index_deep:
    input: rules.glnexus.output.vcf
    output: tbi = pj(current_dir, glnexus_dir[0] + dir_appendix, "{region}","{parts}.vcf.gz.tbi")
    conda: CONDA_VCF
    shell: "gatk IndexFeatureFile -I {input}"


if gvcf_caller == "BOTH":
    use rule glnexus as glnexus_2 with:
        input: gvcf_input = expand("{cd}/{DEEPVARIANT}/gVCF/{region}/{sample}.{region}.wg.vcf.gz", cd = current_dir, DEEPVARIANT = DEEPVARIANT, sample = sample_names, allow_missing=True)
        output: vcf=pj(current_dir,glnexus_dir[1] + dir_appendix,"{region}", "{parts}.vcf.gz")
        params: scratch_dir =  temp(current_dir + '/' + tmpdir + "/{region}_{parts}_glnexus_2.DB"),
                bed= region_to_bed_file,
                mem_gb= 7,
                conf_filters= conf_filter
        log: pj(current_dir,LOG,"GLnexus","{region}.{parts}.glnexus_2.log")
    use rule index_deep as index_deep_2 with:
        input: rules.glnexus_2.output.vcf
        output: tbi = pj(current_dir, glnexus_dir[1] + dir_appendix, "{region}","{parts}.vcf.gz.tbi")
