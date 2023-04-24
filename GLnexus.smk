import pandas as pd
import read_stats
import os
import getpass

configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")
gatk = config['gatk']
samtools = config['samtools']
bcftools = config['bcftools']
dragmap = config['dragmap']
verifybamid2 = config['verifybamid2']

ref = config['RES'] + config['ref']
current_dir = os.getcwd()
tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())
tmpdir_alternative = os.path.join(config['tmpdir'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

wildcard_constraints:
    sample="[\w\d_\-@]+",
    chr = "[\w\d]+",
    # readgroup="[\w\d_\-@]+"


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO = load_samplefiles('.', config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

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

mode = config.get("computing_mode", "WES")


gvcf_caller = config.get("caller", "HaplotypeCaller")
glnexus_filtration = config.get("glnexus_filtration", "default")



if gvcf_caller == "HaplotypeCaller":
    use rule * from gVCF
    rule_gvcf_all_input = rules.gVCF_all.input
    gvcf_input = expand("{cd}/{gvcfs}/{chr}/{sample}.{chr}.{sex}.g.vcf.gz",cd = current_dir, gvcfs=config['gVCF'],sample=sample_names,mode=mode,allow_missing=True),
    glnexus_dir = ["GLnexus_on_Haplotypecaller"]

elif gvcf_caller == "Deepvariant":
    use rule * from Deepvariant
    rule_gvcf_all_input = rules.Deepvariant_all.input
    gvcf_input = expand("{cd}/{dp}/gVCF/{chr}/{chr}.{sample}.{sex}.g.vcf.gz", cd = current_dir, dp = config['DEEPVARIANT'], sample = sample_names, mode = mode, allow_missing=True)
    glnexus_dir = ["GLnexus_on_Deepvariant"]

elif gvcf_caller == "BOTH":
    use rule * from gVCF
    use rule * from Deepvariant
    rule_gvcf_all_input = [rules.gVCF_all.input, rules.Deepvariant_all.input]
    gvcf_input = expand("{cd}/{gvcfs}/{chr}/{sample}.{chr}.{sex}.g.vcf.gz",cd = current_dir, gvcfs=config['gVCF'],sample=sample_names,mode=[mode],allow_missing=True),
    glnexus_dir = ["GLnexus_on_Haplotypecaller", "GLnexus_on_Deepvariant"]

else:
    raise ValueError(
        "invalid option provided to 'caller'; please choose either 'HaplotypeCaller' or 'Deepvariant'."
    )

if glnexus_filtration == 'default':
    # conf_filters = "DeepVariant{{wildcards.mode}}"
    dir_appendix = "default"
elif glnexus_filtration == 'custom':
    # conf_filters = "/gpfs/home1/gozhegov/short_read_analyzing_pipeline_Snakemake/Glnexus_preset.yml"
    dir_appendix = "custom"

def conf_filter(wildcards):
    if glnexus_filtration == 'default':
        if gvcf_caller == "Deepvariant":
            conf_filters = "DeepVariant" + wildcards.mode
        elif gvcf_caller == "HaplotypeCaller":
            conf_filters = "gatk"
        # dir_appendix = "default"
    elif glnexus_filtration == 'custom':
        conf_filters = "/gpfs/home1/gozhegov/short_read_analyzing_pipeline_Snakemake/Glnexus_preset.yml"
        # dir_appendix = "custom"
    return conf_filters

def get_chrom_capture_kit_bed(wildcards):
    chr = wildcards.chr
    if mode == 'WES':
        capture_kit_chr_path_bed = os.path.join(config['RES'], config['kit_folder'], config['MERGED_CAPTURE_KIT'] + '_hg38', config['MERGED_CAPTURE_KIT'] + '_hg38_' + chr + '.interval_list.bed')
    elif mode == 'WGS':
        capture_kit_chr_path_bed = chr
    return capture_kit_chr_path_bed

rule GLnexus_all:
    input:
        expand("{cur_dir}/{types_of_gl}{appendix}/{chr}/{chr}_{sex}.vcf.gz", cur_dir = current_dir, mode = mode, chr = main_chrs, types_of_gl = glnexus_dir, appendix = dir_appendix),
        expand("{cur_dir}/{types_of_gl}{appendix}/{chr}/{chr}_{sex}.vcf.gz.tbi", cur_dir = current_dir, mode = mode, chr = main_chrs, types_of_gl = glnexus_dir, appendix = dir_appendix),
        rule_gvcf_all_input
    default_target: True

rule glnexus:
    input: gvcf_input
    # input: gvcf = expand("{cd}/{dp}/gVCF/{chr}.{sample}.{sex}.g.vcf.gz", cd = current_dir, dp = config['DEEPVARIANT'], sample = sample_names, mode = mode, allow_missing=True)
    output: vcf = os.path.join(current_dir, glnexus_dir[0] + dir_appendix, "{chr}", "{chr}_{sex}.vcf.gz")
    container: "docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    params: bed = get_chrom_capture_kit_bed,
            mem_gb = int(config['glnexus']['mem'] // 1024),
            scratch_dir =  temp(current_dir + '/' + tmpdir + "/{chr}_{sex}_glnexus.DB"),
            conf_filters = conf_filter
    log: os.path.join(current_dir,config['LOG'],"{chr}_{sex}.glnexus.log")
    benchmark:
        os.path.join(current_dir,config['BENCH'],"{chr}_{sex}.glnexus.txt")
    threads: config['glnexus']['n']
    resources: mem_mb = config['glnexus']['mem']
    shell:
        """
        glnexus_cli  --dir {params.scratch_dir} --bed {params.bed} --threads {threads} --mem-gbytes {params.mem_gb} --config {params.conf_filters}  {input} | bcftools view - | bgzip -@ {threads} -c > {output} 2> {log}
        """
rule index_deep:
    input: rules.glnexus.output.vcf
    output: tbi = os.path.join(current_dir, glnexus_dir[0] + dir_appendix, "{chr}", "{chr}_{sex}.vcf.gz.tbi")
    conda: "envs/preprocess.yaml"
    shell: "gatk IndexFeatureFile -I {input}"

if gvcf_caller == "BOTH":
    use rule glnexus as glnexus_2 with:
        input: expand("{cd}/{dp}/gVCF/{chr}/{chr}.{sample}.{sex}.g.vcf.gz", cd = current_dir, dp = config['DEEPVARIANT'], sample = sample_names, mode = mode, allow_missing=True)
        output: vcf=os.path.join(current_dir,glnexus_dir[1] + dir_appendix,"{chr}","{chr}_{sex}.vcf.gz")
        params: scratch_dir =  temp(current_dir + '/' + tmpdir + "/{chr}_{sex}_glnexus_2.DB"),
                bed= get_chrom_capture_kit_bed,
                mem_gb= int(config['glnexus']['mem'] // 1024),
                conf_filters= conf_filter
        benchmark:
            os.path.join(current_dir,config['BENCH'],"{chr}_{sex}.glnexus_2.txt")
    use rule index_deep as index_deep_2 with:
        input: rules.glnexus_2.output.vcf
        output: tbi = os.path.join(current_dir, glnexus_dir[1] + dir_appendix, "{chr}", "{chr}_{sex}.vcf.gz.tbi")

# rule norma_gln:
#     input: rules.glnexus.output.vcf
#     output: normVCF = os.path.join(current_dir, "glnexus_norm", "{chr}", "{chr}_{sex}.vcf.gz")
#     log: os.path.join(current_dir,config['LOG'],"{chr}_{sex}.norma_gln.log")
#     conda: "envs/preprocess.yaml"
#     benchmark:
#         os.path.join(current_dir,config['BENCH'],"{chr}_{sex}.norma_gln.txt")
#     shell: "bcftools norm -f {ref} {input} -m -both -O v | bcftools norm --check-ref ws -d exact -f {ref} -O z > {output.normVCF} 2> {log}"
