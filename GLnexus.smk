import os
from common import *
current_dir = os.getcwd()


wildcard_constraints:
    sample="[\w\d_\-@]+",
    chr = "[\w\d]+",
    # readgroup="[\w\d_\-@]+"

gvcf_caller = config.get("caller", "BOTH")
glnexus_filtration = config.get("glnexus_filtration", "custom")
sample_types = config.get("sample_types","WES")


print(f"Caller: {gvcf_caller}")
print(f"Sample type: {sample_types}")
print(f"Filtration setting: {glnexus_filtration}")
parts = level2_regions_diploid

def generate_gvcf_input(gvcf_folder):
    res = []
    for samplefile in SAMPLE_FILES:
        sample_names = SAMPLEFILE_TO_SAMPLES[samplefile]
        samplefile_folder = get_samplefile_folder(samplefile)
        gvcf_input = []
        for sample in sample_names:
            # Determine if it is WGS or WES
            if SAMPLEINFO[sample]["sample_type"] == "WGS":
                # Check the sex of the sample
                sex_file = pj(samplefile_folder, KMER, SAMPLEINFO[sample]["sample"] + ".result.yaml")
                with open(sex_file) as f:
                    xsample = yaml.load(f,Loader=yaml.FullLoader)
                    if  xsample['sex'] == 'M' or not parts.startswith('Y'):
                        region = convert_to_level1(parts)
                    else:
                        continue
            else:  # WES
                sex_file = pj(samplefile_folder, KMER,SAMPLEINFO[sample]["sample"] + ".result.yaml")
                with open(sex_file) as f:
                    xsample = yaml.load(f,Loader=yaml.FullLoader)
                    if xsample['sex'] == 'M' or not parts.startswith('Y'):
                        region = convert_to_level0(parts)
                    else:
                        continue
            filename = expand("{cd}/{GVCF}/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=gvcf_folder,region = region, sample=sample_names,allow_missing=True)
            gvcf_input.append(filename)
        res.extend(gvcf_input)
    return res

if gvcf_caller == "HaplotypeCaller":
    gvcf_input = generate_gvcf_input(GVCF + '/exome_gatk')
    glnexus_dir = ["GLnexus_on_Haplotypecaller"]

elif gvcf_caller == "Deepvariant":
    gvcf_input = generate_gvcf_input(DEEPVARIANT + '/gVCF/exomes')
    glnexus_dir = ["GLnexus_on_Deepvariant"]

elif gvcf_caller == "BOTH":
    gvcf_input = generate_gvcf_input(GVCF + '/exome_gatk')
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
        "Invalid option provided to 'glnexus_filtration'; \n"
        "please choose either 'default' for default GLnexus filtration \n "
        "or 'custom' (default) for absence of hard filters \n"
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
    return region_to_file(region,wgs=sample_types == 'WGS',extension='bed')



rule GLnexus_all:
    input:
        expand("{cur_dir}/{types_of_gl}{appendix}/{region}/{parts}.vcf.gz", cur_dir = current_dir, region = ["F"], parts = parts, types_of_gl = glnexus_dir, appendix = dir_appendix),
        expand("{cur_dir}/{types_of_gl}{appendix}/{region}/{parts}.vcf.gz.tbi", cur_dir = current_dir, region = ["F"], parts = parts, types_of_gl = glnexus_dir, appendix = dir_appendix),
    default_target: True

rule glnexus:
    input: gvcf_input
    output: vcf = pj(current_dir, glnexus_dir[0] + dir_appendix, "{region}", "{parts}.vcf.gz")
    container: "docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    params: bed = region_to_bed_file,
            mem_gb = 7,
            scratch_dir =  temp(current_dir + '/' + tmpdir + "/{region}_{parts}_glnexus.DB"),
            conf_filters = conf_filter
    threads: 4
    resources: 
        n = "5",
        mem_mb = 7000
    shell:
        """
        rm -rf {params.scratch_dir} &&
        glnexus_cli  --dir {params.scratch_dir} --bed {params.bed} --threads {threads} --mem-gbytes {params.mem_gb} --config {params.conf_filters}  {input}  | bcftools view -  | bgzip -@ {threads} -c > {output} 
        """
rule index_deep:
    input: rules.glnexus.output.vcf
    output: tbi = pj(current_dir, glnexus_dir[0] + dir_appendix, "{region}","{parts}.vcf.gz.tbi")
    conda: CONDA_VCF
    shell: "gatk IndexFeatureFile -I {input}"


if gvcf_caller == "BOTH":
    use rule glnexus as glnexus_2 with:
        input: gvcf_input = generate_gvcf_input(DEEPVARIANT + '/gVCF/exomes')
        output: vcf=pj(current_dir,glnexus_dir[1] + dir_appendix,"{region}", "{parts}.vcf.gz")
        params: scratch_dir =  temp(current_dir + '/' + tmpdir + "/{region}_{parts}_glnexus_2.DB"),
                bed= region_to_bed_file,
                mem_gb= 7,
                conf_filters= conf_filter
    use rule index_deep as index_deep_2 with:
        input: rules.glnexus_2.output.vcf
        output: tbi = pj(current_dir, glnexus_dir[1] + dir_appendix, "{region}","{parts}.vcf.gz.tbi")
