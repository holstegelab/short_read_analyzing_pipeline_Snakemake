import os
from common import *
current_dir = os.getcwd()


wildcard_constraints:
    sample="[\w\d_\-@]+",
    chr = "[\w\d]+",
    region = "[\w\d]+",
    # readgroup="[\w\d_\-@]+"

gvcf_caller = config.get("caller", "BOTH")
glnexus_filtration = config.get("glnexus_filtration", "custom")
sample_types = config.get("sample_types","WES")
genotype_mode = config.get("genotype_mode", "WES") #or WGS

print(f"Caller: {gvcf_caller}")
print(f"Sample type: {sample_types}")
print(f"Filtration setting: {glnexus_filtration}")
parts = level2_regions_diploid

def generate_gvcf_input_HC(wildcards):
    region = wildcards['region']
    res = []
    for samplefile in SAMPLE_FILES:
        sample_names = SAMPLEFILE_TO_SAMPLES[samplefile]
        samplefile_folder = get_samplefile_folder(samplefile)
        sample_sex = read_sexchrom(pj(samplefile_folder, samplefile + '.sex_chrom.tab'))
        gvcf_input = []
        for sample in sample_names:
            if sample_sex[sample] == 'F' and (region.startswith('Y') or region.endswith('H')):
                continue
            # Determine if it is WGS or WES
            if 'wgs' in SAMPLEINFO[sample]["sample_type"] or "WGS" in SAMPLEINFO[sample]["sample_type"] :
                chunk = convert_to_level1(region)
                if genotype_mode == 'WES':
                    filenames = expand("{cd}/{GVCF}/exome_extract/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=GVCF,region = chunk, sample=sample,allow_missing=True)
                else:
                    filenames = expand("{cd}/{GVCF}/reblock/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=GVCF,region = chunk, sample=sample,allow_missing=True)
            else:  # WES
                chunk = convert_to_level0(region)
                filenames = expand("{cd}/{GVCF}/reblock/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=GVCF,region = chunk, sample=sample,allow_missing=True)
            gvcf_input.extend(filenames)
        res.extend(gvcf_input)
    return res

def generate_gvcf_input_DV(wildcards):
    region = wildcards['region']
    res = []
    for samplefile in SAMPLE_FILES:
        sample_names = SAMPLEFILE_TO_SAMPLES[samplefile]
        samplefile_folder = get_samplefile_folder(samplefile)
        sample_sex = read_sexchrom(pj(samplefile_folder, samplefile + '.sex_chrom.tab'))
        gvcf_input = []
        for sample in sample_names:
            if sample_sex[sample] == 'F' and (region.startswith('Y') or region.endswith('H')):
                continue
            # Determine if it is WGS or WES
            if 'wgs' in SAMPLEINFO[sample]["sample_type"] or "WGS" in SAMPLEINFO[sample]["sample_type"] :
                chunk = convert_to_level1(region)
                if genotype_mode == 'WES':
                    filenames = expand("{cd}/{GVCF}/gVCF/exome_extract/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=DEEPVARIANT,region = chunk, sample=sample,allow_missing=True)
                else:
                    filenames = expand("{cd}/{GVCF}/gVCF/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=DEEPVARIANT,region = chunk, sample=sample,allow_missing=True)
            else:  # WES
                chunk = convert_to_level0(region)
                filenames = expand("{cd}/{GVCF}/gVCF/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=DEEPVARIANT,region = chunk, sample=sample,allow_missing=True)
            gvcf_input.extend(filenames)
        res.extend(gvcf_input)
    return res

if gvcf_caller == "HaplotypeCaller":
    glnexus_dir = ["GLnexus_on_Haplotypecaller"]

elif gvcf_caller == "Deepvariant":
    glnexus_dir = ["GLnexus_on_Deepvariant"]

elif gvcf_caller == "BOTH":
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
    region = wildcards['region']
    return region_to_file(region,wgs=genotype_mode == 'WGS',extension='bed')



rule GLnexus_all:
    input:
        expand("{cur_dir}/{genotype_mode}_{types_of_gl}{appendix}/{region}.vcf.gz", genotype_mode = genotype_mode, cur_dir = current_dir, region  = parts, types_of_gl = glnexus_dir, appendix = dir_appendix),
        expand("{cur_dir}/{genotype_mode}_{types_of_gl}{appendix}/{region}.vcf.gz.tbi", genotype_mode = genotype_mode, cur_dir = current_dir, region = parts, types_of_gl = glnexus_dir, appendix = dir_appendix),
        expand("{cur_dir}/{genotype_mode}_{types_of_gl}{appendix}/ANNOTATED/{region}_annotated.hg38_multianno.vcf.gz", genotype_mode= genotype_mode, cur_dir = current_dir, region  = parts, types_of_gl = glnexus_dir, appendix = dir_appendix),
    default_target: True

rule glnexus_HC:
    input: generate_gvcf_input_HC
    output: vcf = pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Haplotypecaller" + dir_appendix, "{region}.vcf.gz")
    container: "docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    params: bed = region_to_bed_file,
            mem_gb = 9,
            scratch_dir =  temp(current_dir + '/' + tmpdir + "/{genotype_mode}_{region}_glnexus.DB"),
            conf_filters = conf_filter
    threads: 6
    log: pj(current_dir, "logs", "glnexus", "glnexus_HC_{region}.{genotype_mode}.log")
    resources: 
        n = "6",
        mem_mb = 14000
    shell:
        """
        rm -rf {params.scratch_dir} &&
        glnexus_cli  --dir {params.scratch_dir} --bed {params.bed} --threads 5 --mem-gbytes {params.mem_gb} --config {params.conf_filters}  {input} 2> {log}  |  bcftools view -  | bgzip -@ 6 -c > {output} 2>> {log}
        """
rule index_deep:
    input: rules.glnexus_HC.output.vcf
    output: tbi = pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Haplotypecaller" + dir_appendix, "{region}.vcf.gz.tbi")
    conda: CONDA_VCF
    resources: n = "1"
    shell: "gatk IndexFeatureFile -I {input}"



use rule glnexus_HC as glnexus_DV with:
    input: gvcf_input = generate_gvcf_input_DV
    output: vcf= pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Deepvariant" + dir_appendix, "{region}.vcf.gz")
    log: pj(current_dir,"logs","glnexus","glnexus_DV_{region}.{genotype_mode}.log")
    params: scratch_dir =  temp(current_dir + '/' + tmpdir + "/{genotype_mode}_{region}_glnexus_2.DB"),
            bed= region_to_bed_file,
            mem_gb= 9,
            conf_filters= conf_filter
use rule index_deep as index_deep_2 with:
    input: rules.glnexus_DV.output.vcf
    output: tbi = pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Deepvariant" + dir_appendix, "{region}.vcf.gz.tbi")

rule annotate_genes:
    input: vcf = pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix +  "/{region}.vcf.gz"),
            tbi = pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix +  "/{region}.vcf.gz.tbi")
    output: vcf_annotated = pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED" , "{region}_annotated.hg38_multianno.vcf.gz"),
    conda: CONDA_ANNOVAR
    resources: n = "2"
    params: temp_vcf = pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp" , "{region}_annotated.vcf.gz"),
            out = pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED" , "{region}_annotated")
    shell:
        """
        bcftools annotate -a {REVEL} -h {REVEL_header} -c CHROM,POS,REF,ALT,REVEL {input.vcf} -O v -o {params.temp_vcf} --threads {resources.n} &&
        
        perl {annovar} {params.temp_vcf} {annovar_db} -out {params.out} -protocol ensGene,refGene -operation g,g -vcfinput -buildver hg38 -thread {resources.n}
        
        rm -rf {params.temp_vcf} 
        """