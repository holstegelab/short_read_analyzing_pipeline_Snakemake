import os
from common import *
current_dir = os.getcwd()
import concurrent.futures
# onsuccess: "rm -rf "

wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    chr = r"[\w\d]+",
    region = r"[\w\d]+",
    genotype_mode = r"WES|WGS",
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

def generate_gvcf_input_HC_divided(wildcards):
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
            filenames = expand("{cd}/{region}_gvcfs_HC/{i}.vcf", cd = current_dir, i = sample, region = region , allow_missing=True)
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
                # 'gVCF', 'DIVIDED', '{region}', '{sample}.{region}.wg.vcf.gz'
                filenames = expand("{cd}/{GVCF}/gVCF/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=DEEPVARIANT,region = chunk, sample=sample,allow_missing=True)
            gvcf_input.extend(filenames)
        res.extend(gvcf_input)
    return res

def generate_gvcf_input_DV_divided(wildcards):
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
            filenames = expand("{cd}/{region}_gvcfs_DV/{i}.vcf", cd = current_dir, i = sample, region = region ,  allow_missing=True)
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
        expand("{cur_dir}/{genotype_mode}_{types_of_gl}{appendix}/ANNOTATED/{region}.annotated.vcf.gz", genotype_mode= genotype_mode, cur_dir = current_dir, region  = parts, types_of_gl = glnexus_dir, appendix = dir_appendix),
    default_target: True




def run_bcftools_HC(i, bed, region):
    filename = i.split("/")[-1]
    samplename = filename.split(".")[0]
    cmd = f"bcftools view -R {bed} {i} -O v -o {region}_gvcfs_HC/{samplename}.vcf"
    return cmd

rule glnexus_HC:
    input: generate_gvcf_input_HC
    output: vcf = pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Haplotypecaller" + dir_appendix, "{region}.vcf.gz"),
    log: pj(current_dir,"logs","glnexus","glnexus_HC_{region}.{genotype_mode}.log")
    container: "docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    params: bed = region_to_bed_file,
            mem_gb = 160,
            scratch_dir =  temp(current_dir + '/' + tmpdir + "/{genotype_mode}_{region}_glnexus_HC.DB"),
            conf_filters = conf_filter,
            generate_gvcf_input_HC_divided = generate_gvcf_input_HC_divided
    threads: 64
    resources:
        n = "64",
        mem_mb = 160000,
        active_use_add= 100
    run:
        shell("mkdir -p {wildcards.region}_gvcfs_HC")
        cmds = [run_bcftools_HC(i, params.bed, wildcards.region) for i in input]
        with concurrent.futures.ProcessPoolExecutor(max_workers=64) as executor:
            executor.map(shell, cmds)

        shell("""
        rm -rf {params.scratch_dir} &&
        glnexus_cli  --dir {params.scratch_dir} --bed {params.bed} --threads 63 --mem-gbytes {params.mem_gb} --config {params.conf_filters}  {params.generate_gvcf_input_HC_divided} 2> {log}  |  bcftools view --threads 64 -  | bgzip -@ 64 -c > {output} 2>> {log}
        rm -rf {wildcards.region}_gvcfs_HC
        """)
rule index_deep:
    input: rules.glnexus_HC.output.vcf
    output: tbi = pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Haplotypecaller" + dir_appendix, "{region}.vcf.gz.tbi")
    conda: CONDA_VCF
    resources: n = "1",
            active_use_remove = 500,
            mem_mb = 2500
    shell: "gatk IndexFeatureFile -I {input}"


def run_bcftools_DV(i, bed, region):
    filename = i.split("/")[-1]
    samplename = filename.split(".")[0]
    cmd = f"bcftools view -R {bed} {i} -O v -o {region}_gvcfs_DV/{samplename}.vcf"
    return cmd

rule glnexus_DV:
    input: generate_gvcf_input_DV
    output: vcf= pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Deepvariant" + dir_appendix, "{region}.vcf.gz"),
    log: pj(current_dir,"logs","glnexus","glnexus_DV_{region}.{genotype_mode}.log")
    params: bed = region_to_bed_file,
            mem_gb = 120,
            scratch_dir =  temp(current_dir + '/' + tmpdir + "/{genotype_mode}_{region}_glnexus_DV.DB"),
            conf_filters = conf_filter,
            generate_gvcf_input_DV_divided = generate_gvcf_input_DV_divided
    container: "docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    threads: 64
    resources:
        n = "64",
        mem_mb = 120000,
        active_use_add= 100
    run:
        shell("mkdir -p {wildcards.region}_gvcfs_DV")
        cmds = [run_bcftools_DV(i, params.bed, wildcards.region) for i in input]
        with concurrent.futures.ProcessPoolExecutor(max_workers=64) as executor:
            executor.map(shell, cmds)
        shell("""
        rm -rf {params.scratch_dir} &&
        glnexus_cli  --dir {params.scratch_dir} --bed {params.bed} --threads 63 --mem-gbytes {params.mem_gb} --config {params.conf_filters}  {params.generate_gvcf_input_DV_divided} 2> {log}  |  bcftools view --threads 64 -  | bgzip -@ 64 -c > {output} 2>> {log}
        rm -rf {wildcards.region}_gvcfs_DV
        """)
use rule index_deep as index_deep_2 with:
    input: rules.glnexus_DV.output.vcf
    output: tbi = pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Deepvariant" + dir_appendix, "{region}.vcf.gz.tbi")

rule check_glnexus_lof_file:
    input: vcf = pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Haplotypecaller" + dir_appendix, "{region}.vcf.gz"),
            log= pj(current_dir,"logs","glnexus","glnexus_HC_{region}.{genotype_mode}.log")
    output: pj(current_dir, "{genotype_mode}_GLnexus_on_Haplotypecaller", "{region}.vcf_is_ok")# file to staret annotation {wildcard.region}
    shell: #check if 'Finish' message is in log file, then write output, otherwise delete input vcf file
        """
        if grep -q "genotyping complete!" {input.log}; then
            touch {output}
        else
            mkdir error_vcfs
            mv {input.vcf} error_vcfs/
        fi
         """

use rule check_glnexus_lof_file as check_glnexus_lof_file_2 with:
    input: vcf = pj(current_dir, "{genotype_mode}_" + "GLnexus_on_Deepvariant" + dir_appendix, "{region}.vcf.gz"),
            log= pj(current_dir,"logs","glnexus","glnexus_DV_{region}.{genotype_mode}.log")
    output: pj(current_dir,"{genotype_mode}_GLnexus_on_Deepvariant","{region}.vcf_is_ok")  # file to staret annotation {wildcard.region}


rule split_multiallelic:
    input: vcf = pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix +  "/{region}.vcf.gz"),
            tbi= pj(current_dir,"{genotype_mode}_{types_of_gl}" + dir_appendix + "/{region}.vcf.gz.tbi"),
            logcheck=pj(current_dir, "{genotype_mode}_{types_of_gl}", "{region}.vcf_is_ok")
    output:
        temp_tsv = temp(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "{region}.split.tsv.gz")),
            vcf = temp(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "{region}.split.vcf.gz"))
    conda: CONDA_MAIN
    resources:
        mem_mb=24000,
        n=12
    shell:
        """
        {bcftools_patched} view -U -c1:nonmajor {input.vcf}  --threads {resources.n} -Ou | {bcftools_patched} query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' | awk '{{l2=split($4,a,","); if(l2 > 1 && $6 > 0) {{split($5,ac,","); ref=1.0; for (i in ac) ref-= (ac[i]/$6); failcount=0; for (i in ac) failcount+=(ref<(ac[i]/$6)); OFS="\t"; print $1, $2, $3, $4, $4, l2, $5, ref, failcount}}}}' | bgzip --threads {resources.n}  > {output.temp_tsv}
        tabix -s1 -b2 -e2 {output.temp_tsv}
        {bcftools_patched} view -U -c1:nonmajor {input.vcf} --threads {resources.n} -Ou | {bcftools_patched} annotate -a {output.temp_tsv} -h {multiallelic_hdr} -c CHROM,POS,REF,ALT,MA_ALT,NMA_ALT,MA_ALT_AC,MA_REF_AF,MA_REFLOW_COUNT --threads {resources.n} -Ov | {bcftools_patched} norm -f {REF} -m- -c w --thread {resources.n} -o {output.vcf} --output-type z
        """




rule extract_positions:
    input: vcf = pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix +  "/{region}.split.vcf.gz"),
            # tbi= pj(current_dir,"{genotype_mode}_{types_of_gl}" + dir_appendix + "/{region}.split.vcf.gz.tbi"),
            logcheck=pj(current_dir, "{genotype_mode}_{types_of_gl}", "{region}.vcf_is_ok")
    output: vcfgz = temp(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp" , "{region}_pos_only.vcf.gz"))
    conda: CONDA_MAIN
    params: vcf = temp(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp" , "{region}_pos_only.vcf"))
    shell:
        """
        tabix -fp vcf {input.vcf}
        
        bcftools view --drop-genotypes -O v -o {params.vcf} {input.vcf}
        bgzip {params.vcf}
        tabix -fp vcf {output.vcfgz}
        """

rule annotate_revel:
    input: vcf = pj(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp" , "{region}_pos_only.vcf.gz")),
    output: vcf_annotated = temp(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp" , "{region}.annotated_pos_only.vcf"))
    conda: CONDA_MAIN
    resources: n = "6",
            mem_mb = 8000
    log: pj(current_dir,"logs","glnexus","annotate_revel_{region}.{genotype_mode}.{types_of_gl}.log")
    params: temp_vcf = temp(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp" , "{region}_annotated_temp.vcf.gz")),
            temp_vcf_2 = temp(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp" , "{region}_annotated_2_temp.vcf.gz")),
            temp_vcf_3 = temp(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp" , "{region}_annotated_3_temp.vcf.gz")),
            temp_dir = pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp"),
    shell:
        """
        mkdir -p {params.temp_dir} &&

        bcftools annotate -a {GNOMAD_4} -c "INFO/Gnomad4_AF:=INFO/AF" -O z -o {params.temp_vcf} {input.vcf} --threads {resources.n} 2> {log}
        tabix -fp vcf {params.temp_vcf}
        bcftools annotate -a {CLINVAR} -c INFO -O z -o {params.temp_vcf_2} {params.temp_vcf} --threads {resources.n} 2>> {log}
        tabix -fp vcf {params.temp_vcf_2}
        bcftools annotate -a {GNOMAD_2} -c INFO/non_neuro_AF -O z -o {params.temp_vcf_3} {params.temp_vcf_2} --threads {resources.n} 2>> {log}
        tabix -fp vcf {params.temp_vcf_3}
        bcftools annotate -a {REVEL} -h {REVEL_header} -c CHROM,POS,REF,ALT,REVEL {params.temp_vcf_3} -O v -o {output.vcf_annotated} --threads {resources.n} 2>> {log}
        
        """
# INFO/non_neuro_AF
rule annotate_gene:
    input: temp_vcf = pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp" , "{region}.annotated_pos_only.vcf")
    output: vcf_annotated=temp(pj(current_dir,"{genotype_mode}_" + "{types_of_gl}" + dir_appendix,"ANNOTATED_temp","{region}.annotated.hg38_multianno.vcf")),
    conda: CONDA_ANNOVAR
    params:
        out=pj(current_dir,"{genotype_mode}_" + "{types_of_gl}" + dir_appendix,"ANNOTATED_temp","{region}.annotated"),
    log: pj(current_dir,"logs","glnexus","annotate_gene_{region}.{genotype_mode}.{types_of_gl}.log")
    resources: n = "2",
                mem_mb = 5000
    shell:
        """
        perl {annovar} {input.temp_vcf} {annovar_db} -out {params.out} -protocol ensGene,refGene,ljb26_all,dbnsfp42c,avsnp150 -operation g,g,f,f,f -vcfinput -buildver hg38 -thread {resources.n} 2> {log}
        """

rule bring_anno_to_samples:
    input: vcf_annotated = pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED_temp", "{region}.annotated.hg38_multianno.vcf"),
            samples_vcf = pj(current_dir, "{genotype_mode}_{types_of_gl}" + dir_appendix +  "/{region}.split.vcf.gz"),
    output: vcf_anno_samples = pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz"),
            # tbi = ensure(pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz.tbi"), non_empty=True)
    conda: CONDA_MAIN
    resources:
        n = "2",
        mem_mb = 6000
    shell:
        """
        bgzip --force {input.vcf_annotated}
        tabix -fp vcf {input.vcf_annotated}.gz
        bcftools annotate -a {input.vcf_annotated}.gz -c INFO -O z -o {output.vcf_anno_samples} {input.samples_vcf}  
        """


rule bgzip:
    input: vcf_annotated = pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz")
    output: # vcf_annotated_gz = pj(current_dir, "{genotype_mode}_" + "{types_of_gl}" + dir_appendix, "ANNOTATED", "{region}.annotated.vcf.gz"),
            vcf_annotated_gz_tbi= pj(current_dir,"{genotype_mode}_" + "{types_of_gl}" + dir_appendix,"ANNOTATED","{region}.annotated.vcf.gz.tbi")
    conda: CONDA_MAIN
    shell:
        """
        tabix -fp vcf {input.vcf_annotated}
        rm -rf {wildcards.genotype_mode}_{wildcards.types_of_gl}{dir_appendix}/ANNOTATED_temp/{wildcards.region}*
        """