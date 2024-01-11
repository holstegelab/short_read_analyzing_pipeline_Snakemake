wildcard_constraints:
    sample="[\w\d_\-@]+",
    region = "[\w\d]+",
    # readgroup="[\w\d_\-@]+"
onsuccess: shell("rm -fr logs/Deepvariant/*")
from common import *


module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner
module Tools:
    snakefile: 'Tools.smk'
    config: config
use rule * from Tools
mode = config.get("computing_mode", "WES")

rule DeepVariant_all:
    input:
        expand("{dv}/{sample}.done",sample=sample_names, dv = DEEPVARIANT),
        rules.Aligner_all.input
    default_target: True


def get_deepvariant_files(wildcards):#{{{
    sample = wildcards['sample']
    regions = level1_regions if 'wgs' in SAMPLEINFO[sample]['sample_type'] else level0_regions
    return [pj(DEEPVARIANT,  'gVCF', region, f'{sample}.{region}.wg.vcf.gz') for region in regions]#}}}

rule deepvariant_sample_done:
    input:
        get_deepvariant_files
    output:
        temp(touch(pj(DEEPVARIANT, "{sample}.done")))
    resources:
        mem_mb = 100,
        n = "1.0"



def get_sequencing_mode(wildcards):#{{{
    return "WGS" if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else "WES"#}}}

def get_mem_mb_deepvariant(wildcards, attempt):#{{{
    res = 10000 if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else 9000
    return (attempt - 1) * 0.5 * res + res#}}}


def region_to_bed_file(wildcards):#{{{
    """Converts a region to a bed file location (see common.py and Tools.smk)"""
    sample = wildcards['sample']
    region = wildcards['region']
    return region_to_file(region, wgs='wgs' in SAMPLEINFO[sample]['sample_type'], extension='bed')#}}}

rule deepvariant:
    input:
        bam = rules.markdup.output.mdbams,
        bai = rules.markdup.output.mdbams_bai,
        bed = region_to_bed_file,
        validated_sex = rules.get_validated_sex.output.yaml
    output:
        vcf = temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz")),
        vcf_tbi = temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz.tbi")),
        gvcf = temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz")),
        gvcf_tbi = temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz.tbi")),
    params:
            cd = current_dir + '/',
            mode=get_sequencing_mode,
            haploid_contigs=lambda wildcards: 'chrX,chrX_KI270880v1_alt,chrX_KI270881v1_alt,chrX_KI270913v1_alt,chrY,chrY_KI270740v1_random' if wildcards['region'].endswith("H") else 'chrNONE',
            skipsex = lambda wildcards, input: int(get_validated_sex_file(input) == 'female' and wildcards['region'].startswith('Y')),
            inter_dir = pj(DEEPVARIANT,'DV_intermediate')
    container: 'docker://google/deepvariant:1.6.0'
    resources:
        n="2.5",
        nshards=4,
        mem_mb=get_mem_mb_deepvariant
    shell:
        """
        if [ {params.skipsex} -eq 0 ]
        then
            mkdir -p "{params.inter_dir}/{wildcards.sample}.{wildcards.region}"
            /opt/deepvariant/bin/run_deepvariant --make_examples_extra_args="normalize_reads=true" --call_variants_extra_args config_string="device_count:{{key:'CPU' value:4}} inter_op_parallelism_threads:4 intra_op_parallelism_threads:4" --num_shards={resources.nshards} --model_type={params.mode} --regions={input.bed} --ref={REF_MALE} --reads={params.cd}{input.bam} --output_vcf={output.vcf} --output_gvcf={output.gvcf} --haploid_contigs {params.haploid_contigs} --intermediate_results_dir "{params.inter_dir}/{wildcards.sample}.{wildcards.region}" --postprocess_cpus 4
            rm -rf "{params.inter_dir}/{wildcards.sample}.{wildcards.region}"
        else
            touch {output.vcf}
            touch {output.vcf_tbi}
            touch {output.gvcf}
            touch {output.gvcf_tbi}
        fi
        """


rule DVWhatshapPhasingMerge:
    """Phase VCF with Whatshap and merge into the gVCF"""
    input:
        vcf = rules.deepvariant.output.vcf,
        vcf_tbi = rules.deepvariant.output.vcf_tbi,
        gvcf = rules.deepvariant.output.gvcf,
        gvcf_tbi = rules.deepvariant.output.gvcf_tbi,
        bams = rules.markdup.output.mdbams,
        bai = rules.markdup.output.mdbams_bai,
        validated_sex = rules.get_validated_sex.output.yaml
    output:
        vcf = pj(DEEPVARIANT, "VCF/{region}/{sample}.{region}.w.vcf.gz"),
        wstats = pj(STAT, "whatshap_dvphasing/{sample}.{region}.stats"),
        mwstats = pj(STAT, "whatshap_dvphasing/{sample}.{region}.merge_stats"),
        tmp_gvcf= temp(pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf")),
        gvcf= pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz"),
        gvcf_tbi = pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz.tbi"),
    params:
        merge_script=srcdir(MERGEPHASE),
        ploidy=lambda wildcards: 1 if wildcards["region"].endswith("H") else 2,
        skipsex = lambda wildcards, input: int(get_validated_sex_file(input) == 'female' and wildcards['region'].startswith('Y'))
    resources: 
        n="1.0",
        mem_mb = 900
    conda: CONDA_VCF
    shell: """
        mkdir -p `dirname {output.wstats}`
        if [ {params.ploidy} -eq 2 ]
        then 
            whatshap phase  --ignore-read-groups --reference {REF} {input.vcf} {input.bams} -o {output.vcf}        
            whatshap stats {output.vcf} > {output.wstats}
            python {params.merge_script} {input.gvcf} {output.vcf} {output.tmp_gvcf} {output.mwstats}
            bcftools view {output.tmp_gvcf} -o {output.gvcf}
            bcftools index --tbi {output.gvcf}
        else
            touch {output.tmp_gvcf}
            touch {output.vcf}
            touch {output.vcf}.tbi
            touch {output.wstats}
            touch {output.mwstats}
            touch {output.tmp_gvcf}
            if [ {params.skipsex} -eq 0 ]
            then
                bcftools view {input.gvcf} -o {output.gvcf}
                bcftools index --tbi {output.gvcf}
            else
                touch {output.gvcf}
                touch {output.gvcf_tbi}
            fi
        fi
        """


rule extract_exomes:
    input:
        gvcf = pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz"),
        tbi = pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz.tbi"),
    output:
        gvcf_exome = ensure(pj(DEEPVARIANT, "gVCF/exomes/{region}/{sample}.{region}.wg.vcf.gz"), non_empty=True),
        tbi = ensure(pj(DEEPVARIANT, "gVCF/exomes/{region}/{sample}.{region}.wg.vcf.gz.tbi"), non_empty=True),
    conda: CONDA_VCF
    params: java_options=DEFAULT_JAVA_OPTIONS,
            interval = lambda wildcards: region_to_file(region = wildcards.region, extension="interval_list"),
            padding = 500,
    resources: n= "1.0",
               mem_mb= 1500,
    run:
        if SAMPLEINFO[wildcards.sample]["sample_type"] == "wgs":
            shell(
                """
                    gatk --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" SelectVariants \
                    -V {input.gvcf} -O {output.gvcf_exome} \
                    -L {params.interval} -ip {params.padding} --seconds-between-progress-updates 120 \
                    -G StandardAnnotation -G AS_StandardAnnotation
                """),
        else:
            shell(
                """
                    cp {input.gvcf} {output.gvcf_exome}
                    cp {input.tbi} {output.tbi}
                """)
