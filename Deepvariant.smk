wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    region = r"[\w\d]+",
    # readgroup="[\w\d_\-@]+"
onsuccess: shell("rm -fr logs/Deepvariant/*")
from common import *


module Tools:
    snakefile: 'Tools.smk'
    config: config
use rule * from Tools
mode = config.get("computing_mode", "WES")

rule DeepVariant_all:
    input:
        expand("{dv}/{sample}.done",sample=sample_names, dv = DEEPVARIANT),
    default_target: True


def get_deepvariant_files(wildcards):#{{{
    sample = wildcards['sample']
    if 'wgs' in SAMPLEINFO[sample]['sample_type']:
        return [pj(DEEPVARIANT,  'gVCF', 'exome_extract', region, f'{sample}.{region}.wg.vcf.gz') for region in level1_regions]
    else:
        return [pj(DEEPVARIANT,  'gVCF', region, f'{sample}.{region}.wg.vcf.gz') for region in level0_regions]
#}}}


rule deepvariant_sample_done:
    input:
        get_deepvariant_files
    output:
        touch(pj(DEEPVARIANT, "{sample}.done"))
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
        bed = region_to_bed_file,
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        validated_sex=pj(KMER,"{sample}.result.yaml"),
    output:
        vcf = ensure(temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz")), non_empty=True),
        vcf_tbi = ensure(temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz.tbi")),non_empty=True),
        gvcf = ensure(temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz")),non_empty=True),
        gvcf_tbi = ensure(temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz.tbi")),non_empty=True),
    log: pj(LOG, "Deepvariant", "{sample}.{region}.log"),
    params:
            cd = current_dir + '/',
            mode=get_sequencing_mode,
            haploid_contigs=lambda wildcards: 'chrX,chrX_KI270880v1_alt,chrX_KI270881v1_alt,chrX_KI270913v1_alt,chrY,chrY_KI270740v1_random' if wildcards['region'].endswith("H") else 'chrNONE',
            skipsex = lambda wildcards, input: int(get_validated_sex_file(input) == 'female' and wildcards['region'].startswith('Y')),
            inter_dir = pj(current_dir, DEEPVARIANT,'DV_intermediate'),
            # check = CHECKEMPTY
    container: 'docker://google/deepvariant:1.6.1'
    resources:
        n="2.5",
        nshards=4,
        mem_mb=get_mem_mb_deepvariant,
        time = 6600
    shell:
        #         mkdir -p "{params.inter_dir}/{wildcards.sample}.{wildcards.region}"
        """
        mkdir -p "{params.cd}{params.inter_dir}/{wildcards.sample}.{wildcards.region}"
        if [ {params.skipsex} -eq 0 ]
        then

            (/opt/deepvariant/bin/run_deepvariant --make_examples_extra_args="normalize_reads=true" --call_variants_extra_args config_string="device_count:{{key:'CPU' value:4}} inter_op_parallelism_threads:4 intra_op_parallelism_threads:4" --num_shards={resources.nshards} --model_type={params.mode} --regions={input.bed} --ref={REF_MALE} --reads={params.cd}{input.bam} --output_vcf={params.cd}{output.vcf} --output_gvcf={params.cd}{output.gvcf} --haploid_contigs {params.haploid_contigs} --intermediate_results_dir "{params.cd}{params.inter_dir}/{wildcards.sample}.{wildcards.region}" --postprocess_cpus 4 )
            rm -rf "{params.cd}{params.inter_dir}/{wildcards.sample}.{wildcards.region}"

        else
            touch {output.vcf}
            touch {output.vcf_tbi}
            touch {output.gvcf}
            touch {output.gvcf_tbi}
        fi
        """

# python {params.check} {output.vcf}
# python {params.check} {output.gvcf}

rule DVWhatshapPhasingMerge:
    """Phase VCF with Whatshap and merge into the gVCF"""
    input:
        vcf = rules.deepvariant.output.vcf,
        vcf_tbi = rules.deepvariant.output.vcf_tbi,
        gvcf = rules.deepvariant.output.gvcf,
        gvcf_tbi = rules.deepvariant.output.gvcf_tbi,
        bams=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        validated_sex=pj(KMER,"{sample}.result.yaml"),
    output:
        vcf = pj(DEEPVARIANT, "VCF/{region}/{sample}.{region}.w.vcf.gz"),
        wstats = pj(STAT, "whatshap_dvphasing/{sample}.{region}.stats"),
        mwstats = pj(STAT, "whatshap_dvphasing/{sample}.{region}.merge_stats"),
        tmp_gvcf= temp(pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf")),
        gvcf= pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz"),
        gvcf_tbi = pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz.tbi"),
    params:
        merge_script=srcdir(MERGEPHASEDIRECT),
        ploidy=lambda wildcards: 1 if wildcards["region"].endswith("H") else 2,
        skipsex = lambda wildcards, input: int(get_validated_sex_file(input) == 'female' and wildcards['region'].startswith('Y'))
    log: pj(LOG, "Deepvariant", "{sample}.{region}.whatshap.log"),
    resources: 
        n="1.0",
        mem_mb = 900
    conda: CONDA_VCF
    shell: """
        mkdir -p `dirname {output.wstats}`
        if [ {params.ploidy} -eq 2 ] && [ {params.skipsex} -eq 0 ]
        then 
            if [ {params.skipsex} -eq 0 ]
            then 
                whatshap phase  --ignore-read-groups --reference {REF} {input.vcf} {input.bams} -o {output.vcf}        
                whatshap stats {output.vcf} > {output.wstats}
                python {params.merge_script} {input.gvcf} {output.vcf} {output.tmp_gvcf} {output.mwstats}
                bcftools view {output.tmp_gvcf} -o {output.gvcf}
                bcftools index --tbi {output.gvcf}
            else
                touch {output.vcf}
                touch {output.vcf}.tbi
                touch {output.wstats}
                touch {output.mwstats}
                touch {output.tmp_gvcf}
                touch {output.gvcf}
                touch {output.gvcf_tbi}
            fi
        else
            cp {input.vcf} {output.vcf}
            cp {input.vcf_tbi} {output.vcf}.tbi
            cp {input.gvcf} {output.gvcf}
            cp {input.gvcf_tbi} {output.gvcf_tbi}

            touch {output.tmp_gvcf}
            touch {output.wstats}
            touch {output.mwstats}
            touch {output.tmp_gvcf}
        fi
        """


rule extract_exomes_dv:
    input:
        gvcf = rules.DVWhatshapPhasingMerge.output.gvcf,
        tbi = rules.DVWhatshapPhasingMerge.output.gvcf_tbi,
        validated_sex=pj(KMER,"{sample}.result.yaml"),
    output:
        gvcf_exome = pj(DEEPVARIANT, "gVCF/exome_extract/{region}/{sample}.{region}.wg.vcf.gz"),
        tbi = pj(DEEPVARIANT, "gVCF/exome_extract/{region}/{sample}.{region}.wg.vcf.gz.tbi"),
    conda: CONDA_VCF
    params: java_options=DEFAULT_JAVA_OPTIONS,
            interval = lambda wildcards: region_to_file(region = wildcards.region, extension="interval_list"),
            padding = 500,
            skipsex= lambda wildcards,input: int(get_validated_sex_file(input) == 'female' and wildcards['region'].startswith('Y')),
            check = CHECKEMPTY
    resources: n= "1.0",
               mem_mb= 1500,

    shell:
        """
        if [ {params.skipsex} -eq 0 ]
        then
            gatk --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" SelectVariants \
            -V {input.gvcf} -O {output.gvcf_exome} \
            -L {params.interval} -ip {params.padding} --seconds-between-progress-updates 120 
            python {params.check} {output.gvcf_exome} 
            python {params.check} {output.tbi} 
        else
            touch {output.gvcf_exome}
            touch {output.tbi}
        fi
        """
