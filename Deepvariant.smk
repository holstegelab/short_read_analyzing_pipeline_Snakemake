wildcard_constraints:
    sample="[\w\d_\-@]+",
    region = "[\w\d]+",
    # readgroup="[\w\d_\-@]+"
onsuccess: shell("rm -fr logs/*")
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
    return [pj(DEEPVARIANT,  'gVCF', region, f'{sample}.{region}.wg.vcf.gz') for region in regions if not region.endswith('H')]#}}}
    # return [pj(DEEPVARIANT,'gVCF',region,f'{sample}.{region}.wg.vcf.gz') for region in regions]  #}}}

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
    res = 3500 if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else 3000
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
        bed = region_to_bed_file
    output:
        vcf = temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz")),
        vcf_tbi = temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz.tbi")),
        gvcf = temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz")),
        gvcf_tbi = temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz.tbi")),
        inter_dir = temp(directory(pj(DEEPVARIANT,'DV_intermediate', "{sample}.{region}")))
    params: 
            cd = current_dir + '/',
            mode=get_sequencing_mode
    container: 'docker://google/deepvariant:1.5.0'
    benchmark:
        pj(BENCH,"{sample}.{region}.wholedeepvariant.txt")
    resources:
        n="1.5", #set in profile using singularity-args. Waiting for rule-specific args. 
        nshards=2,
        mem_mb=get_mem_mb_deepvariant
    log: pj(LOG,"{sample}.{region}.wholedeepvariant.log")
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant --make_examples_extra_args="normalize_reads=true" --call_variants_extra_args config_string="device_count:{{key:'CPU' value:4}} inter_op_parallelism_threads:4 intra_op_parallelism_threads:4" --num_shards={resources.nshards} --model_type={params.mode} --regions={input.bed} --ref={REF_MALE} --reads={params.cd}{input.bam} --output_vcf={output.vcf} --output_gvcf={output.gvcf} --intermediate_results_dir "{output.inter_dir}"  2> {log}
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
        merge_script=srcdir(MERGEPHASE)
    resources: 
        n="1.0",
        mem_mb = 600
    conda: CONDA_VCF
    shell: """
        mkdir -p `dirname {output.wstats}`
        whatshap phase  --ignore-read-groups --reference {REF} {input.vcf} {input.bams} -o {output.vcf}        
        whatshap stats {output.vcf} > {output.wstats}
        python {params.merge_script} {input.gvcf} {output.vcf} {output.tmp_gvcf} {output.mwstats}
        bcftools view {output.tmp_gvcf} -o {output.gvcf}
        bcftools index --tbi {output.gvcf}
        """        


# for i in {0..1}; do singularity run -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:"1.4.0" /opt/deepvariant/bin/make_examples --ref=/gpfs/work3/0/qtholstg/hg38_res/Ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --reads=/gpfs/work3/0/qtholstg/Georgii_tests/10_additional_samples_gVCF/bams/NL_VUMC_KG-013832.markdup.bam --gvcf=/gpfs/work3/0/qtholstg/Georgii_tests/deepvariant_test/NL_VUMC_KG-01382_deepvariant_call_variant_test.gvcf.tfrecord_chr20@2.gz --mode calling --regions chr20 --examples /gpfs/work3/0/qtholstg/Georgii_tests/deepvariant_test/NL_VUMC_KG-01382_deepvariant_call_variant_test_examples_chr20@2.gz --task $i --channels="insert_size"    --gvcf_gq_binsize 3; done
#
# singularity run -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:"1.4.0" /opt/deepvariant/bin/call_variants --examples /gpfs/work3/0/qtholstg/Georgii_tests/deepvariant_test/NL_VUMC_KG-01382_deepvariant_call_variant_test_examples_chr20@2.gz --outfile /gpfs/work3/0/qtholstg/Georgii_tests/deepvariant_test/NL_VUMC_KG-01382_deepvariant_step2_out --checkpoint /opt/models/wgs/model.ckpt
#
# singularity run -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:"1.4.0" /opt/deepvariant/bin/postprocess_variants --ref=/gpfs/work3/0/qtholstg/hg38_res/Ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --infile /gpfs/work3/0/qtholstg/Georgii_tests/deepvariant_test/NL_VUMC_KG-01382_deepvariant_step2_out --outfile /gpfs/work3/0/qtholstg/Georgii_tests/deepvariant_test/NL_VUMC_KG-01382_deepvariant_step3.vcf --gvcf_outfile /gpfs/work3/0/qtholstg/Georgii_tests/deepvariant_test/NL_VUMC_KG-01382_deepvariant_step3.g.vcf --nonvariant_site_tfrecord_path /gpfs/work3/0/qtholstg/Georgii_tests/deepvariant_test/NL_VUMC_KG-01382_deepvariant_call_variant_test.gvcf.tfrecord_chr20@2.gz
#

#def get_model(wildcards):
#     mode =  "WGS" if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else "WES"
#     if mode == 'WES':
#         model = pj("/opt/models/wes/model.ckpt")
#     elif mode == 'WGS':
#         model = pj("/opt/models/wgs/model.ckpt")
#     return model


#rule deepvariant2:
#    input:
#        bam = rules.markdup.output.mdbams,
#        bai = rules.markdup.output.mdbams_bai,
#        bed = region_to_bed_file
#    output:
#        vcf = temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz")),
#        vcf_tbi = temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz.tbi")),
#        gvcf = temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz")),
#        gvcf_tbi = temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz.tbi")),
#        pregvcf = temp(pj(DEEPVARIANT,'intermediate', "{sample}.{region}.gvcf.trfrecord.gz")),
#        examples = temp(pj(DEEPVARIANT,'intermediate', "{sample}.{region}.examples.gz")),
#        precall = temp(pj(DEEPVARIANT,'intermediate', "{sample}.{region}.intermediate_out.gz")),
#    params: 
#            model = get_model,
#            mode=get_sequencing_mode
#    container: 'docker://google/deepvariant:1.5.0'
#    benchmark:
#        pj(BENCH,"{sample}.{region}.wholedeepvariant.txt")
#    resources:
#        n="1.5", #set in profile using singularity-args. Waiting for rule-specific args. 
#        nshards=2,
#        mem_mb=get_mem_mb_deepvariant
#    log: pj(LOG,"{sample}.{region}.wholedeepvariant.log")
#    shell:
#        """
#        mkdir -p `dirname {output.pregvcf}`
#
#        /opt/deepvariant/bin/make_examples --ref={REF_MALE} --reads={input.bam} --gvcf={output.pregvcf} --mode calling --regions {input.bed} --examples={output.examples} --channels="insert_size" --normalize_reads tru --gvcf_gq_binsize 1
#        /opt/deepvariant/bin/call_variants --examples {output.examples} --outfile {output.precall} --checkpoint {params.model} --config_string "device_count:{{key:'CPU' value:4}} inter_op_parallelism_threads:4 intra_op_parallelism_threads:4"
#        /opt/deepvariant/bin/postprocess_variants --ref={REF_MALE}  --infile {output.precall} --outfile {output.vcf} --gvcf_outfile {output.gvcf} --nonvariant_site_tfrecord_path {output.pregvcf} 
#        """



#         
# rule make_examples:
#     input:
#         bam =  pj(current_dir,rules.markdup.output.mdbams),
#         bai = pj(current_dir,rules.markdup_index.output.mdbams_bai)
#     output:
#         pregvcf = pj(current_dir, config['DEEPVARIANT'],'intermediate', "{region}.{sample}.{mode}.gvcf.trfrecord.gz"),
#         examples = pj(current_dir, config['DEEPVARIANT'],'intermediate', "{region}.{sample}.{mode}.examples.gz"),
#     threads: config['make_examples']['n']
#     params: inter_dir = pj(current_dir, config['DEEPVARIANT'],'intermediate', "{region}.{sample}.{mode}"),
#             # change to merged bed capture kit divided per BINs
#             bed= get_chrom_capture_kit_bed
#     container: 'docker://google/deepvariant:1.4.0'
#     benchmark:
#         pj(current_dir, config['BENCH'],"{region}.{sample}.{mode}.makeexamples.txt")
#     log: pj(current_dir, config['LOG'],"{region}.{sample}.{mode}.makeexamples.log")
#     shell:
#         # singularity run -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:"1.4.0" \
#         """
#          """
#
#
# rule call_variants:
#     input:
#         examples =  pj(current_dir,rules.make_examples.output.examples),
#     output:
#         precall = pj(current_dir, config['DEEPVARIANT'],'intermediate', "{region}.{sample}.{mode}.intermediate_out.gz"),
#     params: model = get_model
#     container: 'docker://google/deepvariant:1.4.0'
#     threads: config['call_variants']['n']
#     benchmark:
#         pj(current_dir, config['BENCH'],"{region}_{sample}.{mode}.call_variants.txt")
#     log: pj(current_dir, config['LOG'],"{region}_{sample}.{mode}.call_variants.log")
#     shell:
#         """
#         /opt/deepvariant/bin/call_variants --examples {input.examples} --outfile {output.precall} --checkpoint {params.model} 2> {log}
#         """
#
# rule postprocess_variants:
#     input:
#         non_variant_records =  pj(current_dir,rules.make_examples.output.pregvcf),
#         precall= pj(current_dir,rules.call_variants.output.precall),
#     output:
#         vcf = pj(current_dir, config['DEEPVARIANT'],'VCF', "{region}.{sample}.{mode}.vcf.gz"),
#         gvcf = pj(current_dir, config['DEEPVARIANT'],'gVCF', "{region}.{sample}.{mode}.g.vcf.gz"),
#     params: model = get_model
#     container: 'docker://google/deepvariant:1.4.0'
#     benchmark:
#         pj(current_dir, config['BENCH'],"{region}_{sample}.{mode}.postprocess_variants.txt")
#     log: pj(current_dir, config['LOG'],"{region}_{sample}.{mode}.postprocess_variants.log")
#     threads: config['postprocess_variants']['n']
#     shell:
#             """
#             /opt/deepvariant/bin/postprocess_variants --ref={ref}  --infile {input.precall} --outfile {output.vcf} --gvcf_outfile {output.gvcf} --nonvariant_site_tfrecord_path {input.non_variant_records} 2> {log}
#             """
