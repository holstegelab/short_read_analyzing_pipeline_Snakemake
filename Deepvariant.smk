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
    chrom = "[\w\d]+",
    # readgroup="[\w\d_\-@]+"


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

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
        expand("{dv}/{sample}.done",sample=sample_names, dv = config['DEEPVARIANT']),
        rules.Aligner_all.input
    default_target: True


def get_deepvariant_files(wildcards):
    sample = wildcards['sample']
    res = []
    for chrom in main_chrs:
        res.append(os.path.join(config['DEEPVARIANT'], 'gVCF', chrom, sample + '.' + chrom + '.wg.vcf.gz'))
    return res

rule deepvariant_sample_done:
    input:
        get_deepvariant_files
    output:
        touch(os.path.join(config['DEEPVARIANT'], "{sample}.done"))



def get_chrom_capture_kit_bed(wildcards):
    chrom = wildcards.chrom
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        capture_kit_chr_path_bed = []
    else:
        capture_kit_chr_path_bed = [os.path.join(config['RES'], config['kit_folder'], config['MERGED_CAPTURE_KIT'], config['MERGED_CAPTURE_KIT'] + '_chrom_' + chrom + '.bed')]
    return capture_kit_chr_path_bed

def get_sequencing_mode(wildcards):
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        mode = 'WGS'
    else:
        mode = 'WES'
    return mode        

def get_mem_mb_deepvariant(wildcards, attempt):
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        res = 3500
    else:
        res = 3000

    return (attempt - 1) * 0.5 * res + res

rule deepvariant:
    input:
        bam =  rules.markdup.output.mdbams,
        bai = rules.markdup.output.mdbams_bai,
        bed = get_chrom_capture_kit_bed
    output:
        vcf = temp(os.path.join(config['DEEPVARIANT'],'VCF', "{chrom}","{sample}.{chrom}.vcf.gz")),
        vcf_tbi = temp(os.path.join(config['DEEPVARIANT'],'VCF', "{chrom}","{sample}.{chrom}.vcf.gz.tbi")),
        gvcf = temp(os.path.join(config['DEEPVARIANT'],'gVCF', "{chrom}","{sample}.{chrom}.g.vcf.gz")),
        gvcf_tbi = temp(os.path.join(config['DEEPVARIANT'],'gVCF', "{chrom}","{sample}.{chrom}.g.vcf.gz.tbi")),
        inter_dir = temp(directory(os.path.join(config['DEEPVARIANT'],'DV_intermediate', "{sample}.{chrom}")))
    params: 
            interval = lambda wildcards, input: wildcards['chrom'].replace('YH', 'Y').replace('XH','X') if len(input.bed) == 0 else input.bed[0],
            cd = current_dir + '/',
            mode=get_sequencing_mode
    container: 'docker://google/deepvariant:1.5.0'
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{chrom}.wholedeepvariant.txt")
    resources:
        n="1.5", #set in profile using singularity-args. Waiting for rule-specific args. 
        nshards=2,
        mem_mb=get_mem_mb_deepvariant 
    log: os.path.join(config['LOG'],"{sample}.{chrom}.wholedeepvariant.log")
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant --call_variants_extra_args config_string="device_count:{{key:'CPU' value:4}} inter_op_parallelism_threads:4 intra_op_parallelism_threads:4" --num_shards={resources.nshards} --model_type={params.mode} --regions={params.interval} --ref={ref} --reads={params.cd}{input.bam} --output_vcf={output.vcf} --output_gvcf={output.gvcf} --intermediate_results_dir "{output.inter_dir}"  2> {log}
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
        vcf = os.path.join(config['DEEPVARIANT'], "VCF/{chrom}/{sample}.{chrom}.w.vcf.gz"),
        wstats = os.path.join(config['STAT'], "whatshap_dvphasing/{sample}.{chrom}.stats"),
        mwstats = os.path.join(config['STAT'], "whatshap_dvphasing/{sample}.{chrom}.merge_stats"),
        tmp_gvcf= temp(os.path.join(config['DEEPVARIANT'], "gVCF/{chrom}/{sample}.{chrom}.wg.vcf")), 
        gvcf= os.path.join(config['DEEPVARIANT'], "gVCF/{chrom}/{sample}.{chrom}.wg.vcf.gz"), 
        gvcf_tbi = os.path.join(config['DEEPVARIANT'], "gVCF/{chrom}/{sample}.{chrom}.wg.vcf.gz.tbi"), 
    params:
        merge_script=srcdir("scripts/merge_phasing.py")
    resources: 
        n=1,
        mem_mb = 600
    conda: "envs/whatshap.yaml"
    shell: """
        mkdir -p `dirname {output.wstats}`
        whatshap phase  --ignore-read-groups --reference {ref} {input.vcf} {input.bams} -o {output.vcf}        
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





#         
# rule make_examples:
#     input:
#         bam =  os.path.join(current_dir,rules.markdup.output.mdbams),
#         bai = os.path.join(current_dir,rules.markdup_index.output.mdbams_bai)
#     output:
#         pregvcf = os.path.join(current_dir, config['DEEPVARIANT'],'intermediate', "{chrom}.{sample}.{mode}.gvcf.trfrecord.gz"),
#         examples = os.path.join(current_dir, config['DEEPVARIANT'],'intermediate', "{chrom}.{sample}.{mode}.examples.gz"),
#     threads: config['make_examples']['n']
#     params: inter_dir = os.path.join(current_dir, config['DEEPVARIANT'],'intermediate', "{chrom}.{sample}.{mode}"),
#             # change to merged bed capture kit divided per BINs
#             bed= get_chrom_capture_kit_bed
#     container: 'docker://google/deepvariant:1.4.0'
#     benchmark:
#         os.path.join(current_dir, config['BENCH'],"{chrom}.{sample}.{mode}.makeexamples.txt")
#     log: os.path.join(current_dir, config['LOG'],"{chrom}.{sample}.{mode}.makeexamples.log")
#     shell:
#         # singularity run -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:"1.4.0" \
#         """
#         /opt/deepvariant/bin/make_examples --ref={ref} --reads={input.bam} --gvcf={output.pregvcf} --mode calling --regions {input.bed} --examples={output.examples} --channels="insert_size" --gvcf_gq_binsize 1 2> {log}
#          """
#
# def get_model(wildcards):
#     mode = wildcards.mode
#     if mode == 'WES':
#         model = os.path.join("/opt/models/wes/model.ckpt")
#     elif mode == 'WGS':
#         model = os.path.join("/opt/models/wgs/model.ckpt")
#     return model
#
# rule call_variants:
#     input:
#         examples =  os.path.join(current_dir,rules.make_examples.output.examples),
#     output:
#         precall = os.path.join(current_dir, config['DEEPVARIANT'],'intermediate', "{chrom}.{sample}.{mode}.intermediate_out.gz"),
#     params: model = get_model
#     container: 'docker://google/deepvariant:1.4.0'
#     threads: config['call_variants']['n']
#     benchmark:
#         os.path.join(current_dir, config['BENCH'],"{chrom}_{sample}.{mode}.call_variants.txt")
#     log: os.path.join(current_dir, config['LOG'],"{chrom}_{sample}.{mode}.call_variants.log")
#     shell:
#         """
#         /opt/deepvariant/bin/call_variants --examples {input.examples} --outfile {output.precall} --checkpoint {params.model} 2> {log}
#         """
#
# rule postprocess_variants:
#     input:
#         non_variant_records =  os.path.join(current_dir,rules.make_examples.output.pregvcf),
#         precall= os.path.join(current_dir,rules.call_variants.output.precall),
#     output:
#         vcf = os.path.join(current_dir, config['DEEPVARIANT'],'VCF', "{chrom}.{sample}.{mode}.vcf.gz"),
#         gvcf = os.path.join(current_dir, config['DEEPVARIANT'],'gVCF', "{chrom}.{sample}.{mode}.g.vcf.gz"),
#     params: model = get_model
#     container: 'docker://google/deepvariant:1.4.0'
#     benchmark:
#         os.path.join(current_dir, config['BENCH'],"{chrom}_{sample}.{mode}.postprocess_variants.txt")
#     log: os.path.join(current_dir, config['LOG'],"{chrom}_{sample}.{mode}.postprocess_variants.log")
#     threads: config['postprocess_variants']['n']
#     shell:
#             """
#             /opt/deepvariant/bin/postprocess_variants --ref={ref}  --infile {input.precall} --outfile {output.vcf} --gvcf_outfile {output.gvcf} --nonvariant_site_tfrecord_path {input.non_variant_records} 2> {log}
#             """
