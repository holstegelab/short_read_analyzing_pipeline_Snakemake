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

mode = config.get("computing_mode", "WES")

rule Deepvariant_all:
    input:
        expand("{cur_dir}/{deepvarinat}/gVCF/{chr}/{chr}.{sample}.{mode}.g.vcf.gz", cur_dir = current_dir, deepvarinat = config['DEEPVARIANT'], sample = sample_names, mode = mode, chr = main_chrs),
        rules.Aligner_all.input
    default_target: True

def get_chrom_capture_kit_bed(wildcards):
    chr = wildcards.chr
    if mode == 'WES':
        capture_kit_chr_path_bed = config['RES'] + config['kit_folder'] + config['MERGED_CAPTURE_KIT'] + '_hg38/' + config['MERGED_CAPTURE_KIT'] + '_hg38_' + chr + '.interval_list.bed'
    elif mode == 'WGS':
        capture_kit_chr_path_bed = chr
    return capture_kit_chr_path_bed

rule deepvariant:
    input:
        bam =  rules.markdup.output.mdbams,
        bai = rules.markdup_index.output.mdbams_bai
    output:
        vcf = os.path.join(current_dir, config['DEEPVARIANT'],'VCF', "{chr}","{chr}.{sample}.{mode}.vcf.gz"),
        gvcf = os.path.join(current_dir, config['DEEPVARIANT'],'gVCF', "{chr}","{chr}.{sample}.{mode}.g.vcf.gz"),
    threads: config['deepvariant']['n']
    params: inter_dir = os.path.join(current_dir, config['DEEPVARIANT'],'DV_intermediate', "{chr}.{sample}.{mode}"),
            bed=get_chrom_capture_kit_bed,
            cd = current_dir + '/'
    #conda: "envs/deepvariant.yaml"
    container: 'docker://google/deepvariant:1.4.0'
    benchmark:
        os.path.join(current_dir, config['BENCH'],"{chr}_{sample}.{mode}.wholedeepvariant.txt")
    log: os.path.join(current_dir, config['LOG'],"{chr}_{sample}.{mode}.wholedeepvariant.log")
    shell:
        # singularity run -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:"1.4.0" \
        """
        /opt/deepvariant/bin/run_deepvariant --model_type={wildcards.mode} --regions={params.bed} --ref={ref} --reads={params.cd}{input.bam} --output_vcf={output.vcf} --output_gvcf={output.gvcf} --intermediate_results_dir "{params.inter_dir}" --num_shards={threads} 2> {log}
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
#         pregvcf = os.path.join(current_dir, config['DEEPVARIANT'],'intermediate', "{chr}.{sample}.{mode}.gvcf.trfrecord.gz"),
#         examples = os.path.join(current_dir, config['DEEPVARIANT'],'intermediate', "{chr}.{sample}.{mode}.examples.gz"),
#     threads: config['make_examples']['n']
#     params: inter_dir = os.path.join(current_dir, config['DEEPVARIANT'],'intermediate', "{chr}.{sample}.{mode}"),
#             # change to merged bed capture kit divided per BINs
#             bed= get_chrom_capture_kit_bed
#     container: 'docker://google/deepvariant:1.4.0'
#     benchmark:
#         os.path.join(current_dir, config['BENCH'],"{chr}.{sample}.{mode}.makeexamples.txt")
#     log: os.path.join(current_dir, config['LOG'],"{chr}.{sample}.{mode}.makeexamples.log")
#     shell:
#         # singularity run -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:"1.4.0" \
#         """
#         /opt/deepvariant/bin/make_examples --ref={ref} --reads={input.bam} --gvcf={output.pregvcf} --mode calling --regions {params.bed} --examples={output.examples} --channels="insert_size" --gvcf_gq_binsize 1 2> {log}
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
#         precall = os.path.join(current_dir, config['DEEPVARIANT'],'intermediate', "{chr}.{sample}.{mode}.intermediate_out.gz"),
#     params: model = get_model
#     container: 'docker://google/deepvariant:1.4.0'
#     threads: config['call_variants']['n']
#     benchmark:
#         os.path.join(current_dir, config['BENCH'],"{chr}_{sample}.{mode}.call_variants.txt")
#     log: os.path.join(current_dir, config['LOG'],"{chr}_{sample}.{mode}.call_variants.log")
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
#         vcf = os.path.join(current_dir, config['DEEPVARIANT'],'VCF', "{chr}.{sample}.{mode}.vcf.gz"),
#         gvcf = os.path.join(current_dir, config['DEEPVARIANT'],'gVCF', "{chr}.{sample}.{mode}.g.vcf.gz"),
#     params: model = get_model
#     container: 'docker://google/deepvariant:1.4.0'
#     benchmark:
#         os.path.join(current_dir, config['BENCH'],"{chr}_{sample}.{mode}.postprocess_variants.txt")
#     log: os.path.join(current_dir, config['LOG'],"{chr}_{sample}.{mode}.postprocess_variants.log")
#     threads: config['postprocess_variants']['n']
#     shell:
#             """
#             /opt/deepvariant/bin/postprocess_variants --ref={ref}  --infile {input.precall} --outfile {output.vcf} --gvcf_outfile {output.gvcf} --nonvariant_site_tfrecord_path {input.non_variant_records} 2> {log}
#             """
