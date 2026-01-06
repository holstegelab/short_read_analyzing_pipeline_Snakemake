import pandas as pd
import read_stats
import os
configfile: srcdir("Snakefile.cluster.json")

gatk = config['gatk']
samtools = config['samtools']
bcftools = config['bcftools']
dragmap = config['dragmap']
verifybamid2 = config['verifybamid2']

ref = config['RES'] + config['ref']

wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    readgroup=r"[\w\d_\-@]+"


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)
# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
module gVCF:
    snakefile: 'gVCF.smk'
    config: config
module DBImport:
    snakefile: 'DBImport.smk'
    config: config
module Genotype:
    snakefile: 'Genotype.smk'
    config: config

use rule * from Genotype

mode = config.get("computing_mode", "WES")
rule VQSR_all:
    input:
        rules.Genotype_all.input,
        expand("{vcf}/VQSR/{chr}/Merged_after_VQSR_{chr}_{mode}.vcf",vcf=config['VCF_Final'], chr = main_chr, mode = mode),
    default_target: True



# VQSR
#select SNPs for VQSR
# SNPs and INDELs require different options
rule SelectSNPs:
    input:
        vcf = rules.merge_per_chr.output.per_chr_vcfs,
        tbi= rules.index_per_chr.output.vcfidx
    output:
        SNP_vcf = temp(config['VCF_Final'] + "/{chr}/Merged_SNPs_{chr}_{mode}.vcf")
        # SNP_vcf=temp(config['VCF'] + "/Merged_SNPs.vcf")
    priority: 56
    log: config['LOG'] + '/' + "SelectSNPs_{chr}_{mode}.log"
    benchmark: config['BENCH'] + "/SelectSNPs_{chr}_{mode}.txt"
    conda: CONDA_VCF
    shell:
        """
        {gatk} SelectVariants \
                --select-type-to-include SNP \
                -V {input.vcf} -O {output} 2> {log}
        """

# 1st step of VQSR - calculate scores
rule VQSR_SNP:
    input:
        rules.SelectSNPs.output.SNP_vcf
    output:
        recal_snp=temp(config['VCF_Final'] + "/VQSR/{chr}/SNPs_vqsr_{chr}_{mode}.recal"),
        tranches_file_snp=temp(config['VCF_Final'] + "/VQSR/{chr}/SNPs_vqsr_{chr}_{mode}.tranches"),
        # recal_snp=temp(config['VCF'] + "/SNPs_vqsr.recal"),
        # tranches_file_snp=temp(config['VCF'] + "/SNPs_vqsr.tranches"),
        r_snp=config['STAT'] + "/VQSR/R/{chr}/SNPs_vqsr_plots_{mode}.R"
    log: config['LOG'] + '/' + "VQSR_SNP_{chr}_{mode}.log"
    benchmark: config['BENCH'] + "/VQSR_SNP_{chr}_{mode}.txt"
    params:
        hapmap = config['RES'] + config['hapmap'],
        omni = config['RES'] + config['omni'],
        kilo_g = config['RES'] + config['kilo_g'],
        dbsnp = config['RES'] + config['dbsnp']
    conda: CONDA_VCF
    priority: 58
    shell:
        # -an InbreedingCoeff if 10+
        """
        {gatk} VariantRecalibrator\
         -R {ref} -V {input} \
         -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
         -resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
         -resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.kilo_g} \
         -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
         -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP --trust-all-polymorphic -AS TRUE -an InbreedingCoeff \
         -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
         -O {output.recal_snp} \
         --tranches-file {output.tranches_file_snp} \
         --rscript-file {output.r_snp} 2> {log}
        """

rule ApplyVQSR_SNPs:
    input:
        recal_snp=rules.VQSR_SNP.output.recal_snp,
        tranches_snp=rules.VQSR_SNP.output.tranches_file_snp,
        snps_variants=rules.SelectSNPs.output.SNP_vcf
    output:
        recal_vcf_snp = (config['VCF_Final'] + "/VQSR/{chr}/SNPs_recal_apply_vqsr_{chr}_{mode}.vcf.gz"),
        tbi = ensure((config['VCF_Final'] + "/VQSR/{chr}/SNPs_recal_apply_vqsr_{chr}_{mode}.vcf.gz.tbi"), non_empty=True)
        # recal_vcf_snp=temp(config['VCF'] + "/SNPs_recal_apply_vqsr.vcf")
    log: config['LOG'] + '/' + "Apply_VQSR_SNP_{chr}_{mode}.log"
    benchmark: config['BENCH'] + "/Apply_VQSR_SNP_{chr}_{mode}.txt"
    params:
        ts_level='99.0'  #ts-filter-level show the "stregnth" of VQSR could be from 90 to 100
    priority: 60
    conda: CONDA_VCF
    shell:
        """
        {gatk} ApplyVQSR -R {ref} -mode SNP \
        --recal-file {input.recal_snp} --tranches-file {input.tranches_snp} \
        -O {output.recal_vcf_snp} -V {input.snps_variants} -ts-filter-level {params.ts_level} -AS TRUE --create-output-variant-index true 2> {log}
        """

# select INDELs for VQSR
rule SelectINDELs:
    input:
        vcf = rules.merge_per_chr.output.per_chr_vcfs,
        tbi= rules.index_per_chr.output.vcfidx
    output:
        INDEL_vcf=temp(config['VCF_Final'] + "/{chr}/Merged_INDELs_{chr}_{mode}.vcf")
    log: config['LOG'] + '/' + "SelectINDELS_{chr}_{mode}.log"
    benchmark: config['BENCH'] + "/SelectINDELs_{chr}_{mode}.txt"
    priority: 56
    conda: CONDA_VCF
    shell:
        """
        {gatk} SelectVariants \
                --select-type-to-include INDEL \
                -V {input.vcf} -O {output} 2> {log}
        """

# 1st step of VQSR - calculate scores
rule VQSR_INDEL:
    input:
        rules.SelectINDELs.output.INDEL_vcf
    output:
        recal_indel=temp(config['VCF_Final'] + "/VQSR/{chr}/INDELs_vqsr_{chr}_{mode}.recal"),
        tranches_file_indel=temp(config['VCF_Final'] + "/VQSR/{chr}/INDELs_vqsr_{chr}_{mode}.tranches"),
        r_indel=config['STAT'] + "/VQSR/R/INDELs_vqsr_plots_{chr}_{mode}.R"
    log: config['LOG'] + '/' + "VQSR_INDEL_{chr}_{mode}.log"
    benchmark: config['BENCH'] + "/VQSR_INDEL_{chr}_{mode}.txt"
    priority: 58
    params:
        mills = config['RES'] + config['mills'],
        dbsnp_indel = config['RES'] + config['dbsnp_indel']
    conda: CONDA_VCF
    shell:
        # -an InbreedingCoeff if 10+
        """
        {gatk} VariantRecalibrator -R {ref} -V {input} \
        -O {output.recal_indel} --tranches-file {output.tranches_file_indel} --rscript-file {output.r_indel} \
         --max-gaussians 4 --trust-all-polymorphic -AS TRUE\
         -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
         -resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mills} \
         -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp_indel} \
         -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode INDEL -an InbreedingCoeff \
         --trust-all-polymorphic 2> {log}
        """

rule ApplyVQSR_INDEs:
    input:
        recal_indel=rules.VQSR_INDEL.output.recal_indel,
        tranches_indel=rules.VQSR_INDEL.output.tranches_file_indel,
        indel_variants=rules.SelectINDELs.output.INDEL_vcf
    log: config['LOG'] + '/' + "ApplyVQSR_INDELs_{chr}_{mode}.log"
    benchmark: config['BENCH'] + "/ApplyVQSR_INDELs_{chr}_{mode}.txt"
    output:
        recal_vcf_indel=temp(config['VCF_Final'] + "/VQSR/{chr}/INDELs_recal_apply_vqsr_{chr}_{mode}.vcf")
    params:
        ts_level='97.0'  #ts-filter-level show the "stregnth" of VQSR could be from 90 to 100
    conda: CONDA_VCF
    priority: 60
    shell:
        """
        {gatk} ApplyVQSR -R {ref} -mode INDEL \
        --recal-file {input.recal_indel} --tranches-file {input.tranches_indel} \
        -O {output} -V {input.indel_variants} -ts-filter-level {params.ts_level} -AS TRUE 2> {log}
        """

rule select_norsnp_nor_indels:
    input:
        vcf = rules.merge_per_chr.output.per_chr_vcfs,
        tbi= rules.index_per_chr.output.vcfidx
    output:
        MIX_vcf=temp(config['VCF_Final'] + "/{chr}/merged_mix_{chr}_{mode}.vcf")
    log: config['LOG'] + '/' + "SelectMix_{chr}_{mode}.log"
    benchmark: config['BENCH'] + "/SelectMix_{chr}_{mode}.txt"
    priority: 59
    conda: CONDA_VCF
    shell:
        """
        {gatk} SelectVariants \
                --select-type-to-exclude INDEL --select-type-to-exclude SNP --sites-only-vcf-output TRUE\
                -V {input.vcf} -O {output} 2> {log}
        """


#combine filtr results
rule combine:
    input:
        snps=rules.ApplyVQSR_SNPs.output.recal_vcf_snp,
        indel=rules.ApplyVQSR_INDEs.output.recal_vcf_indel,
        mix = rules.select_norsnp_nor_indels.output.MIX_vcf
    log: config['LOG'] + "/combine_{chr}_{mode}.log"
    benchmark: config['BENCH'] + "/combine_{chr}_{mode}.txt"
    output:
        filtrVCF=(config['VCF_Final'] + "/VQSR/{chr}/Merged_after_VQSR_{chr}_{mode}.vcf")
    priority: 70
    conda: CONDA_VCF
    shell:
        "{gatk} MergeVcfs -I {input.snps} -I {input.indel} -O {output} 2> {log}"

# # normalization with bcftools
# rule norm:
#     input:
#         rules.combine.output.filtrVCF
#     output:CollectVariantCallingMetrics
#         normVCF=config['VCF_Final'] + "/Merged_after_VQSR_norm.vcf",
#         idx=config['VCF_Final'] + "/Merged_after_VQSR_norm.vcf.idx"
#     log: config['LOG'] + '/' + "normalization.log"
#     benchmark: config['BENCH'] + "/normalization.txt"
#     priority: 80
#     conda: "envs/preprocess.yaml"
#     shell:
#         "{bcftools} norm -f {ref} {input} -m -both -O v | {bcftools} norm -d exact -f {ref} > {output.normVCF} 2> {log} && {gatk} IndexFeatureFile -I {output.normVCF} -O {output.idx} "
