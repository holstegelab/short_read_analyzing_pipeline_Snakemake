import pandas as pd
import read_stats
import os
configfile: "Snakefile.cluster.json"
configfile: "Snakefile.paths.yaml"
gatk = config['miniconda'] + config['gatk']
samtools = config['miniconda'] + config['samtools']
bcftools = config['miniconda'] + config['bcftools']
dragmap = config['miniconda'] + config['dragmap']
cutadapt = config['miniconda'] + config['cutadapt']
verifybamid2 = config['miniconda'] + config['verifybamid2']
ref = config['RES'] + config['ref']

wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+"

# main chromosomes from GRCh38 splitted into 99 bins
# bins = config['RES'] + config['bin_file_ref']
# import csv
# chrs = []
# with open(bins) as file:
#     tsv_file = csv.reader(file, delimiter="\t")
#     for line in tsv_file:
#         if line[1] == str('0'):
#             out = line[0] + ':' + str('1') + '-' + line[2]
#         else:
#             out = line[0] + ':' + line[1] + '-' + line[2]
#         chrs.append(out)
# print(chrs)

from read_samples import *
from common import *

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

rule all:
    input:
        config['STAT'] + "/BASIC.variant_calling_detail_metrics",
        expand("{stats}/{sample}_hs_metrics",sample=sample_names, stats = config['STAT']),
        expand("{stats}/{sample}.OXOG", sample=sample_names, stats = config['STAT']),
        expand("{stats}/{sample}_samtools.stat", sample=sample_names, stats = config['STAT']),
        expand("{samplefile}.oxo_quality.tab", samplefile = SAMPLE_FILES),
        expand("{samplefile}.bam_quality.tab", samplefile = SAMPLE_FILES),



# basic stats
# include hom-het ratio, titv ratio, etc.
rule Basic_stats:
    input:
        rules.norm.output.normVCF
    output:
        config['STAT'] + "/BASIC.variant_calling_detail_metrics",
        config['STAT'] + "/BASIC.variant_calling_summary_metrics"
    priority: 90
    log: config['LOG'] + '/' + "VCF_stats.log"
    benchmark: config['BENCH'] + "/VCF_stats.txt"
    params: dbsnp = config['RES'] + config['dbsnp']
    threads: config['Stat_Basic_stats']['n']
    shell:
        "{gatk} CollectVariantCallingMetrics \
        -R {ref} -I {input} -O stats/BASIC \
        --DBSNP {params.dbsnp} --THREAD_COUNT {threads} 2> {log}"

# return interval_list file instead of bed file
def get_capture_kit_interval_list(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    capture_kit_path = config['RES'] + config['kit_folder'] + capture_kit + '_hg38.interval_list'
    return capture_kit_path


#hsmetrics
#include off-target metrics
rule HS_stats:
    input:
        bam = rules.CalibrateDragstrModel.input.bam
    output:
        HS_metrics=config['STAT'] + "/{sample}_hs_metrics"
    log: config['LOG'] + '/' + "HS_stats_{sample}.log"
    benchmark: config['BENCH'] + "/HS_stats_{sample}.txt"
    priority: 99
    params:
        interval = get_capture_kit_interval_list,
        #minimum Base Quality for a base to contribute cov
        #def is 20
        Q=10,
        #minimin Mapping Quality for a read to contribute cov
        #def is 20
        MQ=10,
    shell:
        "{gatk} CollectHsMetrics \
            -I {input.bam} -R {ref} -BI {params.interval} -TI {params.interval} \
            -Q {params.Q} -MQ {params.MQ} \
            --PER_TARGET_COVERAGE stats/{wildcards.sample}_per_targ_cov \
            -O stats/{wildcards.sample}_hs_metrics 2> {log}"

rule Artifact_stats:
    input:
        bam = rules.CalibrateDragstrModel.input.bam
    output:
        Bait_bias = config['STAT'] + '/{sample}.bait_bias.bait_bias_summary_metrics',
        Pre_adapter = config['STAT'] + '/{sample}.bait_bias.pre_adapter_summary_metrics',
        Bait_bias_det = config['STAT'] + '/{sample}.bait_bias.bait_bias_detail_metrics',
        Pre_adapter_det = config['STAT'] + '/{sample}.bait_bias.pre_adapter_detail_metrics',
        # Artifact_matrics = config['STAT'] + "/{sample}.bait_bias.bait_bias_detail_metrics"
    priority: 99
    log: config['LOG'] + '/' + "Artifact_stats_{sample}.log"
    benchmark: config['BENCH'] + "/Artifact_stats_{sample}.txt"
    params:
        interval = get_capture_kit_interval_list,
        # output define prefix, not full filename
        # params.out define prefix and output define whole outputs' filename
        out = config['STAT'] + "/{sample}.bait_bias",
        dbsnp = config['RES'] + config['dbsnp']
    shell:
        "{gatk} CollectSequencingArtifactMetrics -I {input.bam} -O {params.out} \
        -R {ref} --DB_SNP {params.dbsnp} --INTERVALS {params.interval} 2> log"

rule OXOG_metrics:
    input:
        bam = rules.CalibrateDragstrModel.input.bam,
    output:
        Artifact_matrics = config['STAT'] + "/{sample}.OXOG"
    priority: 99
    log: config['LOG'] + '/' + "OXOG_stats_{sample}.log"
    benchmark: config['BENCH'] + "/OxoG_{sample}.txt"
    params:
        interval = get_capture_kit_interval_list,
        dbsnp = config['RES'] + config['dbsnp']
    shell:
        "{gatk} CollectOxoGMetrics -I {input.bam} -O {output} -R {ref} \
         --DB_SNP {params.dbsnp} --INTERVALS {params.interval} 2> {log}"

rule samtools_stat:
    input:
        bam = rules.CalibrateDragstrModel.input.bam
    output: samtools_stat = config['STAT'] + "/{sample}_samtools.stat"
    priority: 99
    log: config['LOG'] + '/' + "samtools_{sample}.log"
    benchmark: config['BENCH'] + "/samtools_stat_{sample}.txt"
    threads: config['Stat_samtools_stat']['n']
    shell:
        "{samtools} stat -@ {threads} -r {ref} {input.bam} > {output}"

# extract info about capture kit from SAMPLEFILE
# assume that all kits bed and interval_list files are existing and download to res folder
def get_capture_kit_bed(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    # if capture_kit.strip() == '':
    #     capture_kit = os.path.basename(MERGED_CAPTURE_KIT)[:-4]
    capture_kit_path = config['RES'] + config['kit_folder'] + capture_kit + '_hg38.bed'
    return capture_kit_path

rule samtools_stat_exome:
    input:
        bam = rules.CalibrateDragstrModel.input.bam,
    output: samtools_stat_exome = config['STAT'] + "/{sample}_samtools.exome.stat"
    priority: 99
    params:
        bed_interval = get_capture_kit_bed
    log: config['LOG'] + '/' + "samtools_exome_{sample}.log"
    benchmark: config['BENCH'] + "/samtools_stat_exome_{sample}.txt"
    threads: config['Stat_samtools_stat']['n']
    shell:
        "{samtools} stat -@ {threads} -t {params.bed_interval} -r {ref} {input.bam} > {output}"

rule bamstats_exome:
    input:
        # very annoying bug here
        # if in input 2 functions and one of them is check chekpoint (check_supp and get_capture_kit_bed here was as example)
        # first command has not been executed
        # and in shell wildcard (instead of iutput of function) has putted
        bam = rules.CalibrateDragstrModel.input.bam,
    output:
        All_exome_stats = config['STAT'] + '/{sample}.bam_exome.tsv'
    threads: config['Stat_bamstats_exome']['n']
    params:
        py_stats = config['BAMSTATS'],
        bed_interval= get_capture_kit_bed,
    shell:
        "samtools view -s 0.05 -h {input.bam} --threads {threads} -L {params.bed_interval} | python3 {params.py_stats} stats > {output}"

def get_quality_stats(wildcards):
    sampleinfo = SAMPLES_BY_FILE[os.path.basename(wildcards['samplefile'])]
    files = []
    samples = list(sampleinfo.keys())
    samples.sort()
    for sample in samples:
        files.append(config['STAT'] + '/' + sample + "_samtools.stat")
        files.append(config['STAT'] + '/' + sample + '_samtools.exome.stat')
        files.append(config['STAT'] + '/contam/' + sample + '_verifybamid.pca2.selfSM')
        files.append(check_supp_stats(sample))
        files.append(config['STAT'] + '/' + sample + '.bam_exome.tsv')
        files.append(config['STAT'] + '/' + sample + '.bait_bias.pre_adapter_summary_metrics')
        files.append(config['STAT'] + '/' + sample + '.bait_bias.bait_bias_summary_metrics')
    return files

rule gatherstats:
    # keep in mind, that samtools_stat create file even if it it's finished with error or you force to stop it
    # if you force to stop samtools_stat delete all output to prevent errors
    # rm -r stats/*samtools*
    input:
        get_quality_stats
    output:
        '{samplefile}.bam_quality.tab'
    run:
        sampleinfo = SAMPLES_BY_FILE[os.path.basename(wildcards['samplefile'])]
        samples = list(sampleinfo.keys())
        samples.sort()

        rinput = numpy.array(input).reshape(len(samples),int(len(input) / len(samples)))
        stats, exome_stats, vpca2, bam_extra_all, bam_extra_exome, pre_adapter, bait_bias = zip(*rinput)

        header, data = read_stats.combine_quality_stats(samples,stats,exome_stats,vpca2,bam_extra_all,bam_extra_exome,pre_adapter,bait_bias)
        read_stats.write_tsv(str(output),header,data)

def get_oxo_stats(wildcards):
    sampleinfo = SAMPLES_BY_FILE[os.path.basename(wildcards['samplefile'])]
    files = []
    samples = list(sampleinfo.keys())
    samples.sort()
    for sample in samples:
        files.append(config['STAT'] + '/' + sample + '.bait_bias.pre_adapter_detail_metrics')
        files.append(config['STAT'] + '/' + sample + '.bait_bias.bait_bias_detail_metrics')
    return files

rule gatherosostats:
    input:
        get_oxo_stats
    output:
        '{samplefile}.oxo_quality.tab'
    run:
        sampleinfo = SAMPLES_BY_FILE[os.path.basename(wildcards['samplefile'])]
        samples = list(sampleinfo.keys())
        samples.sort()

        rinput = numpy.array(input).reshape(len(samples),int(len(input) / len(samples)))
        pre_adapter, bait_bias = zip(*rinput)

        header, data = read_stats.combine_oxo_stats(samples,pre_adapter,bait_bias)
        read_stats.write_tsv(str(output),header,data)

