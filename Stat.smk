import pandas as pd
import numpy
import read_stats
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
MERGED_CAPTURE_KIT = config['MERGED_CAPTURE_KIT']
ref = config['RES'] + config['ref']

wildcard_constraints:
    sample="[\w\d_\-@]+"


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO = load_samplefiles('.', config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()


tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner
# module VCF:
#     snakefile: 'gVCF.smk'

rule Stat_all:
    input:
        expand("{stats}/{sample}_hs_metrics",sample=sample_names, stats = config['STAT']),
        expand("{stats}/{sample}.OXOG", sample=sample_names, stats = config['STAT']),
        expand("{stats}/{sample}_samtools.stat", sample=sample_names, stats = config['STAT']),
        expand("{samplefile}.oxo_quality.tab", samplefile = SAMPLE_FILES),
        expand("{samplefile}.bam_quality.tab", samplefile = SAMPLE_FILES),
    default_target: True

# return interval_list file instead of bed file
def get_capture_kit_interval_list(wildcards):
    if SAMPLEINFO[wildcards['sample']]['sample_type'].startswith('illumina_wgs'):
        capture_kit = MERGED_CAPTURE_KIT
    else:
        capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    capture_kit_path = os.path.join(config['RES'], config['kit_folder'], capture_kit + '_hg38.interval_list')
    return capture_kit_path

def get_mem_mb_hs_stats(wildcrads, attempt):
    return (attempt * int(config['hs_stats']['mem']))
#hsmetrics
#include off-target metrics
rule hs_stats:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output:
        HS_metrics=os.path.join(config['STAT'], "{sample}_hs_metrics")
    log: os.path.join(config['LOG'], "HS_stats_{sample}.log")
    benchmark: os.path.join(config['BENCH'],  "HS_stats_{sample}.txt")
    priority: 99
    params:
        interval = get_capture_kit_interval_list,
        #minimum Base Quality for a base to contribute cov
        #def is 20
        Q=10,
        #minimin Mapping Quality for a read to contribute cov
        #def is 20
        MQ=10,
    conda: "envs/preprocess.yaml"
    resources: mem_mb = get_mem_mb_hs_stats,
                tmpdir = tmpdir
    shell:
        """
            gatk  CollectHsMetrics --java-options "-Xmx{resources.mem_mb}m"  --TMP_DIR {resources.tmpdir} \
            -I {input.bam} -R {ref} -BI {params.interval} -TI {params.interval} \
            -Q {params.Q} -MQ {params.MQ} \
            --PER_TARGET_COVERAGE stats/{wildcards.sample}_per_targ_cov \
            -O stats/{wildcards.sample}_hs_metrics 2> {log}"""

def get_mem_mb_Artifact_stats(wildcrads, attempt):
    return (attempt * int(config['Artifact_stats']['mem']))

rule Artifact_stats:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output:
        Bait_bias = os.path.join(config['STAT'], '{sample}.bait_bias.bait_bias_summary_metrics'),
        Pre_adapter = os.path.join(config['STAT'], '{sample}.bait_bias.pre_adapter_summary_metrics'),
        Bait_bias_det = os.path.join(config['STAT'],'{sample}.bait_bias.bait_bias_detail_metrics'),
        Pre_adapter_det = os.path.join(config['STAT'], '{sample}.bait_bias.pre_adapter_detail_metrics'),
        # Artifact_matrics = config['STAT'] + "/{sample}.bait_bias.bait_bias_detail_metrics"
    priority: 99
    log: os.path.join(config['LOG'], "Artifact_stats_{sample}.log")
    benchmark: os.path.join(config['BENCH'], "Artifact_stats_{sample}.txt")
    params:
        interval = get_capture_kit_interval_list,
        # output define prefix, not full filename
        # params.out define prefix and output define whole outputs' filename
        out = os.path.join(config['STAT'], "{sample}.bait_bias"),
        dbsnp = os.path.join(config['RES'], config['dbsnp'])
    conda: "envs/preprocess.yaml"
    resources: mem_mb = get_mem_mb_Artifact_stats,
                tmpdir= tmpdir
    shell:
        """
            gatk  CollectSequencingArtifactMetrics --java-options "-Xmx{resources.mem_mb}m" --TMP_DIR {resources.tmpdir} -I {input.bam} -O {params.out} \
        -R {ref} --DB_SNP {params.dbsnp} --INTERVALS {params.interval} 2> log"""

def get_mem_mb_OXOG_metrics(wildcrads, attempt):
    return (attempt * int(config['OXOG_metrics']['mem']))
rule OXOG_metrics:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output:
        Artifact_matrics = os.path.join(config['STAT'], "{sample}.OXOG")
    priority: 99
    log: os.path.join(config['LOG'], "OXOG_stats_{sample}.log")
    benchmark: os.path.join(config['BENCH'], "OxoG_{sample}.txt")
    params:
        interval = get_capture_kit_interval_list,
        dbsnp = os.path.join(config['RES'], config['dbsnp'])
    conda: "envs/preprocess.yaml"
    resources:
                mem_mb = get_mem_mb_OXOG_metrics,
                tmpdir= tmpdir
    shell:
        """gatk CollectOxoGMetrics --java-options "-Xmx{resources.mem_mb}m" -I {input.bam} -O {output} -R {ref} \
         --DB_SNP {params.dbsnp} --INTERVALS {params.interval} 2> {log}"""

rule samtools_stat:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output: samtools_stat = os.path.join(config['STAT'], "{sample}_samtools.stat")
    priority: 99
    log: os.path.join(config['LOG'], "samtools_{sample}.log")
    benchmark: os.path.join(config['BENCH'], "samtools_stat_{sample}.txt")
    threads: config['samtools_stat']['n']
    conda: "envs/preprocess.yaml"
    shell:
        "samtools stat -@ {threads} -r {ref} {input.bam} > {output}"

# extract info about capture kit from SAMPLEFILE
# assume that all kits bed and interval_list files are existing and download to res folder
def get_capture_kit_bed(wildcards):
    if SAMPLEINFO[wildcards['sample']]['sample_type'] == 'illumina_wgs':
        capture_kit = MERGED_CAPTURE_KIT
    elif SAMPLEINFO[wildcards['sample']]['capture_kit'] == '':
        capture_kit = MERGED_CAPTURE_KIT
    else:
        capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    #
    # if capture_kit.strip() == '':
    #     capture_kit = os.path.basename(MERGED_CAPTURE_KIT)[:-4]
    capture_kit_path = os.path.join(config['RES'], config['kit_folder'], f'{capture_kit}_hg38.bed')
    return capture_kit_path

rule samtools_stat_exome:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output: samtools_stat_exome = os.path.join(config['STAT'], "{sample}_samtools.exome.stat")
    priority: 99
    params:
        bed_interval = get_capture_kit_bed
    log: os.path.join(config['LOG'], "samtools_exome_{sample}.log")
    benchmark: os.path.join(config['BENCH'],"samtools_stat_exome_{sample}.txt")
    threads: config['samtools_stat']['n']
    conda: "envs/preprocess.yaml"
    shell:
        "samtools stat -@ {threads} -t {params.bed_interval} -r {ref} {input.bam} > {output}"

rule bamstats_all:
    input:
        # very annoying bug here
        # if in input 2 functions and one of them is check chekpoint (rules.markdup.output.mdbams and get_capture_kit_bed here was as example)
        # first command has not been executed
        # and in shell wildcard (instead of iutput of function) has putted
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output:
        All_exome_stats = os.path.join(config['STAT'],'{sample}.bam_all.tsv')
    benchmark: os.path.join(config['BENCH'],"bamstats_all_{sample}.txt")
    threads: config['bamstats_all']['n']
    params:
        py_stats = srcdir(config['BAMSTATS'])
    conda: "envs/pypy.yaml"
    shell:
        "samtools view -s 0.05 -h {input.bam} --threads {threads}  | pypy {params.py_stats} stats > {output}"

rule bamstats_exome:
    input:
        # very annoying bug here
        # if in input 2 functions and one of them is check chekpoint (rules.markdup.output.mdbams and get_capture_kit_bed here was as example)
        # first command has not been executed
        # and in shell wildcard (instead of iutput of function) has putted
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output:
        All_exome_stats = os.path.join(config['STAT'], '{sample}.bam_exome.tsv')
    benchmark: os.path.join(config['BENCH'],"bamstats_exome_{sample}.txt")
    threads: config['bamstats_exome']['n']
    params:
        py_stats = srcdir(config['BAMSTATS']),
        bed_interval= get_capture_kit_bed,
    conda: "envs/pypy.yaml"
    shell:
        "samtools view -s 0.05 -h {input.bam} --threads {threads} -L {params.bed_interval} | pypy {params.py_stats} stats > {output}"


def get_quality_stats(wildcards):
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    files = []
    samples = list(sampleinfo.keys())
    samples.sort()
    for sample in samples:
        files.append(os.path.join(config['STAT'],  f"{sample}_samtools.stat"))
        files.append(os.path.join(config['STAT'],  f'{sample}_samtools.exome.stat'))
        files.append(os.path.join(config['STAT'], 'contam',  f'{sample}_verifybamid.pca2.selfSM'))
        files.append(os.path.join(config['STAT'], f'{sample}.bam_all.tsv'))
        files.append(os.path.join(config['STAT'], f'{sample}.bam_exome.tsv'))
        files.append(os.path.join(config['STAT'], f'{sample}.bait_bias.pre_adapter_summary_metrics'))
        files.append(os.path.join(config['STAT'], f'{sample}.bait_bias.bait_bias_summary_metrics'))
    return files

rule gatherstats:
    # keep in mind, that samtools_stat create file even if it it's finished with error or you force to stop it
    # if you force to stop samtools_stat delete all output to prevent errors
    # rm -r stats/*samtools*
    input:
        get_quality_stats
    benchmark: os.path.join(config['BENCH'],"{samplefile}_gatherstat.txt")
    output:
        '{samplefile}.bam_quality.tab'
    run:
        sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
        samples = list(sampleinfo.keys())
        samples.sort()

        rinput = numpy.array(input).reshape(len(samples),int(len(input) / len(samples)))
        stats, exome_stats, vpca2, bam_extra_all, bam_extra_exome, pre_adapter, bait_bias = zip(*rinput)

        header, data = read_stats.combine_quality_stats(samples,stats,exome_stats,vpca2,bam_extra_all,bam_extra_exome,pre_adapter,bait_bias)
        read_stats.write_tsv(str(output),header,data)

def get_oxo_stats(wildcards):
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    files = []
    samples = list(sampleinfo.keys())
    samples.sort()
    for sample in samples:
        files.append(os.path.join(config['STAT'], f'{sample}.bait_bias.pre_adapter_detail_metrics'))
        files.append(os.path.join(config['STAT'],f'{sample}.bait_bias.bait_bias_detail_metrics'))
    return files

rule gatherosostats:
    input:
        get_oxo_stats
    output:
        '{samplefile}.oxo_quality.tab'
    benchmark: os.path.join(config['BENCH'],"{samplefile}_gatherOXOstat.txt")
    run:
        sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
        samples = list(sampleinfo.keys())
        samples.sort()

        rinput = numpy.array(input).reshape(len(samples),int(len(input) / len(samples)))
        pre_adapter, bait_bias = zip(*rinput)

        header, data = read_stats.combine_oxo_stats(samples,pre_adapter,bait_bias)
        read_stats.write_tsv(str(output),header,data)

