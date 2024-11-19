import pandas as pd
import read_stats
import os
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
    sample=r"[\w\d_\-@]+",
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

# find a variant to determine samples with close coverage profile
# cov2bed


rule CNV_with_cnvkit_Module_all:
    input:
        rules.Aligner_all.input,
        config['CNVKIT'] + '/stats/Sex_sample_table.tsv',
        expand(config['CNVKIT'] + '/plots/diagram/{sample}-diagram.pdf', sample=sample_names),
        expand(config['CNVKIT'] + '/plots/scatter/{sample}-scatter.pdf', sample=sample_names),
        # expand('{cnvkit}/plots/scatter/{sample}-scatter.pdf', sample=sample_names, cnvkit = config['CNVKIT']),
        expand('{cnvkit}/VCF/{sample}_cnv.vcf', sample=sample_names, cnvkit = config['CNVKIT']),
        expand('{cnvkit}/call/{sample}.cal.cns', sample = sample_names, cnvkit = config['CNVKIT']),
        expand('{cnvkit}/stats/genemetric/{sample}_genemetric.stat', sample = sample_names, cnvkit = config['CNVKIT']),
        expand('{cnvkit}/stats/Metrics_table.tsv', cnvkit = config['CNVKIT']),
        expand('{cnvkit}/breaks/{sample}_breaks', sample = sample_names, cnvkit = config['CNVKIT'])
    default_target: True

# bams that align only on "main" chromosomes
rule neewBams:
    input:
        bam_files = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output:
        NBams = config['BAM'] + '/extracted_for_CNV/{sample}.extracted.bam'
    params: interval = config['RES'] + config['main_chr_bed']
    benchmark: os.path.join(config['BENCH'],'{sample}_newBaams.txt')
    conda: 'envs/preprocess.yaml'
    threads: config['neewBams']['n']

    shell:
        "{samtools} view -@ {threads} -L {params.interval} -b -o {output} {input.bam_files}"

def get_capture_kit_antitarget(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    antitarget = capture_kit + '.antitarget.bed'
    return antitarget

def get_capture_kit_bed(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    capture_kit_path = os.path.join(config['RES'], config['kit_folder'], capture_kit + '_hg38.bed')
    return capture_kit_path

rule autobin:
    input:
        bam = expand("{bam}/extracted_for_CNV/{sample_name}.extracted.bam", bam = config['BAM'], sample_name = sample_names),
    output:
        target = config['CNVKIT'] + '/target/{sample}_covered.bed',
    params:
            inputs = expand(" {bam}/extracted_for_CNV/{sample_name}.extracted.bam", bam = config['BAM'], sample_name = sample_names),
            iBED = get_capture_kit_bed,
            antitarget = get_capture_kit_antitarget,
            anno = config['RES'] + config['refflat']
    benchmark: os.path.join(config['BENCH'], '{sample}_autobin.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py autobin {params.inputs} --annotate {params.anno} --target-output-bed {output.target} -f {ref} -t {params.iBED} --antitarget-output-bed {params.antitarget}"
#
rule coverage_target:
    input:
        bam = rules.neewBams.output.NBams,
    output:
        target_cov = config['CNVKIT'] + '/coverage/{sample}.targetcoverage.cnn'
    params:
        iBED = get_capture_kit_bed,
    benchmark: os.path.join(config['BENCH'],'{sample}_coverage_target.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py coverage {input.bam} {params.iBED} -o {output.target_cov}"

rule coverage_antitarget:
    input:
        bam = rules.neewBams.output.NBams,
        ord = expand('{cnvkit}/target/{sample}_covered.bed', sample=sample_names, cnvkit = config['CNVKIT'])
    output:
        antitarget_cov = config['CNVKIT'] + '/coverage/{sample}.antitargetcoverage.cnn'
    params:
        anti_BED = get_capture_kit_antitarget
    benchmark: os.path.join(config['BENCH'],'{sample}_coverage_antitarget.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py coverage {input.bam} {params.anti_BED} -o {output.antitarget_cov}"

rule reference:
    input:
        expand('{cnvkit}/coverage/{sample}.targetcoverage.cnn', sample = sample_names, cnvkit = config['CNVKIT']),
        expand('{cnvkit}/coverage/{sample}.antitargetcoverage.cnn', sample=sample_names, cnvkit = config['CNVKIT'])
    output:
        cnv_ref = config['CNVKIT'] + '/Reference.cnn'
    params: cnvkit = config['CNVKIT']
    benchmark: os.path.join(config['BENCH'],'ReferenceCNN.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py reference {params.cnvkit}/coverage/*coverage.cnn -f {ref} -o {output}"

rule fix:
    input:
        target_cov = rules.coverage_target.output.target_cov,
        antitarget_cov = rules.coverage_antitarget.output.antitarget_cov,
        Ref = rules.reference.output.cnv_ref
    output:
        corrected = config['CNVKIT'] + '/fix/{sample}.cnr'
    benchmark: os.path.join(config['BENCH'],'{sample}_fixcnn.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py fix {input.target_cov} {input.antitarget_cov} {input.Ref} -o {output}"

rule segment:
    input:
        rules.fix.output.corrected
    output:
        segmeted = config['CNVKIT'] + '/seg/{sample}.cns'
    benchmark: os.path.join(config['BENCH'],'{sample}_seg.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py segment {input} -o {output}"

rule plot_scatter:
    input:
        cnr = rules.fix.output.corrected,
        cns = rules.segment.output.segmeted
    output:
        scatter = config['CNVKIT'] + '/plots/scatter/{sample}-scatter.pdf'
    benchmark: os.path.join(config['BENCH'],'{sample}_plotscatter.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py scatter {input.cnr} -s {input.cns} -o {output}"

rule plot_diagram:
    input:
        cnr = rules.fix.output.corrected,
        cns = rules.segment.output.segmeted
    output:
        scatter = config['CNVKIT'] + '/plots/diagram/{sample}-diagram.pdf'
    benchmark: os.path.join(config['BENCH'],'{sample}_plotdiagram.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py diagram {input.cnr} -s {input.cns} -o {output}"

rule call:
    input:
        rules.segment.output.segmeted
    output:
        call = config['CNVKIT'] + '/call/{sample}.cal.cns'
    benchmark: os.path.join(config['BENCH'],'{sample}_call_cnv.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py call {input} -o {output}"


rule define_sex:
    input:
        ref = config['CNVKIT'] + '/Reference.cnn',
        target = expand('{cnvkit}/coverage/{sample}.targetcoverage.cnn',sample=sample_names, cnvkit = config['CNVKIT']),
        antitarget = expand('{cnvkit}/coverage/{sample}.antitargetcoverage.cnn',sample=sample_names, cnvkit = config['CNVKIT']),
    output:
        config['CNVKIT'] + '/stats/Sex_sample_table.tsv'
    benchmark: os.path.join(config['BENCH'],'sexdefine.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py sex {input.ref} {input.target} {input.antitarget} -o {output}"

rule annotate_cns:
    input:
        rules.call.output.call
    output:
        a_cns = 'cnvkit/annotate/{sample}_annotate_call.cns'
    params:
        anno = config['RES'] + config['refflat']
    benchmark: os.path.join(config['BENCH'],'{sample}_annotatecns.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnv_annotate.py {params.anno} {input} -o {output}"

rule annotate_cnr:
    input:
        rules.fix.output.corrected
    output:
        a_cnr = config['CNVKIT'] + '/annotate/{sample}_annotate_fix.cnr'
    params:
        anno = config['RES'] + config['refflat']
    benchmark: os.path.join(config['BENCH'],'{sample}_annotatecnr.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnv_annotate.py {params.anno} {input} -o {output}"


rule genemetrics:
    input:
        cns = rules.annotate_cns.output.a_cns,
        cnr = rules.annotate_cnr.output.a_cnr
    output:
        config['CNVKIT'] + '/stats/genemetric/{sample}_genemetric.stat'
    benchmark: os.path.join(config['BENCH'],'{sample}_genemetrics.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py genemetrics {input.cnr} -s {input.cns} -o {output}"

rule metrics:
    input:
        cnr = expand('{cnvkit}/fix/{sample}.cnr',sample=sample_names, cnvkit = config['CNVKIT']),
        cns = expand('{cnvkit}/call/{sample}.cal.cns', sample = sample_names, cnvkit = config['CNVKIT']),
    output:
        config['CNVKIT'] + '/stats/Metrics_table.tsv'
    benchmark: os.path.join(config['BENCH'],'cnvmetrics.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        "cnvkit.py metrics {input.cnr} -s {input.cns} -o {output}"

rule breaks:
    input:
        cns = rules.annotate_cns.output.a_cns,
        cnr = rules.annotate_cnr.output.a_cnr
    output:
        config['CNVKIT'] + '/breaks/{sample}_breaks'
    conda: 'envs/preprocess.yaml'
    benchmark: os.path.join(config['BENCH'],'{sample}_breaks.txt')
    shell:
        "cnvkit.py breaks {input.cnr} {input.cns} -o {output}"

rule export_as_vcf:
    input:
        cns = rules.annotate_cns.output.a_cns,
    output:
        config['CNVKIT'] + '/VCF/{sample}_cnv.vcf'
    benchmark: os.path.join(config['BENCH'],'{sample}_export_as_vcf.txt')
    conda: 'envs/preprocess.yaml'
    shell:
        """cnvkit.py export vcf {input.cns} -i "{wildcards.sample}" -o {output}"""



#         PON creation:
#     cov2bed create few files with different coverage, lets say 1, 10, 15, 20, 30, 50
#     intersect bed files across 
