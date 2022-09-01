configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")
gatk = config['miniconda'] + config['gatk']
samtools = config['miniconda'] + config['samtools']
ref = config['RES'] + config['ref']

wildcard_constraints:
    sample="[\w\d_\-@]+"
from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO = load_samplefiles('.', config)
# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()
rule all:
    input:
        'cnvkit/Reference.cnn',
        expand('cnvkit/target/{sample}_covered.bed',sample=sample_names),
        expand('cnvkit/coverage/{sample}.targetcoverage.cnn',sample=sample_names),
        expand('cnvkit/coverage/{sample}.antitargetcoverage.cnn',sample=sample_names),
        expand('cnvkit/fix/{sample}.cnr',sample=sample_names),
        expand('cnvkit/plots/diagram/{sample}-diagram.pdf', sample=sample_names),
        expand('cnvkit/plots/scatter/{sample}-scatter.pdf', sample=sample_names),
        expand('cnvkit/call/{sample}.cal.cns', sample = sample_names),
        'cnvkit/stats/Sex_sample_table.tsv',
        expand('cnvkit/stats/genemetric/{sample}_genemetric.stat', sample = sample_names),
        'cnvkit/stats/Metrics_table.tsv',
        expand('cnvkit/breaks/{sample}_breaks', sample = sample_names)


module Aligner:
    snakefile: 'Aligner.smk'

# bams that align only on "main" chromosomes
rule neewBams:
    input:
        bam_files = rules.markdup.output.mdbams
    output:
        NBams = config['BAM'] + '/extracted_for_CNV/{sample}.extracted.bam'
    params: interval = config['RES'] + config['main_chr_bed']
    threads: config['neewBams']['n']

    shell:
        "{samtools} view -@ {threads} -L {params.interval} -b -o {output} {input}"


rule autobin:
    input:
        bam = rules.neewBams.output.NBams
    output:
        target = config['TARGET'] + '/{sample}_covered.bed'
    params:
            inputs = list(map(" {}/extracted_for_CNV/{}.extracted.bam".format,config['BAM'], sample_names)),
            iBED = get_capture_kit_bed,
            anno = '/projects/0/qtholstg/hg38_res/refFlat.txt'
    shell:
        "cnvkit.py autobin {params.inputs} --annotate {params.anno} --target-output-bed {output} -f {ref} -t {params.iBED}"
#
rule coverage_target:
    input:
        bam = rules.neewBams.output.NBams,
    output:
        target_cov = 'cnvkit/coverage/{sample}.targetcoverage.cnn'
    params:
        iBED = '/projects/0/qtholstg/hg38_res/intervals/Agilent_V6_hg38.bed',
    benchmark: 'bench/{sample}_target_cov.txt'
    shell:
        "cnvkit.py coverage {input.bam} {params.iBED} -o {output.target_cov}"

rule coverage_antitarget:
    input:
        bam = rules.neewBams.output.NBams,
        ord = expand('cnvkit/target/{sample}_covered.bed', sample=sample_names)
    output:
        antitarget_cov = 'cnvkit/coverage/{sample}.antitargetcoverage.cnn'
    params:
        anti_BED = 'Agilent_V6_hg38.antitarget.bed',
    benchmark: 'bench/{sample}_antitarget_cov.txt'
    shell:
        "cnvkit.py coverage {input.bam} {params.anti_BED} -o {output.antitarget_cov}"

rule reference:
    input:
        expand('cnvkit/coverage/{sample}.targetcoverage.cnn', sample = sample_names),
        expand('cnvkit/coverage/{sample}.antitargetcoverage.cnn', sample=sample_names)
    output:
        cnv_ref = 'cnvkit/Reference.cnn'
    params:
        ref = '/projects/0/qtholstg/hg38_res/Ref/GRCh38_full_analysis_set_plus_decoy_hla.fa',
    shell:
        "cnvkit.py reference cnvkit/coverage/*coverage.cnn -f {params.ref} -o {output}"

rule fix:
    input:
        target_cov = rules.coverage_target.output.target_cov,
        antitarget_cov = rules.coverage_antitarget.output.antitarget_cov,
        Ref = rules.reference.output.cnv_ref
    output:
        corrected = 'cnvkit/fix/{sample}.cnr'
    benchmark: 'bench/{sample}_fix.txt'
    shell:
        "cnvkit.py fix {input.target_cov} {input.antitarget_cov} {input.Ref} -o {output}"

rule segment:
    input:
        rules.fix.output.corrected
    output:
        segmeted = 'cnvkit/seg/{sample}.cns'
    benchmark: 'bench/{sample}_seg.txt'
    shell:
        "cnvkit.py segment {input} -o {output}"

rule plot_scatter:
    input:
        cnr = rules.fix.output.corrected,
        cns = rules.segment.output.segmeted
    output:
        scatter = 'cnvkit/plots/scatter/{sample}-scatter.pdf'
    shell:
        "cnvkit.py scatter {input.cnr} -s {input.cns} -o {output}"

rule plot_diagram:
    input:
        cnr = rules.fix.output.corrected,
        cns = rules.segment.output.segmeted
    output:
        scatter = 'cnvkit/plots/diagram/{sample}-diagram.pdf'
    shell:
        "cnvkit.py diagram {input.cnr} -s {input.cns} -o {output}"

rule call:
    input:
        rules.segment.output.segmeted
    output:
        call = 'cnvkit/call/{sample}.cal.cns'
    benchmark:'bench/{sample}_call.txt'
    shell:
        "cnvkit.py call {input} -o {output}"


rule define_sex:
    input:
        ref = 'cnvkit/Reference.cnn',
        target = expand('cnvkit/coverage/{sample}.targetcoverage.cnn',sample=sample_names),
        antitarget = expand('cnvkit/coverage/{sample}.antitargetcoverage.cnn',sample=sample_names),
    output:
        'cnvkit/stats/Sex_sample_table.tsv'
    shell:
        "cnvkit.py sex {input.ref} {input.target} {input.antitarget} -o {output}"

rule annotate_cns:
    input:
        rules.call.output.call
    output:
        a_cns = 'cnvkit/annotate/{sample}_annotate_call.cns'
    params:
        anno = '/projects/0/qtholstg/hg38_res/refFlat.txt'
    shell:
        "cnv_annotate.py {params.anno} {input} -o {output}"

rule annotate_cnr:
    input:
        rules.fix.output.corrected
    output:
        a_cnr = 'cnvkit/annotate/{sample}_annotate_fix.cnr'
    params:
        anno = '/projects/0/qtholstg/hg38_res/refFlat.txt'
    shell:
        "cnv_annotate.py {params.anno} {input} -o {output}"


rule genemetrics:
    input:
        cns = rules.annotate_cns.output.a_cns,
        cnr = rules.annotate_cnr.output.a_cnr
    output:
        'cnvkit/stats/genemetric/{sample}_genemetric.stat'
    shell:
        "cnvkit.py genemetrics {input.cnr} -s {input.cns} -o {output}"

rule metrics:
    input:
        cnr = expand('cnvkit/fix/{sample}.cnr',sample=sample_names),
        cns = expand('cnvkit/call/{sample}.cal.cns', sample = sample_names),
    output:
        'cnvkit/stats/Metrics_table.tsv'
    shell:
        "cnvkit.py metrics {input.cnr} -s {input.cns} -o {output}"

rule breaks:
    input:
        cns = rules.annotate_cns.output.a_cns,
        cnr = rules.annotate_cnr.output.a_cnr
    output:
        'cnvkit/breaks/{sample}_breaks'
    shell:
        "cnvkit.py breaks {input.cnr} {input.cns} -o {output}"
