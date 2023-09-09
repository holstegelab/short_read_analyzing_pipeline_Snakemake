import read_samples
from common import *
import utils

wildcard_constraints:
    sample="[\w\d_\-@]+",
    extension='sam|bam|cram',
    filetype = 'fq|fastq',
    batchnr='[\d]+',
    readid='R1|R2'
    # readgroup="[\w\d_\-@]+"

rule kraken_all:
    input:
        expand("{kraken}/{sample}.report.tsv",sample=sample_names,kraken =KRAKEN),
        rules.Aligner_all.input
    default_target: True


rule kraken:
    input:
        fastq1=pj(FQ,"{sample}.badmap.R1.fastq.gz"),
        fastq2=pj(FQ,"{sample}.badmap.R2.fastq.gz")
    output:
        report=pj(KRAKEN, "{sample}.report.tsv"),
        read_cls=pj(KRAKEN, "{sample}.read_classification.tsv")
    conda: CONDA_KRAKEN
    resources:
        n="8",
        mem_mb=80000
    shell: """
        kraken2 --db {KRAKEN_DB} --threads {resources.n} --report {output.report} --output {output.read_cls} --paired {input.fastq1} {input.fastq2}
    """

rule viralmap:
    input:
        fastq1=pj(FQ,"{sample}.badmap.R1.fastq.gz"),
        fastq2=pj(FQ,"{sample}.badmap.R2.fastq.gz")
    output:
        bam=pj(KRAKEN, "{sample}.viralmap.bam"),
    params:
        dragmap=pj(SOFTWARE,dragmap),
    log:
        dragmap_log=pj(STAT, "{sample}.viralmap.dragmap.log")
    conda: CONDA_MAIN
    resources:
        n="8",
        use_threads=8,
        mem_mb=10000
    shell:
        "({params.dragmap} -r {params.viralref_dir} -1 {input.fastq1} -2 {input.fastq2} -num-threads {resources.use_threads}  | samtools view -@ 2 -o {output.bam}) 2> {log.dragmap_log} "
