import read_samples
from common import *
import utils
onsuccess: shell("rm -fr logs/*")
wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    extension=r'sam|bam|cram',
    filetype = r'fq|fastq',
    batchnr=r'[\d]+',
    readid=r'R1|R2'
    # readgroup="[\w\d_\-@]+"

rule kraken_all:
    input:
        expand("{kraken}/{sample}.report.tsv",sample=sample_names,kraken =KRAKEN)
    default_target: True


rule kraken:
    input:
        fastq1=pj(FQ_BADMAP,"{sample}.badmap.R1.fastq.gz"),
        fastq2=pj(FQ_BADMAP,"{sample}.badmap.R2.fastq.gz")
    output:
        report=pj(KRAKEN, "{sample}.report.tsv"),
        read_cls=pj(KRAKEN, "{sample}.read_classification.tsv")
    conda: CONDA_KRAKEN
    resources:
        n="1",
        mem_mb=60000
    shell: """
        kraken2 --db {KRAKEN_DB} --threads 8 --report-minimizer-data --report {output.report} --output {output.read_cls} --paired {input.fastq1} {input.fastq2}
    """


rule bracken:
    input:
        report = rules.kraken.output.report
        # report=pj(KRAKEN, "{sample}.report.tsv")
    output:
        report=pj(KRAKEN, "{sample}.bracken_report.tsv")
    resources:
        n="1",
        mem_mb="100"
    conda: CONDA_KRAKEN
    shell:
        """
            bracken -d {KRAKEN_DB} -i {input.report} -o {output.report}
        """

rule viralmap:
    input:
        fastq1=pj(FQ_BADMAP,"{sample}.badmap.R1.fastq.gz"),
        fastq2=pj(FQ_BADMAP,"{sample}.badmap.R2.fastq.gz")
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
