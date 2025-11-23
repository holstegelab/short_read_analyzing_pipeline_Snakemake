import read_samples
from common import *

import utils
from shlex import quote
onsuccess: shell("rm -fr logs/*")
wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    extension=r'sam|bam|cram',
    filetype = r'fq|fastq',
    batchnr=r'[\d]+',
    readid=r'R1|R2'
    # readgroup="[\w\d_\-@]+"

module Aligner:
    snakefile: 'Aligner.smk'
    config: config

use rule * from Aligner

rule kraken_all:
    input:
        expand("{kraken}/{sample}.kraken_summary.tsv", sample=sample_names, kraken=KRAKEN)    


def _samplefile_samples(samplefile):
    return sorted(SAMPLEFILE_TO_SAMPLES[samplefile].keys())


def kraken_report_inputs(wildcards):
    files = []
    for sample in _samplefile_samples(wildcards.samplefile):
        files.extend([
            pj(KRAKEN, f"{sample}.report.tsv"),
            pj(KRAKEN, f"{sample}.bracken_report.tsv"),
            pj(KRAKEN, f"{sample}.kraken_summary.tsv"),
            pj(KRAKEN, f"{sample}.report_bracken_species.tsv"),
        ])
    return files


def kraken_read_classification_inputs(wildcards):
    return [pj(KRAKEN, f"{sample}.read_classification.tsv.gz") for sample in _samplefile_samples(wildcards.samplefile)]


rule tar_kraken_reports:
    input:
        reports=kraken_report_inputs
    output:
        tar=temp(pj(KRAKEN, "{samplefile}.kraken_reports.tar.gz"))
    params:
        files=lambda wildcards, input: " ".join(quote(path) for path in input.reports),
        outdir=lambda wildcards, output: quote(os.path.dirname(output.tar))
    resources:
        mem_mb=2000,
        n="0.5"
    shell:
        """
        mkdir -p {params.outdir}
        tar -czf {output.tar:q} {params.files}
        """


rule tar_kraken_read_classification:
    input:
        reads=kraken_read_classification_inputs
    output:
        tar=temp(pj(KRAKEN, "{samplefile}.kraken_read_classification.tar.gz"))
    params:
        files=lambda wildcards, input: " ".join(quote(path) for path in input.reads),
        outdir=lambda wildcards, output: quote(os.path.dirname(output.tar))
    resources:
        mem_mb=2000,
        n="0.5"
    shell:
        """
        mkdir -p {params.outdir}
        tar -czf {output.tar:q} {params.files}
        """


rule copy_kraken_reports_to_dcache:
    input:
        tar=pj(KRAKEN, "{samplefile}.kraken_reports.tar.gz")
    output:
        copied=temp(pj(KRAKEN, "{samplefile}.kraken_reports.tar.copied")),
        checksum=temp(pj(KRAKEN, "{samplefile}.kraken_reports.tar.ADLER32"))
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1",
        dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        remote_dir = os.path.join(remote_base_for_samplefile(wildcards.samplefile), "kraken")
        remote_name = os.path.basename(input.tar)
        copy_with_checksum(str(input.tar), remote_dir, remote_name, str(output.checksum), AGH_DCACHE_CONFIG, params.ada_script)
        shell(f"touch {quote(str(output.copied))}")


rule copy_kraken_read_classification_to_dcache:
    input:
        tar=pj(KRAKEN, "{samplefile}.kraken_read_classification.tar.gz")
    output:
        copied=temp(pj(KRAKEN, "{samplefile}.kraken_read_classification.tar.copied")),
        checksum=temp(pj(KRAKEN, "{samplefile}.kraken_read_classification.tar.ADLER32"))
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1",
        dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        remote_dir = os.path.join(remote_base_for_samplefile(wildcards.samplefile), "kraken")
        remote_name = os.path.basename(input.tar)
        copy_with_checksum(str(input.tar), remote_dir, remote_name, str(output.checksum), AGH_DCACHE_CONFIG, params.ada_script)
        shell(f"touch {quote(str(output.copied))}")


rule kraken_tar_all:
    input:
        expand(pj(KRAKEN, "{samplefile}.kraken_reports.tar.copied"), samplefile=SAMPLE_FILES),
        expand(pj(KRAKEN, "{samplefile}.kraken_read_classification.tar.copied"), samplefile=SAMPLE_FILES)
    output:
        done=touch(pj(KRAKEN, "kraken_uploads.done"))


##THIS FUNCION IS COPIED ALSO in Aligner.smk and Stat.smk
def sampleinfo(SAMPLEINFO, sample, checkpoint=False):  #{{{
    """If samples are on tape, we do not have sample readgroup info.
    That is, the 'readgroups' field is empty.

    This function first checks if the readgroup info is available on disk,
    in the file SAMPLEINFODIR/<sample>.dat. 

    Alternatively, the function injects a checkpoint rule to load this readgroup info.
    """

    sinfo = SAMPLEINFO[sample]
    if not 'readgroups' in sinfo:
        rgpath = pj(SAMPLEINFODIR,sample + ".dat")
        if os.path.exists(rgpath):
            xsample = utils.load(rgpath)
        elif checkpoint:
            #no readgroup info yet
            filename = Aligner.checkpoints.get_readgroups.get(sample=sample).output[0]
            xsample = utils.load(filename)
        sinfo = sinfo.copy()
        sinfo['readgroups'] = xsample['readgroups']
        sinfo['alternative_names'] = sinfo.get('alternative_names',set()).union(xsample['alternative_names'])
        SAMPLEINFO[sample] = sinfo
    return sinfo  #}}}


def get_merge_stats_inputs(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards.sample, checkpoint=True)
    readgroups = sinfo.get('readgroups', [])
    return sorted(
        pj(STAT, f"{wildcards.sample}.{rg['info']['ID']}.merge_stats.tsv")
        for rg in readgroups
    )


rule kraken:
    input:
        fastq1=pj(FQ_BADMAP,"{sample}.badmap.R1.fastq.gz"),
        fastq2=pj(FQ_BADMAP,"{sample}.badmap.R2.fastq.gz")
    output:
        report=temp(pj(KRAKEN, "{sample}.report.tsv")),
        read_cls=temp(pj(KRAKEN, "{sample}.read_classification.tsv.gz"))
    conda: CONDA_KRAKEN
    resources:
        n="1.5",
        mem_mb=65000
    shell: """
        kraken2 --db {KRAKEN_DB} --threads 8 --report-minimizer-data --report {output.report} --output >(gzip -c > {output.read_cls}) --paired {input.fastq1} {input.fastq2}
    """


rule bracken:
    input:
        report = rules.kraken.output.report
        # report=pj(KRAKEN, "{sample}.report.tsv")
    output:
        report=temp(pj(KRAKEN, "{sample}.bracken_report.tsv")),
        species_report=temp(pj(KRAKEN, "{sample}.report_bracken_species.tsv"))
    resources:
        n="0.5",
        mem_mb="200"
    conda: CONDA_KRAKEN
    shell:
        """
            bracken -d {KRAKEN_DB} -i {input.report} -o {output.report}
            bracken -d {KRAKEN_DB} -i {input.report} -l S -o {output.species_report}
        """


rule kraken_summary:
    input:
        report=rules.kraken.output.report,
        bracken=rules.bracken.output.report,
        merge_stats=get_merge_stats_inputs
    output:
        summary=temp(pj(KRAKEN, "{sample}.kraken_summary.tsv"))
    conda: CONDA_KRAKEN
    resources:
        n="0.5",
        mem_mb=200
    params:
        summary_script=srcdir("scripts/kraken_summary.py"),
        merge_stats=lambda wildcards, input: ",".join(input.merge_stats)
    shell: """
        python {params.summary_script} --sample {wildcards.sample} --report {input.report} --merge-stats {params.merge_stats} --bracken {input.bracken} --output {output.summary}
    """

rule viralmap:
    input:
        fastq1=pj(FQ_BADMAP,"{sample}.badmap.R1.fastq.gz"),
        fastq2=pj(FQ_BADMAP,"{sample}.badmap.R2.fastq.gz")
    output:
        bam=pj(KRAKEN, "{sample}.viralmap.bam"),
        dragmap_log=temp(pj(STAT, "{sample}.viralmap.dragmap.log"))
    params:
        dragmap=pj(SOFTWARE,dragmap),
    conda: CONDA_MAIN
    resources:
        n="8",
        use_threads=8,
        mem_mb=10000
    shell:
        "({params.dragmap} -r {params.viralref_dir} -1 {input.fastq1} -2 {input.fastq2} -num-threads {resources.use_threads}  | samtools view -@ 2 -o {output.bam}) 2> {output.dragmap_log} "
