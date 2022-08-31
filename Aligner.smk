import pandas as pd
import read_stats
import os
configfile: "Snakefile.cluster.json"
configfile: "Snakefile.paths.yaml"
samtools = config['miniconda'] + config['samtools']
dragmap = config['miniconda'] + config['dragmap']
cutadapt = config['miniconda'] + config['cutadapt']
ref = config['RES'] + config['ref']

wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+"


from read_samples import *
from common import *

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

rule all:
    input:
        expand("{bams}/{sample}.merged.bam",sample=sample_names,bams=config['BAM']),
        expand("{cram}/{sample}_mapped_hg38.cram", cram = config['CRAM'], sample=sample_names),


#just alignment and convert to bams

def get_fastqpaired(wildcards):
    # command to extract path to fastq files from samplefile
    sinfo = SAMPLEINFO[wildcards['sample']]  # SMAPLEINFO comes from common.py, it's dict created from samplefile
    # check readgroups
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    # SAMPLEFOLDER it's path folder
    # maybe created symlinks to folder with all fq easier?
    # ln -s
    file1 = os.path.join(config['SAMPLEFOLDER'],readgroup['file1'])
    if file1.endswith('.bz2'):
        file1 = file1[:-4] + '.gz'
    file2 = os.path.join(config['SAMPLEFOLDER'],readgroup['file2'])
    if file2.endswith('.bz2'):
        file2 = file2[:-4] + '.gz'
    return [file1, file2]

# cut adapters from inout
# DRAGMAP doesn;t work well with uBAM, so use fq as input
rule cutadapter:
    input:
        get_fastqpaired
    output:
        forr_f=config['FQ'] + "/{sample}._{readgroup}.cut_1.fq.gz",
        rev_f=config['FQ'] + "/{sample}._{readgroup}.cut_2.fq.gz"
    # log file in this case contain some stats about removed seqs
    log:
        cutadapt_log= config['STAT'] + "/{sample}._{readgroup}.cutadapt.log",
    benchmark:
        config['BENCH'] + "/{sample}._{readgroup}.cutadapt.txt"
    priority: 10
    threads: config["cutadapter"]["n"]
    run:
        sinfo = SAMPLEINFO[wildcards['sample']]
        rgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]['info']
        rgroupid = rgroup['ID']
        rgrouplib = rgroup.get('LB','unknown')
        rgroupplat = rgroup.get('PL','unknown')
        rgrouppu = rgroup.get('PU','unknown')
        rgroupsc = rgroup.get('CN','unknown')
        rgrouprd = rgroup.get('DT','unknown')
        # cut standard Illumina adapters
        cmd="""
        {cutadapt} -j {threads} -m 100 -a AGATCGGAAGAG -A AGATCGGAAGAG -o {output.forr_f} -p {output.rev_f} {input[0]} {input[1]} &> {log.cutadapt_log}
        """
        shell(cmd)

# rule to align reads from cutted fq on hg38 ref
# use dragmap aligner
# samtools fixmate for future step with samtools mark duplicates

def get_mem_mb_align_reads(wildcrads, attempt):
    return attempt*1.5*int(config['align_reads']['mem'])

rule align_reads:
    input:
        # ubam = rules.sort_ubam.output.sortubam
        for_r = rules.cutadapter.output.forr_f,
        rev_r = rules.cutadapter.output.rev_f
        # for_r = rules.uBAM_to_fq.output.f1_masked,
        # rev_r = rules.uBAM_to_fq.output.f2_masked
    output:
        bam=config['BAM'] + "/{sample}._{readgroup}.bam"
    params:
        ref_dir = config['RES'] + config['ref_dir'],
        # mask bed for current reference genome
        mask_bed = config['RES'] + config['mask_bed']
    threads: config["align_reads"]["n"]
    log:
        dragmap_log=config['LOG'] + '/' + "{sample}._{readgroup}_dragmap.log",
        samtools_fixmate=config['LOG'] + '/' + "{sample}._{readgroup}_samtools_fixamte.log",
        samtools_sort=config['LOG'] + '/' + "{sample}._{readgroup}_samtools_sort.log",
        samtools_markdup=config['LOG'] + '/' + "{sample}._{readgroup}_samtools_markdup.log",
        samtools_index = config['LOG'] + '/' + "{sample}._{readgroup}_samtools_index.log"
    benchmark:
        config['BENCH'] + "/{sample}._{readgroup}.dragmap.txt"
    priority: 15
    resources:
        mem_mb = get_mem_mb_align_reads
    shell:
        # "{dragmap} -r {params.ref_dir} -b {input.ubam} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --ht-mask-bed {params.mask_bed} --num-threads {threads} 2> {log.dragmap_log} |"
        "{dragmap} -r {params.ref_dir} -1 {input.for_r} -2 {input.rev_r} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --ht-mask-bed {params.mask_bed} --num-threads {threads} 2> {log.dragmap_log} | " 
        "{samtools} fixmate -@ {threads} -m - -  2> {log.samtools_fixmate} | "
        "{samtools} sort -T 'sort_temporary' -@ {threads}  -o {output.bam} 2> {log.samtools_sort} && "
        "{samtools} index -@ {threads} {output.bam} 2> {log.samtools_index}"

# function to get information about reaadgroups
# needed if sample contain more than 1 fastq files
def get_readgroups(wildcards):
    readgroups = SAMPLEINFO[wildcards['sample']]['readgroups']
    files = []
    for readgroup in readgroups:
        files.append(os.path.join(config['BAM'] + '/' + wildcards['sample'] + '._' + readgroup['info']['ID'] + '.bam'))
    return files

# merge different readgroups bam files for same sample
rule merge_rgs:
    input:
        get_readgroups
    output:
        mer_bam = (config['BAM'] + "/{sample}.merged.bam")
    log: config['LOG'] + '/' + "{sample}.mergereadgroups.log"
    benchmark: "benchmark/{sample}.merge_rgs.txt"
    threads: config['merge_rgs']['n']
    run:
        inputs = ' '.join(f for f in input if f.endswith('.bam'))
        shell("{samtools} merge -@ {threads} -o {output} {inputs} 2> {log}")

rule markdup:
    input:
        rules.merge_rgs.output.mer_bam
    output:
        mdbams = config['BAM'] + "/{sample}.markdup.bam",
        MD_stat = config['STAT'] + "/{sample}.markdup.stat"
    benchmark: "benchmark/{sample}.markdup.txt"
    params:
        machine = 2500 #change to function
    # 100 for HiSeq and 2500 for NovaSeq
    # if fastq header not in known machines?
    # if not Illumina?
    log:
        samtools_markdup = config['LOG'] + '/' + "{sample}.markdup.log",
        samtools_index_md = config['LOG'] + '/' + "{sample}.markdup_index.log"
    threads: config['markdup']['n']
    shell:
        "{samtools} markdup -f {output.MD_stat} -S -d {params.machine} -@ {threads} {input} {output.mdbams} 2> {log.samtools_markdup} && "
        "{samtools} index -@ {threads} {output.mdbams} 2> {log.samtools_index_md}"


# checkpoint cause we ned to check supplemetary ratio
# if supp_ratio is too high run additional clean process
checkpoint bamstats_all:
    input:
        rules.markdup.output.mdbams
    output:
        All_stats = config['STAT'] + '/{sample}.bam_all.tsv'
    threads: config['bamstats_all']['n']
    params: py_stats = config['BAMSTATS']
    shell:
        "{samtools} view -s 0.05 -h {input} --threads {threads} | python3 {params.py_stats} stats > {output}"


# this rule triggers in case of high supp_ratio
# resort bam file before additional cleanup
rule resort_by_readname:
    input:
        rules.markdup.output.mdbams
    output: resort_bams = temp(config['BAM'] + '/{sample}_resort.bam')
    threads: config['resort_by_readname']['n']
    shell: "{samtools} sort -n -@ {threads} -o  {output} {input}"
# additional cleanup with custom script
rule declip:
    input:
        rules.resort_by_readname.output.resort_bams
    output: declip_bam = temp(config['BAM'] + '/{sample}_declip.bam')
    threads: config['declip']['n']
    params: declip = config['DECLIP']
    shell:
        "{samtools} view -s 0.05 -h {input} --threads {threads} | python3 {params.declip} > {output}"
# back to original sort order after cleanup
rule sort_back:
    input:
        rules.declip.output.declip_bam,
    output:
        ready_bams = config['BAM'] + '/{sample}.DeClipped.bam',
        All_stats= config['STAT'] + '/{sample}.bam_all.additional_cleanup.tsv'
    threads: config['sort_back']['n']
    params:
        py_stats= config['BAMSTATS']
    shell:
        "{samtools} sort -@ {threads} -o {output.ready_bams} {input} &&"
        "{samtools} index -@ {threads} {output.ready_bams} &&"
        "{samtools} view -s 0.05 -h {input} --threads {threads} | python3 {params.py_stats} stats > {output.All_stats}"

# mapped cram
rule mCRAM:
    input:
        rules.markdup.output.mdbams
    output:
        CRAM = config['CRAM'] + "/{sample}_mapped_hg38.cram"
    threads: config['mCRAM']['n']
    benchmark: config['BENCH'] + '/{sample}_mCRAM.txt'
    shell:
        "{samtools} view --cram -T {ref} -@ {threads} -o {output.CRAM} {input}"
