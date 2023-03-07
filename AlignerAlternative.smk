import pandas as pd
import read_stats
import os
configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")

gatk = os.path.join(config['miniconda'], config['gatk'])
samtools = os.path.join(config['miniconda'], config['samtools'])
bcftools = os.path.join(config['miniconda'], config['bcftools'])
dragmap = os.path.join(config['miniconda'], config['dragmap'])
cutadapt = os.path.join(config['miniconda'], config['cutadapt'])
verifybamid2 = os.path.join(config['miniconda'], config['verifybamid2'])
ref = os.path.join(config['RES'], config['ref'])

wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+"


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

rule Aligner_all:
    input:
        expand("{bams}/{sample}.merged.bam", sample=sample_names, bams=config['BAM']),
        expand("{cram}/{sample}_mapped_hg38.cram", cram = config['CRAM'], sample=sample_names),
    default_target: True



#just alignment and convert to bams

def get_fastqpaired(wildcards):
    # command to extract path to fastq files from samplefile
    sinfo = SAMPLEINFO[wildcards['sample']]  # SMAPLEINFO comes from common.py, it's dict created from samplefile
    # check readgroups
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    file1 = os.path.join(readgroup['prefix'], readgroup['file1'])
    if file1.endswith('.bz2'):
        file1 = file1[:-4] + '.gz'
    file2 = os.path.join(readgroup['prefix'], readgroup['file2'])
    if file2.endswith('.bz2'):
        file2 = file2[:-4] + '.gz'
    return [file1, file2]

# cut adapters from inout
# DRAGMAP doesn;t work well with uBAM, so use fq as input
rule cutadapter:
    input:
        get_fastqpaired
    output:
        forr_f=os.path.join(config['FQ'], "{sample}._{readgroup}.cut_1.fq.gz"),
        rev_f=os.path.join(config['FQ'], "{sample}._{readgroup}.cut_2.fq.gz")
    # log file in this case contain some stats about removed seqs
    log:
        cutadapt_log= os.path.join(config['STAT'], "{sample}._{readgroup}.cutadapt.log"),
    benchmark:
        os.path.join(config['BENCH'], "{sample}._{readgroup}.cutadapt.txt")
    priority: 10
    conda: "preprocess"
    threads: config["cutadapter"]["n"]
    shell:
        "{cutadapt} -j {threads} -m 100 -a AGATCGGAAGAG -A AGATCGGAAGAG -o {output.forr_f} -p {output.rev_f} {input[0]} {input[1]} &> {log.cutadapt_log}"
    # run:
    #     sinfo = SAMPLEINFO[wildcards['sample']]
    #     rgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]['info']
    #     rgroupid = rgroup['ID']
    #     rgrouplib = rgroup.get('LB','unknown')
    #     rgroupplat = rgroup.get('PL','unknown')
    #     rgrouppu = rgroup.get('PU','unknown')
    #     rgroupsc = rgroup.get('CN','unknown')
    #     rgrouprd = rgroup.get('DT','unknown')
    #     # cut standard Illumina adapters
    #     cmd="""
    #     {cutadapt} -j {threads} -m 100 -a AGATCGGAAGAG -A AGATCGGAAGAG -o {output.forr_f} -p {output.rev_f} {input[0]} {input[1]} &> {log.cutadapt_log}
    #     """
    #     shell(cmd)

# rule to align reads from cutted fq on hg38 ref
# use dragmap aligner
# samtools fixmate for future step with samtools mark duplicates

def get_mem_mb_align_reads(wildcrads, attempt):
    return attempt*int(config['align_reads']['mem'])

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
        ref_dir = os.path.join(config['RES'], config['ref_dir']),
        # mask bed for current reference genome
        mask_bed = os.path.join(config['RES'], config['mask_bed'])
    conda: "preprocess"
    threads: config["align_reads"]["n"]
    log:
        dragmap_log=os.path.join(config['LOG'], "{sample}._{readgroup}_dragmap.log"),
        samtools_fixmate=os.path.join(config['LOG'], "{sample}._{readgroup}_samtools_fixamte.log"),
        samtools_sort=os.path.join(config['LOG'], "{sample}._{readgroup}_samtools_sort.log"),
        samtools_markdup=os.path.join(config['LOG'], "{sample}._{readgroup}_samtools_markdup.log"),
        samtools_index = os.path.join(config['LOG'], "{sample}._{readgroup}_samtools_index.log")
    benchmark:
        os.path.join(config['BENCH'], "{sample}._{readgroup}.dragmap.txt")
    priority: 15
    resources:
        mem_mb = get_mem_mb_align_reads
    shell:
        # "{dragmap} -r {params.ref_dir} -b {input.ubam} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --ht-mask-bed {params.mask_bed} --num-threads {threads} 2> {log.dragmap_log} |"
        "{dragmap} -r {params.ref_dir} -1 {input.for_r} -2 {input.rev_r} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --ht-mask-bed {params.mask_bed} --num-threads {threads} 2> {log.dragmap_log} | " 
        "{samtools} fixmate -@ {threads} -m - -  2> {log.samtools_fixmate} | "
        "{samtools} sort -T 'sort_temporary' -@ {threads}  -o {output.bam} 2> {log.samtools_sort} && "
        "{samtools} index -@ {threads} {output.bam} 2> {log.samtools_index}"

# # function to get information about reaadgroups
# # needed if sample contain more than 1 fastq files
def get_readgroups_bam(wildcards):
    readgroups_b = SAMPLEINFO[wildcards['sample']]['readgroups']
    files = []
    for readgroup in readgroups_b:
        files.append(os.path.join(config['BAM'],  wildcards['sample'] + '._' + readgroup['info']['ID'] + '.bam'))
    return files

# merge different readgroups bam files for same sample
rule merge_rgs:
    input:
        get_readgroups_bam
    output:
        mer_bam = os.path.join(config['BAM'],  "{sample}.merged.bam")
    log: os.path.join(config['LOG'], "{sample}.mergereadgroups.log")
    benchmark: "benchmark/{sample}.merge_rgs.txt"
    threads: config['merge_rgs']['n']
    conda: "preprocess"
    shell: "{samtools} merge -@ {threads} -o {output} {input} 2> {log}"

    # run:
    #     inputs = ' '.join(f for f in input if f.endswith('.bam'))
    #     shell("{samtools} merge -@ {threads} -o {output} {inputs} 2> {log}")

rule markdup:
    input:
        rules.merge_rgs.output.mer_bam
    output:
        mdbams = os.path.join(config['BAM'], "{sample}.markdup.bam"),
        MD_stat = os.path.join(config['STAT'], "{sample}.markdup.stat")
    benchmark: "benchmark/{sample}.markdup.txt"
    params:
        machine = 2500 #change to function
    # 100 for HiSeq and 2500 for NovaSeq
    # if fastq header not in known machines?
    # if not Illumina?
    log:
        samtools_markdup = os.path.join(config['LOG'], "{sample}.markdup.log"),
        samtools_index_md = os.path.join(config['LOG'], "{sample}.markdup_index.log")
    threads: config['markdup']['n']
    conda: "preprocess"
    shell:
        "{samtools} markdup -f {output.MD_stat} -S -d {params.machine} -@ {threads} {input} {output.mdbams} 2> {log.samtools_markdup} && "
        "{samtools} index -@ {threads} {output.mdbams} 2> {log.samtools_index_md}"


# checkpoint cause we ned to check supplemetary ratio
# if supp_ratio is too high run additional clean process
checkpoint bamstats_all:
    input:
        rules.markdup.output.mdbams
    output:
        All_stats = os.path.join(config['STAT'],  '{sample}.bam_all.tsv')
    threads: config['bamstats_all']['n']
    params: py_stats = config['BAMSTATS']
    conda: "preprocess"
    shell:
        "{samtools} view -s 0.05 -h {input} --threads {threads} | python3 {params.py_stats} stats > {output}"


# this rule triggers in case of high supp_ratio
# resort bam file before additional cleanup
rule resort_by_readname:
    input:
        rules.markdup.output.mdbams
    output: resort_bams = temp(os.path.join(config['BAM'], '{sample}_resort.bam'))
    threads: config['resort_by_readname']['n']
    conda: "preprocess"
    shell: "{samtools} sort -n -@ {threads} -o  {output} {input}"
# additional cleanup with custom script
rule declip:
    input:
        rules.resort_by_readname.output.resort_bams
    output: declip_bam = temp(os.path.join(config['BAM'], '{sample}_declip.bam'))
    threads: config['declip']['n']
    params: declip = config['DECLIP']
    conda: "preprocess"
    shell:
        "{samtools} view -s 0.05 -h {input} --threads {threads} | python3 {params.declip} > {output}"
# back to original sort order after cleanup
rule sort_back:
    input:
        rules.declip.output.declip_bam,
    output:
        ready_bams = os.path.join(config['BAM'], '{sample}.DeClipped.bam'),
        All_stats= os.path.join(config['STAT'],  '{sample}.bam_all.additional_cleanup.tsv')
    threads: config['sort_back']['n']
    params:
        py_stats= config['BAMSTATS']
    conda: "preprocess"
    shell:
        "{samtools} sort -@ {threads} -o {output.ready_bams} {input} &&"
        "{samtools} index -@ {threads} {output.ready_bams} &&"
        "{samtools} view -s 0.05 -h {input} --threads {threads} | python3 {params.py_stats} stats > {output.All_stats}"

# mapped cram
rule mCRAM:
    input:
        rules.markdup.output.mdbams
    output:
        CRAM = os.path.join(config['CRAM'], "{sample}_mapped_hg38.cram")
    threads: config['mCRAM']['n']
    benchmark: os.path.join(config['BENCH'], '{sample}_mCRAM.txt')
    conda: "preprocess"
    shell:
        "{samtools} view --cram -T {ref} -@ {threads} -o {output.CRAM} {input}"
