import pandas as pd
import read_stats
import os
import getpass
configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")

gatk = os.path.join(config['miniconda'], config['gatk'])
samtools = os.path.join(config['miniconda'], config['samtools'])
bcftools = os.path.join(config['miniconda'], config['bcftools'])
dragmap = os.path.join(config['miniconda'], config['dragmap'])
cutadapt = os.path.join(config['miniconda'], config['cutadapt'])
verifybamid2 = os.path.join(config['miniconda'], config['verifybamid2'])
ref = os.path.join(config['RES'], config['ref'])



tmpdir = os.path.join(config['TMPDIR'], getpass.getuser())


os.makedirs(tmpdir, mode=0o700, exist_ok=True)

wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+"


from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO = load_samplefiles('.', config)

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

def get_readgroup_params(wildcards):
    res = [rg for rg in SAMPLEINFO[wildcards['sample']]['readgroups'] if rg['info']['ID'] == wildcards['readgroup']][0]['info']
    
    return {'ID':res['ID'], 'LB':res.get('LB','unknown'), 'PL':res.get('PL','unknown'), 'PU':res.get('PU','unknown'), \
            'CN':res.get('CN','unknown'), 'DT':res.get('DT','unknown')}

rule create_unaligned_bam_fastq:
    input:
        get_fastqpaired
    output:
        bam=os.path.join(config['BAM'], "{sample}.{readgroup}.unaligned.bam")
    log:
        markadapters=os.path.join(config['LOG'], "{sample}.{readgroup}.markilluminaadapers.log")
    benchmark:
        os.path.join(config['BENCH'], "{sample}.{readgroup}.create_unaliged_bam_fastq.txt")
    conda: "picard"
    params: 
        dname=lambda wildcards, output: os.path.dirname(output.bam),
        lb = lambda wildcards: get_readgroup_params(wildcards).get('LB','unknown'),
        pl = lambda wildcards: get_readgroup_params(wildcards).get('PL','unknown'),
        pu = lambda wildcards: get_readgroup_params(wildcards).get('PU','unknown'),
        cn = lambda wildcards: get_readgroup_params(wildcards).get('CN','unknown'),
        dt = lambda wildcards: get_readgroup_params(wildcards).get('DT','unknown')
    threads: config['create_unaligned_bam_fastq']['n']
    resources:
        tmpdir=tmpdir
    shell:
        """
        mkdir -p {params.dname}
        picard -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx16g -Xms8g FastqToSam TMP_DIR={resources.tmpdir} FASTQ="{input[0]}" FASTQ2="{input[1]}" OUTPUT=/dev/stdout READ_GROUP_NAME={wildcards.readgroup} MAX_RECORDS_IN_RAM=10000000 SAMPLE_NAME={wildcards.sample} LIBRARY_NAME={params.lb} PLATFORM={params.pl} PLATFORM_UNIT={params.pu} SEQUENCING_CENTER={params.cn} RUN_DATE={params.dt} COMPRESSION_LEVEL=0 QUIET=true | \
        picard -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx16g -Xms8g MarkIlluminaAdapters METRICS={log.markadapters} COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=SILENT  INPUT=/dev/stdin OUTPUT={output.bam}
        """
# rule to align reads from cutted fq on hg38 ref
# use dragmap aligner
# samtools fixmate for future step with samtools mark duplicates

def get_mem_mb_align_reads(wildcrads, attempt):
    return attempt*int(config['align_reads']['mem'])

rule align_reads:
    input:
        in_bam=rules.create_unaligned_bam_fastq.output.bam
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.aligned.bam")
    params:
        ref_dir = os.path.join(config['RES'], config['ref_dir']),
        # mask bed for current reference genome
        mask_bed = os.path.join(config['RES'], config['mask_bed']),
        temp_sort = os.path.join("sort_temporary_{sample}_{readgroup}")
    conda: "preprocess"
    threads: config["align_reads"]["n"]
    log:
        dragmap_log=os.path.join(config['LOG'], "{sample}.{readgroup}.dragmap.log"),
    benchmark:
        os.path.join(config['BENCH'], "{sample}.{readgroup}.dragmap.txt")
    priority: 15
    resources:
        mem_mb = get_mem_mb_align_reads
    shell:
        "{dragmap} -r {params.ref_dir} -b {input.in_bam} --interleaved 1 --input-qname-suffix-delimiter / --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --ht-mask-bed {params.mask_bed} --num-threads {threads} 2> {log.dragmap_log} | samtools view -@ {threads} -o {output.bam}" 
        #--preserve-map-align-order 1 was tested, so that unaligned and aligned bam have sam read order (thread synchronization). But reduces performance by 1/3.  Better to let mergebamalignment deal with the issue.

rule merge_bam_alignment:
    input:
        in_bam_unaligned = rules.create_unaligned_bam_fastq.output.bam,
        in_bam_aligned = rules.align_reads.output.bam
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.merged.bam")
    params:
        ref_dir = os.path.join(config['RES'], config['ref_dir']),
        # mask bed for current reference genome
        mask_bed = os.path.join(config['RES'], config['mask_bed']),
        temp_sort = os.path.join("sort_temporary_{sample}_{readgroup}")
    conda: "picard"
    threads: config["merge_bam_alignment"]["n"]
    benchmark:
        os.path.join(config['BENCH'], "{sample}.{readgroup}.mergebam.txt")
    priority: 15
    resources:
        mem_mb = get_mem_mb_align_reads,
        tmpdir=tmpdir
    shell:
        # "{dragmap} -r {params.ref_dir} -b {input.ubam} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --ht-mask-bed {params.mask_bed} --num-threads {threads} 2> {log.dragmap_log} |"
       """
            picard -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx16g -Xms8g MergeBamAlignment \
                TMP_DIR={resources.tmpdir} \
                VALIDATION_STRINGENCY=SILENT \
                EXPECTED_ORIENTATIONS=FR \
                ATTRIBUTES_TO_RETAIN=XS \
                ATTRIBUTES_TO_REMOVE=MD \
                ALIGNED_BAM={input.in_bam_aligned} \
                UNMAPPED_BAM={input.in_bam_unaligned} \
                OUTPUT={output.bam} \
                REFERENCE_SEQUENCE={ref} \
                PAIRED_RUN=true \
                COMPRESSION_LEVEL=0 \
                SORT_ORDER="unsorted" \
                IS_BISULFITE_SEQUENCE=false \
                ALIGNED_READS_ONLY=false \
                CLIP_ADAPTERS=true \
                CLIP_OVERLAPPING_READS=true \
                HARD_CLIP_OVERLAPPING_READS=true \
                MAX_RECORDS_IN_RAM=2000000 \
                ADD_MATE_CIGAR=false \
                MAX_INSERTIONS_OR_DELETIONS=-1 \
                PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
                UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
                ALIGNER_PROPER_PAIR_FLAGS=true \
                UNMAP_CONTAMINANT_READS=true \
                ADD_PG_TAG_TO_READS=false
        """

rule sort_bam_alignment:
    input:
        in_bam = rules.merge_bam_alignment.output.bam
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.sorted.bam")
    params:
        # mask bed for current reference genome
        temp_sort = os.path.join("sort_temporary_{sample}_{readgroup}")
    conda: "preprocess"
    threads: config["sort_bam_alignment"]["n"]
    log:
        samtools_fixmate=os.path.join(config['LOG'], "{sample}.{readgroup}.samtools_fixmate.log"),
        samtools_sort=os.path.join(config['LOG'], "{sample}.{readgroup}.samtools_sort.log"),
        samtools_index = os.path.join(config['LOG'], "{sample}.{readgroup}.samtools_index.log")
    benchmark:
        os.path.join(config['BENCH'], "{sample}.{readgroup}.mergebam.txt")
    priority: 15
    resources:
        mem_mb = get_mem_mb_align_reads,
        tmpdir=tmpdir
    shell:
        """
            {samtools} fixmate -@ {threads} -u -m {input.in_bam} -  2> {log.samtools_fixmate} |\
            {samtools} sort -T {resources.tmpdir}/{params.temp_sort} -@ {threads} -l 1 -m 2000M -o {output.bam} 2> {log.samtools_sort} && \
            {samtools} index -@ {threads} {output.bam} 2> {log.samtools_index}
        """



# # function to get information about reaadgroups
# # needed if sample contain more than 1 fastq files
def get_readgroups_bam(wildcards):
    readgroups_b = SAMPLEINFO[wildcards['sample']]['readgroups']
    files = []
    for readgroup in readgroups_b:
        files.append(os.path.join(config['BAM'],  wildcards['sample'] + '.' + readgroup['info']['ID'] + '.sorted.bam'))
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
    run: 
        if len(input) > 1:
            cmd = "{samtools} merge -@ {threads} {output} {input} 2> {log}"
            shell(cmd,conda_env='preprocess')
        else:
            cmd = "ln {input} {output}"
            shell(cmd)

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
    params: temp_sort = "resort_temporary_{sample}"
    threads: config['resort_by_readname']['n']
    conda: "preprocess"
    resources:
        tmpdir=tmpdir
    shell: "{samtools} sort -T {resources.tmpdir}/{params.temp_sort} -n -@ {threads} -o  {output} {input}"
# additional cleanup with custom script
rule declip:
    input:
        rules.resort_by_readname.output.resort_bams
    output: declip_bam = temp(os.path.join(config['BAM'], '{sample}_declip.bam'))
    threads: config['declip']['n']
    params: declip = srcdir(config['DECLIP'])
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
        py_stats= config['BAMSTATS'],
        temp_sort = "resort_back_temporary_{sample}"
    conda: "preprocess"
    resources:
        tmpdir=tmpdir
    shell:
        "{samtools} sort -T {resources.tmpdir}/{params.temp_sort} -@ {threads} -o {output.ready_bams} {input} &&"
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
