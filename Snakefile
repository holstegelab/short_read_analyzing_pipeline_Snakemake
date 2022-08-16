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


# main chromosomes from GRCh38 splitted into 99 bins
bins = config['RES'] + config['bin_file_ref']
import csv
chrs = []
with open(bins) as file:
    tsv_file = csv.reader(file, delimiter="\t")
    for line in tsv_file:
        out = line[0] + ':' + line[1] + '-' + line[2]
        chrs.append(out)
# print(chrs)
# chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

from read_samples import *
from common import *

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()


def get_fastqpaired(wildcards):
    # command to extract path to fastq files from samplefile
    sinfo = SAMPLEINFO[wildcards['sample']] # SMAPLEINFO comes from common.py, it's dict created from samplefile
    # check readgroups
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    # SAMPLEFOLDER it's path folder
    # maybe created symlinks to folder with all fq easier?
    # ln -s
    file1 = os.path.join(config['SAMPLEFOLDER'], readgroup['file1'])
    if file1.endswith('.bz2'):
        file1 = file1[:-4] + '.gz'
    file2 = os.path.join(config['SAMPLEFOLDER'], readgroup['file2'])
    if file2.endswith('.bz2'):
        file2 = file2[:-4] + '.gz'
    return [file1,file2]


# function to get information about reaadgroups
# needed if sample contain more than 1 fastq files
def get_readgroups(wildcards):
    readgroups = SAMPLEINFO[wildcards['sample']]['readgroups']
    files = []
    for readgroup in readgroups:
        files.append(os.path.join(config['BAM'] + '/' + wildcards['sample'] + '.' + readgroup['info']['ID'] + '.bam'))
    return files

# extract info about capture kit from SAMPLEFILE
# assume that all kits bed and interval_list files are existing and download to res folder
def get_capture_kit_bed(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    # if capture_kit.strip() == '':
    #     capture_kit = os.path.basename(MERGED_CAPTURE_KIT)[:-4]
    capture_kit_path = config['RES'] + config['kit_folder'] + capture_kit + '_hg38.bed'
    return capture_kit_path
# return interval_list file instead of bed file
def get_capture_kit_interval_list(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    capture_kit_path = config['RES'] + config['kit_folder'] + capture_kit + '_hg38.interval_list'
    return capture_kit_path

rule all:
    input:
        expand('{stat}/contam/{sample}_verifybamid.selfSM', stat = config['STAT'] , sample = sample_names),
        expand('gvcfs/{sample}.g.vcf.gz', sample = sample_names),
        expand("{bams}/{sample}-dragstr.txt", bams = config['BAM'], sample = sample_names),
        expand("{bams}/{sample}.merged.bam", sample = sample_names, bams = config['BAM']),
        # expand(config['STAT'] + "/{sample}.{readgroup}_adapters.stat", sample = sample_names)
        config['STAT'] + "/BASIC.variant_calling_detail_metrics",
        # # sample_names from sample file accord to all samples
        # # expanded version accord to all samples listed in samplefile
        expand("{stats}/{sample}_hs_metrics",sample=sample_names, stats = config['STAT']),
        # expand("{stats}/{sample}.bait_bias.bait_bias_detail_metrics", sample=sample_names, stats = config['STAT']),
        expand("{stats}/{sample}.OXOG", sample=sample_names, stats = config['STAT']),
        expand("{stats}/{sample}_samtools.stat", sample=sample_names, stats = config['STAT']),
        expand("{vcf}/Merged_raw_DBI_{chrs}.vcf.gz", chrs = chrs, vcf = config['VCF']),
        expand("{vcf}/ALL_chrs.vcf.gz", vcf = config['VCF']),
        expand("{stats}/contam/{sample}_verifybamid.selfSM", sample=sample_names, stats = config['STAT']),
        expand("{stats}/{sample}.bam_all.tsv", sample=sample_names, stats = config['STAT']),
        expand("{samplefile}.oxo_quality.tab", samplefile = SAMPLE_FILES),
        expand("{samplefile}.bam_quality.v4.tab", samplefile = SAMPLE_FILES),
        expand("{cram}/{sample}_unmapped_masked.cram", cram = config['CRAM'], sample=sample_names)

#just alignment and convert to bams


# remove known illumina adapters
# convert to uBAM for masking adapters
rule convert_to_uBAM:
    input:
        get_fastqpaired
    output:
        uBAM = temp(config['uBAM'] + "/{sample}.{readgroup}.unmapped.bam"),
        # maskedBAM = config['uBAM'] + "/{sample}.{readgroup}.unmapped.masked.bam",
        # adapters_stats= config['STAT'] + "/{sample}.{readgroup}_adapters.stat"
    benchmark: config['BENCH'] + "/{sample}.{readgroup}.unmappedBAM.txt"
    log:
        fqtosam = config['LOG'] + '/' + "{sample}.{readgroup}.fq2sam.log",
        # adapters = config['LOG'] + '/' + "{sample}.{readgroup}.adapters.log"
    shell:
        """
        {gatk} FastqToSam --FASTQ {input[0]} --FASTQ2 {input[1]} -O {output.uBAM} -SM {wildcards.sample} -RG {wildcards.readgroup}
        """
# mask adapters in uBAM file
rule mask_adapters:
    input:
        rules.convert_to_uBAM.output.uBAM
    output:
        maskedBAM = config['uBAM'] + "/{sample}.{readgroup}_unmapped_masked.bam",
        adapters_stats= config['STAT'] + "/{sample}.{readgroup}_adapters.stat"
    benchmark: config['BENCH'] + "/{sample}.{readgroup}.maskadapters.txt"
    log:
        adapters = config['LOG'] + '/' + "{sample}.{readgroup}.adapters.log"
    shell:
        """
        {gatk} MarkIlluminaAdapters -I {input} -O {output.maskedBAM} -M {output.adapters_stats} 2> {log.adapters}
        """

def get_readgroups_unmapped(wildcards):
    readgroups = SAMPLEINFO[wildcards['sample']]['readgroups']
    files = []
    for readgroup in readgroups:
        files.append(os.path.join(config['uBAM'] + '/' + wildcards['sample'] + '.' + readgroup['info']['ID'] + '_unmapped_masked.bam'))
    return files
# comdine all uBAM with masked adpters in one file
# extract fq from unmapped bam files with masked adapters
rule merge_ubams:
    input:
        get_readgroups_unmapped
    output:
        merged_bam = config['uBAM']  + "/{sample}_merged_unmapped.masked.bam",
        f1_masked = temp(config['FQ'] + "/{sample}_masked_1.fq.gz"),
        f2_masked = temp(config['FQ'] + "/{sample}_masked_2.fq.gz")
    threads: config['merge_ubams']['n']
    shell:
        "{samtools} merge -@ {threads} -o {output.merged_bam} {input}  &&  "
        "{gatk} SamToFastq -I {output.merged_bam} -F {output.f1_masked} -F2 {output.f2_masked}"

# convert fq with masked adapters to unmapped CRAM file for storage
rule convert_fq_to_uCRAM:
    input:
        f1 = rules.merge_ubams.output.f1_masked,
        f2 = rules.merge_ubams.output.f2_masked
    output:
        uCRAM = config['CRAM'] + "/{sample}_unmapped_masked.cram"
    threads: config['convert_fq_to_uCRAM']['n']
    shell:
        "{samtools} import -@ {threads} -1 {input.f1} -2 {input.f2} -o {output}"


# cut adapters from inout
# DRAGMAP doesn;t work well with uBAM, so use fq as input
rule cutadapter:
    input:
        get_fastqpaired
    output:
        forr_f=config['FQ'] + "/{sample}.{readgroup}.cut_1.fq.gz",
        rev_f=config['FQ'] + "/{sample}.{readgroup}.cut_2.fq.gz"
    # log file in this case contain some stats about removed seqs
    log:
        cutadapt_log= config['STAT'] + "/{sample}.{readgroup}.cutadapt.log",
    benchmark:
        config['BENCH'] + "/{sample}.{readgroup}.cutadapt.txt"
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
rule align_reads:
    input:
        for_r = rules.cutadapter.output.forr_f,
        rev_r = rules.cutadapter.output.rev_f
    output:
        bam=config['BAM'] + "/{sample}.{readgroup}.bam"
    params:
        ref_dir = config['RES'] + config['ref_dir'],
        # mask bed for current reference genome
        mask_bed = config['RES'] + config['mask_bed']
    threads: config["align_reads"]["n"]
    log:
        dragmap_log=config['LOG'] + '/' + "{sample}.{readgroup}_dragmap.log",
        samtools_fixmate=config['LOG'] + '/' + "{sample}.{readgroup}_samtools_fixamte.log",
        samtools_sort=config['LOG'] + '/' + "{sample}.{readgroup}_samtools_sort.log",
        samtools_markdup=config['LOG'] + '/' + "{sample}.{readgroup}_samtools_markdup.log",
        samtools_index = config['LOG'] + '/' + "{sample}.{readgroup}_samtools_index.log"
    benchmark:
        config['BENCH'] + "/{sample}.{readgroup}.dragmap.txt"
    priority: 15
    shell:
        # "{dragmap} -r {params.ref_dir} -b {input} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --ht-mask-bed {params.mask_bed} --num-threads {threads} 2> {log.dragmap_log} |"
        "{dragmap} -r {params.ref_dir} -1 {input.for_r} -2 {input.rev_r} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --ht-mask-bed {params.mask_bed} --num-threads {threads} 2> {log.dragmap_log} |" 
        "{samtools} fixmate -@ {threads} -m - -  &> {log.samtools_fixmate} | "
        "{samtools} sort -T 'sort_temporary' -@ {threads}  -o {output.bam} &> {log.samtools_sort} &&"
        "{samtools} index -@ {threads} {output.bam} &> {log.samtools_index}"

# merge different readgroups files for same sample
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
        shell("{samtools} merge -@ {threads} -o {output} {inputs} &> {log}")

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
        "{samtools} markdup -f {output.MD_stat} -S -d {params.machine} -@ {threads} {input} {output.mdbams} &> {log.samtools_markdup} && "
        "{samtools} index -@ {threads} {output.mdbams} 2> {log.samtools_index_md}"

        # merge bam files here, add rule, option sort order
        #
        # additional ruleto markdup
        # "{samtools} markdup -@ {threads} - {output.bam} 2> {log.samtools_markdup} && "
        # add -S flag, find what it does 
        # -d for optical duplicates
        # optical duplicates depends on machine
        # zcat fastq/90-492165764/SVDL_KG-010048/KG-010048_R1_001.fastq.gz | grep '@' | awk -F ':' '{print $3}' | head
        # extract part of fq file where machine code is coded
        # https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py#L12-L45
        # dictionary of machines
        # "{samtools} index -@ {threads} {output.bam} 2> {log.samtools_index}"
        # "{dragmap} -r {params.ref_dir} -1 {input[0]} -2 {input[1]} --RGID {wildcards.sample} --RGSM {wildcards.sample} --ht-mask-bed {params.mask_bed} --num-threads {threads} 2> {log.dragmap_log} | "
        # "{samtools} fixmate -@ {threads} -m - -  2> {log.samtools_fixmate} | "
        # "{samtools} sort -T 'sort_temorary' -@ {threads}  2> {log.samtools_sort} | "
        # "{samtools} markdup -@ {threads} - {output.bam} 2> {log.samtools_markdup} && "
        # "{samtools} index -@ {threads} {output.bam} 2> {log.samtools_index}"


ruleorder: sort_back > bamstats_all

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


# check amount of supplementary reads
# if value higher than 0.5% - run additional cleanup steps
# from checkpoint step
# trigger this additional steps only in case if these steps necessary
def check_supp(wildcards):
    with checkpoints.bamstats_all.get(sample=wildcards.sample).output[0].open() as f:
        lines = f.readlines()
        # 4th column (3rd if 0-based) is column with supplementary fraction
        if float((lines[1].split()[3])) >= float(0.005):
            # return cleaned bam in case if clean up necessary
            return rules.sort_back.output.ready_bams
        else:
            # return original MD bam if clean up not necessary
            return rules.markdup.output.mdbams

# calibrate model
# step for HaplotypeCaller in dragen mode
# for better resolution in complex regions
rule CalibrateDragstrModel:
    input:
        check_supp
    output:
        dragstr_model = config['BAM'] + "/{sample}-dragstr.txt"
    priority: 16
    params:
        str_ref = config['RES'] + config['str_ref']
    log: config['LOG'] + '/' + "{sample}_calibratedragstr.log"
    benchmark: config['BENCH'] + "/{sample}_calibrate_dragstr.txt"
    shell:
        "{gatk} CalibrateDragstrModel -R {ref} -I {input} -O {output} -str {params.str_ref} 2>{log}"

# verifybamid
# capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']

# verifybamid has some bugs and require samtools lower version
# to avoid any bugs verifybamid step runs in different conda enviroment
#
rule verifybamid:
    input:
        check_supp
    output:
        VBID_stat = config['STAT'] + '/contam/{sample}_verifybamid.selfSM'
    # end of this file hardcoded in Haplotypecaller and read_contam_w
    threads: config['verifybamid']['n']
    priority: 35
    params:
        VBID_prefix = config['STAT'] + '/{sample}_verifybamid',
        SVD = config['RES'] + config['verifybamid_exome']
    shell:
        """
        set +u
        source ~/bin/start_conda
        PS1=''
        conda info --envs
        source activate verifybamid
        set -u
        verifybamid2 --BamFile {input} --SVDPrefix {params.SVD} --Reference {ref} --DisableSanityCheck --NumThread {threads} --Output {params.VBID_prefix}
        """

def read_contam_w(wildcards):
    # filename = rules.verifybamid.output.VBID_stat
    filename = os.path.join(config['STAT'], 'contam', wildcards['sample'] + '_verifybamid.selfSM')
    with open(filename,'r') as f:
        c = csv.reader(f, delimiter='\t')
        c = iter(c)
        headers = c.__next__()
        data = c.__next__()
        freemix = data[6]
    return freemix

#find SNPs from bams
rule HaplotypeCaller:
    input:
        # check what bam file we need to use (with or without additional cleanup)
        bams = check_supp,
        model = rules.CalibrateDragstrModel.output.dragstr_model,
        contam= rules.verifybamid.output.VBID_stat,
        # command to get path to capture_kit interval list from SAMPLEFILE
        interval = get_capture_kit_interval_list,
    output:
        gvcf="gvcfs/{sample}.g.vcf.gz",
    log:
        HaplotypeCaller=config['LOG'] + '/' + "{sample}_haplotypecaller.log"
    benchmark:
        config['BENCH'] + "/{sample}_haplotypecaller.txt"
    params:
        dbsnp = config['RES'] + config['dbsnp'],
        padding=100,  # extend intervals to this bp
        contam_frac = read_contam_w
    priority: 25
    shell:
        "{gatk} HaplotypeCaller \
                 -R {ref} -L {input.interval} -ip {params.padding} -D {params.dbsnp} -ERC GVCF --contamination {params.contam_frac} \
                 -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
                 -I {input.bams} -O {output.gvcf} \
                  --dragen-mode true --dragstr-params-path {input.model} 2> {log.HaplotypeCaller}"

#############################################################################################
#############################################################################################
# downstream analysis require less computational power and don't support multithreading
# change to another snakemake file and run separatly?
#############################################################################################
#############################################################################################

#Genomics DBImport instead CombineGVCFs
# run parts of chrs instead of full chr (at least fisrt 10)
rule GenomicDBImport:
    input:
        expand("gvcfs/{sample}.g.vcf.gz", sample = sample_names)
    log: config['LOG'] + '/' + "GenomicDBImport.{chrs}.log"
    benchmark: config['BENCH'] + "/GenomicDBImport.{chrs}.txt"
    output:
        dbi=directory("genomicsdb_{chrs}")
    threads: config['GenomicDBImport']['n']
    # params:
        # N_intervals=5,
        # threads=16,
        # padding = 100
    priority: 30
    shell:
        "ls gvcfs/*.g.vcf.gz > gvcfs.list && {gatk} GenomicsDBImport --reader-threads {threads}\
        -V gvcfs.list --intervals {wildcards.chrs}  -R {ref} --genomicsdb-workspace-path {output} \
         --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader 2> {log}"
        # "ls gvcfs/*.g.vcf > gvcfs.list && {gatk} GenomicsDBImport -V gvcfs.list --intervals {chrs}  -R {ref} --genomicsdb-workspace-path {output} \
        #      --max-num-intervals-to-import-in-parallel {params.N_intervals} --reader-threads {params.threads}"

# genotype
rule GenotypeDBI:
    input:
        rules.GenomicDBImport.output.dbi
    output:
        raw_vcfDBI=config['VCF'] + "/Merged_raw_DBI_{chrs}.vcf.gz"
    log: config['LOG'] + '/' + "GenotypeDBI.{chrs}.log"
    benchmark: config['BENCH'] + "/GenotypeDBI.{chrs}.txt"
    params:
        dbsnp = config['RES'] + config['dbsnp']
    priority: 40
    shell:
            "{gatk} GenotypeGVCFs -R {ref} -V gendb://{input} -O {output} -D {params.dbsnp} --intervals {wildcards.chrs} 2> {log}"


rule Mergechrs:
    input:
        expand(config['VCF'] + "/Merged_raw_DBI_{chrs}.vcf.gz", chrs = chrs)
    params:
        vcfs = list(map("-I vcfs/Merged_raw_DBI_{}.vcf.gz".format, chrs))
    log: config['LOG'] + '/' + "Mergechrs.log"
    benchmark: config['BENCH'] + "/Mergechrs.txt"
    output:
        vcf = config['VCF'] + "/ALL_chrs.vcf.gz"
    priority: 45
    shell:
        "{gatk} GatherVcfs {params.vcfs} -O {output} -R {ref} 2> {log} && {gatk} IndexFeatureFile -I {output} "


# VQSR
#select SNPs for VQSR
# SNPs and INDELs require different options
rule SelectSNPs:
    input:
        rules.Mergechrs.output.vcf
    output:
        SNP_vcf=temp(config['VCF'] + "/Merged_SNPs.vcf")
    priority: 50
    log: config['LOG'] + '/' + "SelectSNPs.log"
    benchmark: config['BENCH'] + "/SelectSNPs.txt"
    shell:
        """
        {gatk} SelectVariants \
                --select-type-to-include SNP \
                -V {input} -O {output} 2> {log}
        """

# 1st step of VQSR - calculate scores
rule VQSR_SNP:
    input:
        rules.SelectSNPs.output.SNP_vcf
    output:
        recal_snp=temp(config['VCF'] + "/SNPs_vqsr.recal"),
        tranches_file_snp=temp(config['VCF'] + "/SNPs_vqsr.tranches"),
        r_snp=config['STAT'] + "/SNPs_vqsr_plots.R"
    log: config['LOG'] + '/' + "VQSR_SNP.log"
    benchmark: config['BENCH'] + "/VQSR_SNP.txt"
    params:
        hapmap = config['RES'] + config['hapmap'],
        omni = config['RES'] + config['omni'],
        kilo_g = config['RES'] + config['kilo_g'],
        dbsnp = config['RES'] + config['dbsnp']
    priority: 55
    shell:
        # -an InbreedingCoeff if 10+
        """
        {gatk} VariantRecalibrator\
         -R {ref} -V {input} \
         -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
         -resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
         -resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.kilo_g} \
         -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
         -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP -an InbreedingCoeff\
         --trust-all-polymorphic -AS TRUE\
         -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
         -O {output.recal_snp} \
         --tranches-file {output.tranches_file_snp} \
         --rscript-file {output.r_snp} 2> {log}
        """

rule ApplyVQSR_SNPs:
    input:
        recal_snp=rules.VQSR_SNP.output.recal_snp,
        tranches_snp=rules.VQSR_SNP.output.tranches_file_snp,
        snps_variants=rules.SelectSNPs.output.SNP_vcf
    output:
        recal_vcf_snp=temp(config['VCF'] + "/SNPs_recal_apply_vqsr.vcf")
    log: config['LOG'] + '/' + "Apply_VQSR_SNP.log"
    benchmark: config['BENCH'] + "/Apply_VQSR_SNP.txt"
    params:
        ts_level='99.0'  #ts-filter-level show the "stregnth" of VQSR could be from 90 to 100
    priority: 60
    shell:
        """
        {gatk} ApplyVQSR -R {ref} -mode SNP \
        --recal-file {input.recal_snp} --tranches-file {input.tranches_snp} \
        -O {output} -V {input.snps_variants} -ts-filter-level {params.ts_level} -AS TRUE 2> {log}
        """

# select INDELs for VQSR
rule SelectINDELs:
    input:
        rules.Mergechrs.output.vcf
    output:
        INDEL_vcf=temp(config['VCF'] + "/Merged_INDELs.vcf")
    log: config['LOG'] + '/' + "SelectINDELS.log"
    benchmark: config['BENCH'] + "/SelectINDELs.txt"
    priority: 50
    shell:
        """
        {gatk} SelectVariants \
                --select-type-to-include INDEL \
                -V {input} -O {output} 2> {log}
        """

# 1st step of VQSR - calculate scores
rule VQSR_INDEL:
    input:
        rules.SelectINDELs.output.INDEL_vcf
    output:
        recal_indel=temp(config['VCF'] + "/INDELs_vqsr.recal"),
        tranches_file_indel=temp(config['VCF'] + "/INDELs_vqsr.tranches"),
        r_indel=config['STAT'] + "/INDELs_vqsr_plots.R"
    log: config['LOG'] + '/' + "VQSR_INDEL.log"
    benchmark: config['BENCH'] + "/VQSR_INDEL.txt"
    priority: 55
    params:
        mills = config['RES'] + config['mills'],
        dbsnp_indel = config['RES'] + config['dbsnp_indel']
    shell:
        # -an InbreedingCoeff if 10+
        """
        {gatk} VariantRecalibrator -R {ref} -V {input} \
        -O {output.recal_indel} --tranches-file {output.tranches_file_indel} --rscript-file {output.r_indel} \
         --max-gaussians 4 --trust-all-polymorphic -AS TRUE\
         -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
         -resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mills} \
         -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp_indel} \
         -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode INDEL \
         --trust-all-polymorphic 2> {log}
        """

rule ApplyVQSR_INDEs:
    input:
        recal_indel=rules.VQSR_INDEL.output.recal_indel,
        tranches_indel=rules.VQSR_INDEL.output.tranches_file_indel,
        indel_variants=rules.SelectINDELs.output.INDEL_vcf
    log: config['LOG'] + '/' + "ApplyVQSR_INDELs.log"
    benchmark: config['BENCH'] + "/ApplyVQSR_INDELs.txt"
    output:
        recal_vcf_indel=temp(config['VCF'] + "/INDELs_recal_apply_vqsr.vcf")
    params:
        ts_level='97.0'  #ts-filter-level show the "stregnth" of VQSR could be from 90 to 100
    priority: 60
    shell:
        """
        {gatk} ApplyVQSR -R {ref} -mode INDEL \
        --recal-file {input.recal_indel} --tranches-file {input.tranches_indel} \
        -O {output} -V {input.indel_variants} -ts-filter-level {params.ts_level} -AS TRUE 2> {log}
        """

#combine filtr results
rule combine:
    input:
        snps=rules.ApplyVQSR_SNPs.output.recal_vcf_snp,
        indel=rules.ApplyVQSR_INDEs.output.recal_vcf_indel
    log: config['LOG'] + '/' + "combine.log"
    benchmark: config['BENCH'] + "/combine.txt"
    output:
        filtrVCF=config['VCF'] + "/Merged_after_VQSR.vcf"
    priority: 70
    shell:
        "{gatk} MergeVcfs \
                -I {input.snps} -I {input.indel} -O {output} 2> {log}"

# normalization with bcftools
rule norm:
    input:
        rules.combine.output.filtrVCF
    output:
        normVCF=config['VCF'] + "/Merged_after_VQSR_norm.vcf",
        idx=config['VCF'] + "/Merged_after_VQSR_norm.vcf.idx"
    log: config['LOG'] + '/' + "normalization.log"
    benchmark: config['BENCH'] + "/normalization.txt"
    priority: 80
    shell:
        "{bcftools} norm -f {ref} {input} -m -both -O v | {bcftools} norm -d exact -f {ref} > {output.normVCF} 2> {log} && {gatk} IndexFeatureFile -I {output.normVCF} -O {output.idx} "

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
    threads: config['Basic_stats']['n']
    shell:
        "{gatk} CollectVariantCallingMetrics \
        -R {ref} -I {input} -O stats/BASIC \
        --DBSNP {params.dbsnp} --THREAD_COUNT {threads} 2> {log}"

#hsmetrics
#include off-target metrics
rule HS_stats:
    input:
        vcf = check_supp,
        interval = get_capture_kit_interval_list
    output:
        HS_metrics=config['STAT'] + "/{sample}_hs_metrics"
    log: config['LOG'] + '/' + "HS_stats_{sample}.log"
    benchmark: config['BENCH'] + "/HS_stats_{sample}.txt"
    priority: 99
    params:
        #minimum Base Quality for a base to contribute cov
        #def is 20
        Q=10,
        #minimin Mapping Quality for a read to contribute cov
        #def is 20
        MQ=10,
    shell:
        "{gatk} CollectHsMetrics \
            -I {input.vcf} -R {ref} -BI {input.interval} -TI {input.interval} \
            -Q {params.Q} -MQ {params.MQ} \
            --PER_TARGET_COVERAGE stats/{wildcards.sample}_per_targ_cov \
            -O stats/{wildcards.sample}_hs_metrics 2> {log}"

rule Artifact_stats:
    input:
        vcf = check_supp,
        interval= get_capture_kit_interval_list
    output:
        Bait_bias = config['STAT'] + '/{sample}.bait_bias.bait_bias_summary_metrics',
        Pre_adapter = config['STAT'] + '/{sample}.bait_bias.pre_adapter_summary_metrics'
        # Artifact_matrics = config['STAT'] + "/{sample}.bait_bias.bait_bias_detail_metrics"
    priority: 99
    log: config['LOG'] + '/' + "Artifact_stats_{sample}.log"
    benchmark: config['BENCH'] + "/Artifact_stats_{sample}.txt"
    params:
        # output define prefix, not full filename
        # params.out define prefix and output define whole outputs' filename
        out = config['STAT'] + "/{sample}.bait_bias",
        dbsnp = config['RES'] + config['dbsnp']
    shell:
        "{gatk} CollectSequencingArtifactMetrics -I {input.vcf} -O {params.out} \
        -R {ref} --DB_SNP {params.dbsnp} --INTERVALS {input.interval} 2> log"

rule OXOG_metrics:
    input:
        vcf = check_supp,
        interval= get_capture_kit_interval_list
    output:
        Artifact_matrics = config['STAT'] + "/{sample}.OXOG"
    priority: 99
    log: config['LOG'] + '/' + "OXOG_stats_{sample}.log"
    benchmark: config['BENCH'] + "/OxoG_{sample}.txt"
    params:
        dbsnp = config['RES'] + config['dbsnp']
    shell:
        "{gatk} CollectOxoGMetrics -I {input.vcf} -O {output} -R {ref} \
         --DB_SNP {params.dbsnp} --INTERVALS {input.interval} 2> {log}"

rule samtools_stat:
    input:
        check_supp
    output: samtools_stat = config['STAT'] + "/{sample}_samtools.stat"
    priority: 99
    log: config['LOG'] + '/' + "samtools_{sample}.log"
    benchmark: config['BENCH'] + "/samtools_stat_{sample}.txt"
    threads: config['samtools_stat']['n']
    shell:
        "{samtools} stat -@ {threads} -r {ref} {input} > {output}"

rule samtools_stat_exome:
    input:
        bam = check_supp,
        bed_interval = get_capture_kit_bed
    output: samtools_stat_exome = config['STAT'] + "/{sample}_samtools.exome.stat"
    priority: 99
    log: config['LOG'] + '/' + "samtools_exome_{sample}.log"
    benchmark: config['BENCH'] + "/samtools_stat_exome_{sample}.txt"
    threads: config['samtools_stat']['n']
    shell:
        "{samtools} stat -@ {threads} -t {input.bed_interval} -r {ref} {input.bam} > {output}"




rule bamstats_exome:
    input:
        bam = check_supp,
        bed_interval= get_capture_kit_bed
    output:
        All_exome_stats = config['STAT'] + '/{sample}.bam_exome.tsv'
    threads: config['bamstats_exome']['n']
    params:
        py_stats = config['BAMSTATS']
    shell:
        "samtools view -s 0.05 -h {input.bam} --threads {threads} -L {input.bed_interval} | python3 {params.py_stats} stats > {output}"

# which statistic file we need to use depends on additional clean up steps
# if we trigger addition clean up, so we need to provide bam_all.tsv statistic file after cleanup
# other stat files won't chged a lot after cleanup, so we keep them
def check_supp_stats(wildcards):
    with checkpoints.bamstats_all.get(sample=wildcards).output[0].open() as f:
        lines = f.readlines()
        if float((lines[1].split()[3])) >= float(0.005):
            return os.path.join(config['STAT'] + '/' + wildcards + '.bam_all.additional_cleanup.tsv')
        else:
            return os.path.join(config['STAT'] + '/' + wildcards + '.bam_all.tsv')

def get_quality_stats(wildcards):
    sampleinfo = SAMPLES_BY_FILE[os.path.basename(wildcards['samplefile'])]
    files = []
    samples = list(sampleinfo.keys())
    samples.sort()
    for sample in samples:
        files.append(config['STAT'] + '/' + sample + "_samtools.stat")
        files.append(config['STAT'] + '/' + sample + '_samtools.exome.stat')
        files.append(config['STAT'] + '/contam/' + sample + '_verifybamid.selfSM')
        files.append(check_supp_stats(sample))
        files.append(config['STAT'] + '/' + sample + '.bam_exome.tsv')
        files.append(config['STAT'] + '/' + sample + '.bait_bias.pre_adapter_summary_metrics')
        files.append(config['STAT'] + '/' + sample + '.bait_bias.bait_bias_summary_metrics')
    return files


rule gatherstats:
    input:
        get_quality_stats
    output:
        '{samplefile}.bam_quality.v4.tab'
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
        files.append(config['STAT'] + '/' + sample + '.bait_bias.pre_adapter_summary_metrics')
        files.append(config['STAT'] + '/' + sample + '.bait_bias.bait_bias_summary_metrics')
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


