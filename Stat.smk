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
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()
sample_sex_names = get_sample_sex_combi(SAMPLEINFO)

tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner
# module VCF:
#     snakefile: 'gVCF.smk'


   
def get_ref_by_sex(wildcards):
    if wildcards['sex'] == 'female':
        ref_dir=os.path.join(config['RES'],config['ref_female'])
    else:
        ref_dir=os.path.join(config['RES'],config['ref_male'])

    return ref_dir
   


rule Stat_all:
    input:
        expand("{stats}/{sample_sex}.hs_metrics",sample_sex=sample_sex_names, stats = config['STAT']),
        expand("{stats}/{sample_sex}.OXOG", sample_sex=sample_sex_names, stats = config['STAT']),
        expand("{stats}/{sample_sex}.samtools.stat", sample_sex=sample_sex_names, stats = config['STAT']),
        expand("{samplefile}.oxo_quality.tab", samplefile = SAMPLE_FILES),
        expand("{samplefile}.bam_quality.tab", samplefile = SAMPLE_FILES),
    default_target: True

def get_svd(wildcards):
    """Returns the VerifyBamID SVD file for the sample type of the sample"""
    sinfo = SAMPLEINFO[wildcards['sample']]
    if 'wgs' in sinfo['sample_type']:
        SVD =  os.path.join(config['RES'], config['verifybamid_wgs'])
    else:
        SVD = os.path.join(config['RES'],  config['verifybamid_exome'])
    return SVD


rule verifybamid:
    """Estimates contamination in a sample using the verifybamid2 tool"""
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
    output:
        VBID_stat = config['STAT'] + '/contam/{sample}.{sex}.verifybamid.pca2.selfSM'
    # end of this file hardcoded in Haplotypecaller and read_contam_w
    benchmark: config['BENCH'] + "/{sample}.{sex}.verifybamid.txt"
    priority: 27
    params:
        VBID_prefix = os.path.join(config['STAT'], 'contam/{sample}.{sex}.verifybamid.pca2'),
        SVD = get_svd
    conda: 'envs/verifybamid.yaml'
    resources:
        mem_mb=400,
        n=2
    shell:
        """
        verifybamid2 --BamFile {input.bam} --SVDPrefix {params.SVD} --Reference {ref} --DisableSanityCheck --NumThread {resources.n} --Output {params.VBID_prefix}
        """

def get_capture_kit_interval_list(wildcards):
    """Returns the capture kit interval list file for the sample type of the sample"""
    if SAMPLEINFO[wildcards['sample']]['sample_type'].startswith('illumina_wgs'):
        capture_kit = MERGED_CAPTURE_KIT
    else:
        capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    capture_kit_path = os.path.join(config['RES'], config['kit_folder'], capture_kit + '_hg38.interval_list')
    return capture_kit_path

def get_mem_mb_hs_stats(wildcrads, attempt):
    return (attempt * int(8500))

rule hs_stats:
    """Collects HS metrics for a sample using the gatk CollectHsMetrics tool"""
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output:
        HS_metrics=os.path.join(config['STAT'], "{sample}.{sex}.hs_metrics")
    log: os.path.join(config['LOG'], "HS_stats_{sample}.{sex}.log")
    benchmark: os.path.join(config['BENCH'],  "HS_stats_{sample}.{sex}.txt")
    priority: 99
    params:
        ref=get_ref_by_sex,
        interval = get_capture_kit_interval_list,
        #minimum Base Quality for a base to contribute cov (default=20)
        Q=10,
        #minimum Mapping Quality for a read to contribute cov(default=20)
        MQ=10,
        java_options=config['DEFAULT_JAVA_OPTIONS'],
    conda: "envs/preprocess.yaml"
    resources: mem_mb = get_mem_mb_hs_stats,
               tmpdir = tmpdir,  
               n=1            
    shell:
        """
            gatk  --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" CollectHsMetrics  --TMP_DIR {resources.tmpdir} \
            -I {input.bam} -R {params.ref} -BI {params.interval} -TI {params.interval} \
            -Q {params.Q} -MQ {params.MQ} \
            --PER_TARGET_COVERAGE stats/{wildcards.sample}_per_targ_cov \
            -O stats/{wildcards.sample}_hs_metrics 2> {log}"""

def get_mem_mb_Artifact_stats(wildcrads, attempt):
    return (attempt * int(7000))

rule Artifact_stats:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output:
        Bait_bias = os.path.join(config['STAT'], '{sample}.{sex}.bait_bias.bait_bias_summary_metrics'),
        Pre_adapter = os.path.join(config['STAT'], '{sample}.{sex}.bait_bias.pre_adapter_summary_metrics'),
        Bait_bias_det = os.path.join(config['STAT'],'{sample}.{sex}.bait_bias.bait_bias_detail_metrics'),
        Pre_adapter_det = os.path.join(config['STAT'], '{sample}.{sex}.bait_bias.pre_adapter_detail_metrics'),
        # Artifact_matrics = config['STAT'] + "/{sample}.bait_bias.bait_bias_detail_metrics"
    priority: 99
    log: os.path.join(config['LOG'], "Artifact_stats_{sample}.{sex}.log")
    benchmark: os.path.join(config['BENCH'], "Artifact_stats_{sample}.{sex}.txt")
    params:
        interval = get_capture_kit_interval_list,
        # output define prefix, not full filename
        # params.out define prefix and output define whole outputs' filename
        out = os.path.join(config['STAT'], "{sample}.{sex}.bait_bias"),
        dbsnp = os.path.join(config['RES'], config['dbsnp']),
        java_options=config['DEFAULT_JAVA_OPTIONS'],
    conda: "envs/preprocess.yaml"
    resources: mem_mb = get_mem_mb_Artifact_stats,
                tmpdir= tmpdir,
                n=1
    shell:
        """
            gatk --java-options "-Xmx{resources.mem_mb}M {params.java_options}" CollectSequencingArtifactMetrics  --TMP_DIR {resources.tmpdir} -I {input.bam} -O {params.out} \
        -R {ref} --DB_SNP {params.dbsnp} --INTERVALS {params.interval} 2> {log}"""

def get_mem_mb_OXOG_metrics(wildcrads, attempt):
    return (attempt * int(5500))
rule OXOG_metrics:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output:
        Artifact_matrics = os.path.join(config['STAT'], "{sample}.{sex}.OXOG")
    priority: 99
    log: os.path.join(config['LOG'], "OXOG_stats_{sample}.{sex}.log")
    benchmark: os.path.join(config['BENCH'], "OxoG_{sample}.{sex}.txt")
    params:
        interval = get_capture_kit_interval_list,
        java_options=config['DEFAULT_JAVA_OPTIONS'],
    conda: "envs/preprocess.yaml"
    resources:
                mem_mb = get_mem_mb_OXOG_metrics,
                tmpdir= tmpdir,
                n=1
    shell:
        """gatk  --java-options "-Xmx{resources.mem_mb}M {params.java_options}" CollectOxoGMetrics -I {input.bam} -O {output} -R {ref} \
         --INTERVALS {params.interval} 2> {log}"""

rule samtools_stat:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output: samtools_stat = ensure(os.path.join(config['STAT'], "{sample}.{sex}.samtools.stat"), non_empty=True)
    priority: 99
    log: os.path.join(config['LOG'], "samtools_{sample}.{sex}.log")
    benchmark: os.path.join(config['BENCH'], "samtools_stat_{sample}.{sex}.txt")
    resources:
        mem_mb=100,
        n=1
    conda: "envs/preprocess.yaml"
    shell:
        "samtools stat -@ {resources.n} -r {ref} {input.bam} > {output}"

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
        bai= rules.markdup.output.mdbams_bai
    output: samtools_stat_exome = ensure(os.path.join(config['STAT'], "{sample}.{sex}.samtools.exome.stat"), non_empty=True)
    priority: 99
    params:
        bed_interval = get_capture_kit_bed
    log: os.path.join(config['LOG'], "samtools_exome_{sample}.{sex}.log")
    benchmark: os.path.join(config['BENCH'],"samtools_stat_exome_{sample}.{sex}.txt")
    resources:
        mem_mb=100,
        n=1
    conda: "envs/preprocess.yaml"
    shell:
        "samtools stat -@ {resources.n} -t {params.bed_interval} -r {ref} {input.bam} > {output}"

rule bamstats_all:
    input:
        # very annoying bug here
        # if in input 2 functions and one of them is check chekpoint (rules.markdup.output.mdbams and get_capture_kit_bed here was as example)
        # first command has not been executed
        # and in shell wildcard (instead of iutput of function) has putted
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output:
        All_exome_stats = ensure(os.path.join(config['STAT'],'{sample}.{sex}.bam_all.tsv'), non_empty=True)
    benchmark: os.path.join(config['BENCH'],"bamstats_all_{sample}.{sex}.txt")
    params:
        py_stats = srcdir(config['BAMSTATS'])
    resources:
        mem_mb=250,
        n=1        
    conda: "envs/pypy.yaml"
    shell:
        "samtools view -s 0.05 -h {input.bam} --threads {resources.n}  | pypy {params.py_stats} stats > {output}"

rule bamstats_exome:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output:
        All_exome_stats = ensure(os.path.join(config['STAT'], '{sample}.{sex}.bam_exome.tsv'),  non_empty=True)
    benchmark: os.path.join(config['BENCH'],"bamstats_exome_{sample}.{sex}.txt")   
    resources:
        mem_mb=250,
        n=1
    params:
        py_stats = srcdir(config['BAMSTATS']),
        bed_interval= get_capture_kit_bed,
    conda: "envs/pypy.yaml"
    shell:
        "samtools view -s 0.05 -h {input.bam} --threads {resources.n} -L {params.bed_interval} | pypy {params.py_stats} stats > {output}"


def get_quality_stats(wildcards):
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    files = []
    samples_sex = [sample + '.' + ('female' if info['sex'] == 'F' else 'male') for sample, info in sampleinfo.items()]
    
    samples_sex.sort()
    for sample_sex in samples_sex:
        files.append(os.path.join(config['STAT'],  f"{sample_sex}.samtools.stat"))
        files.append(os.path.join(config['STAT'],  f'{sample_sex}.samtools.exome.stat'))
        files.append(os.path.join(config['STAT'], 'contam',  f'{sample_sex}.verifybamid.pca2.selfSM'))
        files.append(os.path.join(config['STAT'], f'{sample_sex}.bam_all.tsv'))
        files.append(os.path.join(config['STAT'], f'{sample_sex}.bam_exome.tsv'))
        files.append(os.path.join(config['STAT'], f'{sample_sex}.bait_bias.pre_adapter_summary_metrics'))
        files.append(os.path.join(config['STAT'], f'{sample_sex}.bait_bias.bait_bias_summary_metrics'))
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
    resources:
        n=1,
        mem_mb=10000        
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
    samples_sex = [sample + '.' + ('female' if info['sex'] == 'F' else 'male') for sample, info in sampleinfo.items()]
    samples_sex.sort()
    for sample_sex in samples_sex:
        files.append(os.path.join(config['STAT'], f'{sample_sex}.bait_bias.pre_adapter_detail_metrics'))
        files.append(os.path.join(config['STAT'],f'{sample_sex}.bait_bias.bait_bias_detail_metrics'))
    return files

rule gatherosostats:
    input:
        get_oxo_stats
    output:
        '{samplefile}.oxo_quality.tab'
    benchmark: os.path.join(config['BENCH'],"{samplefile}_gatherOXOstat.txt")
    resources:
        n=1,
        mem_mb=10000
    run:
        sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
        samples = list(sampleinfo.keys())
        samples.sort()

        rinput = numpy.array(input).reshape(len(samples),int(len(input) / len(samples)))
        pre_adapter, bait_bias = zip(*rinput)

        header, data = read_stats.combine_oxo_stats(samples,pre_adapter,bait_bias)
        read_stats.write_tsv(str(output),header,data)

