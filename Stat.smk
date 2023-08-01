from common import *
import read_stats
import utils
from scripts.pca import *
import pandas as pd

wildcard_constraints:
    sample="[\w\d_\-@]+"


module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

def sampleinfo(SAMPLEINFO, sample, checkpoint=False):
    """If samples are on tape, we do not have sample readgroup info.
    That is, the 'readgroups' field is empty.

    This function first checks if the readgroup info is available on disk,
    in the file SAMPLEINFODIR/<sample>.dat. 

    Alternatively, the function injects a checkpoint rule to load this readgroup info.
    """

    sinfo = SAMPLEINFO[sample]
    if not 'readgroups' in sinfo:
        rgpath = pj(SAMPLEINFODIR, sample + ".dat")
        if os.path.exists(rgpath):
            xsample = utils.load(rgpath)
        elif checkpoint: 
            #no readgroup info yet
            filename = checkpoints.get_readgroups.get(sample=sample).output[0]
            xsample = utils.load(filename)
        sinfo = sinfo.copy()
        sinfo['readgroups'] = xsample['readgroups']
        sinfo['alternative_names'] = sinfo.get('alternative_names',set()).union(xsample['alternative_names'])
        SAMPLEINFO[sample] = sinfo
    return sinfo



module Tools:
    snakefile: 'Tools.smk'
    config: config
use rule BedToIntervalList from Tools
   

rule Stat_all:
    input:
        expand("{samplefile}.oxo_quality.tab", samplefile = SAMPLE_FILES),
        expand("{samplefile}.bam_quality.tab", samplefile = SAMPLE_FILES),
        expand("{samplefile}.bam_rg_quality.tab", samplefile = SAMPLE_FILES),
        expand("{samplefile}.sex_chrom.tab", samplefile = SAMPLE_FILES),
        expand("{samplefile}.coverage.hdf5", samplefile = SAMPLE_FILES),
        expand(pj(STAT, "{sample}.done"), sample = sample_names)
    default_target: True        
    


rule stat_sample_done:
    input:
        pj(STAT, "{sample}.hs_metrics"),
        pj(STAT, "{sample}.samtools.stat"),
        pj(STAT, '{sample}.samtools.exome.stat'),
        pj(STAT, 'contam',  '{sample}.verifybamid.pca2.selfSM'),
        pj(STAT, '{sample}.bam_all.tsv'),
        pj(STAT, '{sample}.bam_exome.tsv'),
        pj(STAT, '{sample}.pre_adapter_summary_metrics'),
        pj(STAT, '{sample}.bait_bias_summary_metrics'),
        pj(STAT, '{sample}.pre_adapter_detail_metrics'),
        pj(STAT,'{sample}.bait_bias_detail_metrics'),
        pj(STAT,'cov', '{sample}.regions.bed.gz')
    output:
        cram = touch(pj(STAT, "{sample}.done"))    
    resources:
        n=1,
        mem_mb=50        



rule coverage:
    """Estimates contamination in a sample using the verifybamid2 tool"""
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai        
    output:
        pj(STAT, 'cov', '{sample}.regions.bed.gz'),
        pj(STAT, 'cov', '{sample}.regions.bed.gz.csi'),
        temp(pj(STAT, 'cov', '{sample}.mosdepth.global.dist.txt')),
        temp(pj(STAT, 'cov', '{sample}.mosdepth.summary.txt'))    
    benchmark: BENCH + "/{sample}.mosdepth.txt"
    priority: 27
    params:
        bed = WINDOWS,
        prefix = pj(STAT, 'cov', '{sample}')
    resources:
        mem_mb=2200,
        n="1.8"
    conda: CONDA_MOSDEPTH
    shell:
        """
            mkdir -p `dirname {output[0]}`
            mosdepth  --threads 2 -b {params.bed} --no-per-base {params.prefix} {input.bam}
        """

def get_regions(wildcards):
    samples = list(SAMPLEFILE_TO_SAMPLES[wildcards['samplefile']])
    samples.sort()
    return [pj(STAT, 'cov', '{sample}.regions.bed.gz'.format(sample=sample)) for sample in samples]

def get_stats_samtools(wildcards):
    samples = list(SAMPLEFILE_TO_SAMPLES[wildcards['samplefile']])
    samples.sort()
    return [pj(STAT, '{sample}.samtools.stat'.format(sample=sample)) for sample in samples]

rule gather_coverage:
    input:
        cov=get_regions,
        mapped=get_stats_samtools
    output:
        hdf5=pj('{samplefile}.coverage.hdf5')
    # resources:
        n=1,
        mem_mb=500 
    run:
        samples = list(SAMPLEFILE_TO_SAMPLES[wildcards['samplefile']])
        samples.sort()
        annotation = WINDOWS_ANNOTATED

        read_stats.write_coverage_to_hdf5(annotation, samples, list(input.mapped), list(input.cov), output.hdf5)


rule clustering_samples:
    input: expand('{samplefile}.coverage.hdf5',  samplefile = SAMPLE_FILES)
    output: pj(STAT, 'PCA_stat.txt')
    conda: CONDA_PCA
    run:
        shell("touch {output}")
        df_with_all_samplefiles = pd.DataFrame()
        for SF in SAMPLE_FILES:
            df_with_all_samplefiles, num_samples = load_hdf5_data('{SF}.coverage.hdf5', df_with_all_samplefiles, 'coverage' )
            with open(pj(STAT, 'PCA_stat.txt'), 'a') as f:
                print(f'Loaded {num_samples} samples from {SF}', file=f)
        PCA_t, cum_exp_var = perform_pca(df_with_all_samplefiles, n_comp= 3)
        best_dict = get_scores_and_labels(combinations, data=PCA_t)
        PCA_t['cluster'] = best_dict['best_labels']
        with open(pj(STAT,'PCA_stat.txt'),'a') as f:
            print(f"Cumulative variance = {cum_exp_var}",file=f)
            print(best_dict, file=f)
        for SF in SAMPLE_FILES:
            add_cluster_info_to_hdf5(SF, clustred_data=PCA_t)



def get_svd(wildcards):
    """Returns the VerifyBamID SVD file for the sample type of the sample"""
    sinfo = SAMPLEINFO[wildcards['sample']]
    return VERIFYBAMID_WGS if 'wgs' in sinfo['sample_type'] else VERIFYBAMID_EXOME



rule verifybamid:
    """Estimates contamination in a sample using the verifybamid2 tool"""
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
        validated_sex = rules.get_validated_sex.output.yaml
    output:
        VBID_stat = STAT + '/contam/{sample}.verifybamid.pca2.selfSM'
    # end of this file hardcoded in Haplotypecaller and read_contam_w
    benchmark: pj(BENCH, "{sample}.verifybamid.txt")
    priority: 27
    params:
        ref=get_ref_by_validated_sex,
        #minimum Base Quality for a base to contribute cov (default=20)
        VBID_prefix = pj(STAT, 'contam/{sample}.verifybamid.pca2'),
        SVD = get_svd
    resources:
        mem_mb=400,
        n="1.5"
    conda: CONDA_VERIFYBAMID
    shell:
        """verifybamid2 --BamFile {input.bam} --SVDPrefix {params.SVD} --Reference {params.ref} --DisableSanityCheck --NumThread 2 --Output {params.VBID_prefix}"""

def get_capture_kit_interval_list(wildcards):
    """Returns the capture kit interval list file for the sample type of the sample"""
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        capture_kit = MERGED_CAPTURE_KIT_IVL
    else:
        capture_kit = pj(INTERVALS_DIR, SAMPLEINFO[wildcards['sample']]['capture_kit'] + '.interval_list')
    return capture_kit

rule hs_stats:
    """Collects HS metrics for a sample using the gatk CollectHsMetrics tool"""
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
        interval = MERGED_CAPTURE_KIT_IVL,
        validated_sex = rules.get_validated_sex.output.yaml,
        targets=TARGETS_IVL,
    output:
        HS_metrics=pj(STAT, "{sample}.hs_metrics")
    benchmark: pj(BENCH,  "HS_stats_{sample}.txt")
    priority: 99
    params:
        ref=get_ref_by_validated_sex,
        #minimum Base Quality for a base to contribute cov (default=20)
        Q=10,
        #minimum Mapping Quality for a read to contribute cov(default=20)
        MQ=10
    resources: mem_mb = lambda wildcards, attempt: attempt * 3100,
               tmpdir = tmpdir,  
               n=1            
    conda: CONDA_VCF               
    shell:
        """gatk  --java-options "-Xmx{resources.mem_mb}M  {DEFAULT_JAVA_OPTIONS}" CollectHsMetrics  --TMP_DIR {resources.tmpdir} \
            -I {input.bam} -R {params.ref} -BI {input.interval} -TI {input.targets} \
            -Q {params.Q} -MQ {params.MQ} \
            -O stats/{wildcards.sample}.hs_metrics"""


rule Artifact_stats:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
        interval = get_capture_kit_interval_list,
        validated_sex = rules.get_validated_sex.output.yaml
    output:
        Bait_bias = pj(STAT, '{sample}.bait_bias_summary_metrics'),
        Pre_adapter = pj(STAT, '{sample}.pre_adapter_summary_metrics'),
        Bait_bias_det = pj(STAT,'{sample}.bait_bias_detail_metrics'),
        Pre_adapter_det = pj(STAT, '{sample}.pre_adapter_detail_metrics'),
    priority: 99
    log: pj(LOG, "Artifact_stats_{sample}.log")
    benchmark: pj(BENCH, "Artifact_stats_{sample}.txt")
    params:
        ref=get_ref_by_validated_sex,
        #minimum Base Quality for a base to contribute cov (default=20)
        # output define prefix, not full filename
        # params.out define prefix and output define whole outputs' filename
        out = pj(STAT, "{sample}")
    resources: mem_mb = lambda wildcards, attempt: attempt * 2750,
                tmpdir= tmpdir,
                n=1
    conda: CONDA_VCF
    shell:
        """gatk --java-options "-Xmx{resources.mem_mb}M {DEFAULT_JAVA_OPTIONS}" CollectSequencingArtifactMetrics  --TMP_DIR {resources.tmpdir} -I {input.bam} -O {params.out} \
                    -R {params.ref} --DB_SNP {DBSNP} --INTERVALS {input.interval} 2> {log}"""


rule OXOG_metrics:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
        interval = get_capture_kit_interval_list,
        validated_sex = rules.get_validated_sex.output.yaml
    output:
        Artifact_matrics = pj(STAT, "{sample}.OXOG")
    priority: 99
    log: pj(LOG, "OXOG_stats_{sample}.log")
    benchmark: pj(BENCH, "OxoG_{sample}.txt")
    params:
        ref=get_ref_by_validated_sex,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2500,
        tmpdir= tmpdir,
        n=1
    conda: CONDA_VCF
    shell:
       """gatk  --java-options "-Xmx{resources.mem_mb}M {DEFAULT_JAVA_OPTIONS}" CollectOxoGMetrics -I {input.bam} -O {output} -R {params.ref} \
         --INTERVALS {input.interval} 2> {log}"""

rule samtools_stat:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
        validated_sex = rules.get_validated_sex.output.yaml
    output: samtools_stat = ensure(pj(STAT, "{sample}.samtools.stat"), non_empty=True)
    priority: 99
    log: pj(LOG, "samtools_{sample}.log")
    benchmark: pj(BENCH, "samtools_stat_{sample}.txt")
    resources:
        mem_mb=130,
        n=1
    conda: CONDA_MAIN        
    params:
            ref=get_ref_by_validated_sex
    shell:
        "samtools stat -@ {resources.n} -r {params.ref} {input.bam} > {output}"


# extract info about capture kit from SAMPLEFILE
# assume that all kits bed and interval_list files are existing and download to res folder
def get_capture_kit_bed(wildcards):
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        capture_kit = MERGED_CAPTURE_KIT_BED
    elif SAMPLEINFO[wildcards['sample']]['capture_kit'] == '':
        capture_kit = MERGED_CAPTURE_KIT_BED
    else:
        capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit'] + '.bed'

    return pj(INTERVALS_DIR, capture_kit)
    
rule samtools_stat_exome:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
        validated_sex = rules.get_validated_sex.output.yaml
    output: samtools_stat_exome = ensure(pj(STAT, "{sample}.samtools.exome.stat"), non_empty=True)
    priority: 99
    params:
        ref=get_ref_by_validated_sex,
        bed_interval = get_capture_kit_bed
    log: pj(LOG, "samtools_exome_{sample}.log")
    benchmark: pj(BENCH,"samtools_stat_exome_{sample}.txt")
    resources:
        mem_mb=130,
        n=1
    conda: CONDA_MAIN
    shell:
        "samtools stat -@ {resources.n} -t {params.bed_interval} -r {params.ref} {input.bam} > {output}"

rule bamstats_all:
    input:
        # very annoying bug here
        # if in input 2 functions and one of them is check chekpoint (rules.markdup.output.mdbams and get_capture_kit_bed here was as example)
        # first command has not been executed
        # and in shell wildcard (instead of iutput of function) has putted
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
    output:
        All_exome_stats = ensure(pj(STAT,'{sample}.bam_all.tsv'), non_empty=True)
    benchmark: pj(BENCH,"bamstats_all_{sample}.txt")
    params:
        py_stats = srcdir(BAMSTATS)
    resources:
        mem_mb=150,
        n=1        
    conda: CONDA_PYPY
    shell:
        "samtools view -s 0.05 -h {input.bam} --threads {resources.n}  | pypy {params.py_stats} stats > {output}"

rule bamstats_exome:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output:
        All_exome_stats = ensure(pj(STAT, '{sample}.bam_exome.tsv'),  non_empty=True)
    benchmark: pj(BENCH,"bamstats_exome_{sample}.txt")   
    resources:
        mem_mb=150,
        n=1
    params:
        py_stats = srcdir(BAMSTATS),
        bed_interval= get_capture_kit_bed,
    conda: CONDA_PYPY
    shell:
        "samtools view -s 0.05 -h {input.bam} --threads {resources.n} -L {params.bed_interval} | pypy {params.py_stats} stats > {output}"


def get_quality_stats(wildcards):
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    return [pj(STAT, f'{sample}.done') for sample in sampleinfo.keys()]

rule gatherstats:
    # keep in mind, that samtools_stat create file even if it it's finished with error or you force to stop it
    # if you force to stop samtools_stat delete all output to prevent errors
    # rm -r stats/*samtools*
    input:
        get_quality_stats
    benchmark: pj(BENCH,"{samplefile}_gatherstat.txt")
    output:
        '{samplefile}.bam_quality.tab'
    resources:
        n=1,
        mem_mb=1000
    run:
        sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
        samples = list(sampleinfo.keys())
        samples.sort()

        stats = [pj(STAT, f"{sample}.samtools.stat") for sample in samples]
        exome_stats = [pj(STAT, f"{sample}.samtools.exome.stat") for sample in samples]
        vpca2 = [pj(STAT, 'contam', f"{sample}.verifybamid.pca2.selfSM") for sample in samples]
        bam_extra_all = [pj(STAT, f"{sample}.bam_all.tsv") for sample in samples]
        bam_extra_exome = [pj(STAT, f"{sample}.bam_exome.tsv") for sample in samples]
        pre_adapter = [pj(STAT, f"{sample}.pre_adapter_summary_metrics") for sample in samples]
        bait_bias = [pj(STAT, f"{sample}.bait_bias_summary_metrics") for sample in samples]
        hs_stats = [pj(STAT, f"{sample}.hs_metrics") for sample in samples]


        header, data = read_stats.combine_quality_stats(samples,stats,exome_stats,vpca2,bam_extra_all,bam_extra_exome,pre_adapter,bait_bias, hs_stats)
        read_stats.write_tsv(str(output),header,data)

rule gather_rg_stats:
    input:
        get_quality_stats
    benchmark: pj(BENCH,"{samplefile}_gatherrgstat.txt")
    output:
        '{samplefile}.bam_rg_quality.tab'
    resources:
        n=1,
        mem_mb=1000
    run:
        smsinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
        sample_readgroups = []
        for sample, sinfo in smsinfo.items():
            x = sampleinfo(SAMPLEINFO, sample)
            for readgroup in x['readgroups']:
                sample_readgroups.append((sample, readgroup['info']['ID']))

        sample_readgroups.sort()

        aremoval = [pj(STAT, f"{s_rg[0]}.{s_rg[1]}.adapter_removal.log") for s_rg in sample_readgroups]
        aidentify= [pj(STAT, f"{s_rg[0]}.{s_rg[1]}.fastq.adapters") for s_rg in sample_readgroups]
        mergestats = [pj(STAT, f"{s_rg[0]}.{s_rg[1]}.merge_stats.tsv") for s_rg in sample_readgroups]
        dragmap_stats = [pj(STAT, f"{s_rg[0]}.{s_rg[1]}.dragmap.log") for s_rg in sample_readgroups]
        dechimer_stats = [pj(STAT, f"{s_rg[0]}.{s_rg[1]}.dechimer_stats.tsv") for s_rg in sample_readgroups]

        header, data = read_stats.combine_rg_quality_stats(sample_readgroups, aremoval, aidentify, mergestats, dragmap_stats, dechimer_stats)
        read_stats.write_tsv(str(output),header,data)


def get_oxo_stats(wildcards):
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    return [pj(STAT, f'{sample}.done') for sample in sampleinfo.keys()]


rule gatheroxostats:
    input:
        get_oxo_stats
    output:
        '{samplefile}.oxo_quality.tab'
    benchmark: pj(BENCH,"{samplefile}_gatherOXOstat.txt")
    resources:
        n=1,
        mem_mb=1000
    run:
        sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
        samples = list(sampleinfo.keys())
        samples.sort()
        

        pre_adapter = [pj(STAT,f'{sample}.pre_adapter_detail_metrics') for sample in samples]
        bait_bias = [pj(STAT,f'{sample}.bait_bias_detail_metrics') for sample in samples]


        header, data = read_stats.combine_oxo_stats(samples,pre_adapter,bait_bias)
        read_stats.write_tsv(str(output),header,data)


def get_sex_stats(wildcards):
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    return [pj(KMER, f'{sample}.result.yaml') for sample in sampleinfo.keys()]

rule gathersexstats:
    input:
        get_sex_stats
    benchmark: pj(BENCH,"{samplefile}_gatherSexstat.txt")
    output:
        '{samplefile}.sex_chrom.tab'
    resources:
        n=1,
        mem_mb=1000
    run:
        sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
        samples = list(sampleinfo.keys())
        samples.sort()

        kmer_stats = [pj(KMER, f"{sample}.result.yaml") for sample in samples]
        sex_reported = [sampleinfo[sample]['sex'] for sample in samples]

        header, data = read_stats.combine_sex_stats(samples,kmer_stats, sex_reported)
        read_stats.write_tsv(str(output),header,data)
