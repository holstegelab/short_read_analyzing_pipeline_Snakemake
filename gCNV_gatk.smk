from common import *
import h5py
import numpy as np

onsuccess: shell(f"rm -fr {pj(LOG, 'gCNV_gatk')}/*")
wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    cohort=r"\d+",
    index=r"\d+",
    scatter=r"\d{4}_of_\d+",

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

module Stat:
    snakefile: 'Stat.smk'
    config: config
use rule * from Stat

module PCA:
    snakefile: 'PCA.smk'
    config: config
use rule * from PCA

module DownloadAndExtractBam:
    snakefile: 'Download_and_extract_bam.smk'
    config: config
use rule download_encrypted_cram_from_tape from DownloadAndExtractBam
use rule decrypt_cram from DownloadAndExtractBam
use rule cram_to_markdup_bam from DownloadAndExtractBam

module Tools:
    snakefile: 'Tools.smk'
    config: config
# use rule * from Tools

ruleorder: cram_to_markdup_bam > markdup

GCNV_WGS_BIN_LENGTH = int(config.get("gcnv_wgs_bin_length", 1000))
GCNV_WGS_SCATTER_COUNT = int(config.get("gcnv_wgs_scatter_count", 300))
GCNV_INTERVALS_DIR = pj(GATK_gCNV, "intervals")
GCNV_WGS_INTERVALS = pj(GCNV_INTERVALS_DIR, f"wgs_{GCNV_WGS_BIN_LENGTH}bp.preprocessed.interval_list")
GCNV_WGS_ANNOTATED_INTERVALS = pj(GCNV_INTERVALS_DIR, f"wgs_{GCNV_WGS_BIN_LENGTH}bp.annotated.tsv")


def get_groups_from_hdf5(hdf5_file):
    unique_groups = set()
    with h5py.File(hdf5_file, 'r') as f:
        groups = f['coverage']['cluster_3_cor_cov']
        for group in groups:
            if group != -1:
                unique_groups.add(group)
    return (unique_groups)

hdf5_files = [f"{SF}.coverage.hdf5" for SF in SAMPLE_FILES]
groups = set()
for hdf5_file in hdf5_files:
    groups.update(get_groups_from_hdf5(hdf5_file))

def get_samples_in_group(cohort, hdf5_files = hdf5_files):
    samples = []
    for hdf5_file in hdf5_files:
        with h5py.File(hdf5_file, 'r') as f:
            data_group = f['coverage']
            sample_list_per_hdf5_file = [s.decode('utf-8') for s in data_group['samples'][()]]
            cohort_mask = np.equal(data_group['cluster_3_cor_cov'], int(cohort))
            cohort_samples = [sample_list_per_hdf5_file[i] for i, mask in enumerate(cohort_mask) if mask]
            samples.extend(cohort_samples)
    samples.sort()
    idx = list(range(len(samples)))
    return samples, idx

paths_output = []
for cohort in groups:
    sample_names, idxs = get_samples_in_group(cohort=cohort)
    for sample_name, idx in zip(sample_names, idxs):
        path = pj(GATK_gCNV, f'GENOTYPED_CALLS_intervals_{cohort}', f'COHORT_{cohort}_SAMPLE_{sample_name}_{idx}.vcf.gz')
        paths_output.append(path)



def input_func(wildcards):
    cohort = wildcards['cohort']
    samples_in_group, _ = get_samples_in_group(cohort=cohort)
    return expand(
        '{gatk_gcnv}/Read_counts_hdf5/{cohort}/{sample}_readcounts.hdf5',
        gatk_gcnv=GATK_gCNV,
        cohort=cohort,
        sample=samples_in_group,
        allow_missing=True
    )


def sample_list_per_cohort(wildcards):
    cohort = wildcards['cohort']
    samples_in_group, _ = get_samples_in_group(cohort)
    return expand(
        ' -I {gatk_gcnv}/Read_counts_hdf5/{cohort}/{sample}_readcounts.hdf5 ',
        gatk_gcnv=GATK_gCNV,
        cohort=cohort,
        sample=samples_in_group,
        allow_missing=True
    )

def original_markdup_bam(wildcards):
    samplefile = os.path.basename(SAMPLEINFO[wildcards.sample]["samplefile"])
    return pj(get_samplefile_folder(samplefile), BAM, f"{wildcards.sample}.markdup.bam")

def original_markdup_bai(wildcards):
    return original_markdup_bam(wildcards) + ".bai"

def use_original_markdup_bam(wildcards):
    bam = original_markdup_bam(wildcards)
    bai = original_markdup_bai(wildcards)
    return os.path.exists(bam) and os.path.exists(bai)

def markdup_bam_input(wildcards):
    if use_original_markdup_bam(wildcards):
        return original_markdup_bam(wildcards)
    return pj(BAM, f"{wildcards.sample}.markdup.bam")

def markdup_bai_input(wildcards):
    if use_original_markdup_bam(wildcards):
        return original_markdup_bai(wildcards)
    return pj(BAM, f"{wildcards.sample}.markdup.bam.bai")

def generate_scatter_list(start, end):
    string_list = [f"{i:04d}_of_{end}" for i in range(int(start), int(end) + 1)]
    return string_list

gcnv_scatter_shards = generate_scatter_list('0001', GCNV_WGS_SCATTER_COUNT)

rule gCNV_gatk_all:
    input:
        paths_output
    default_target: True

rule preprocess_wgs_gcnv_intervals:
    input:
        intervals = MAIN_CHRS_BED
    output:
        intervals = GCNV_WGS_INTERVALS
    params:
        java = java_cnv,
        gatk = gatk_cnv,
        ref = REF_MALE,
        bin_length = GCNV_WGS_BIN_LENGTH,
        output_dir = GCNV_INTERVALS_DIR
    conda: CONDA_GATK_CNV
    log: pj(LOG, "gCNV_gatk", f"preprocess_wgs_{GCNV_WGS_BIN_LENGTH}bp_intervals.log")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.java} -jar {params.gatk} PreprocessIntervals -R {params.ref} -L {input.intervals} -imr OVERLAPPING_ONLY --bin-length {params.bin_length} -O {output.intervals} 2> {log}
        """

rule collect_read_counts:
    input: bam=markdup_bam_input,
            bai=markdup_bai_input,
            intervals = rules.preprocess_wgs_gcnv_intervals.output.intervals
    output: ReadCounts = ensure(pj(GATK_gCNV, 'Read_counts_hdf5', '{cohort}', '{sample}_readcounts.hdf5'), non_empty=True)
    params: java = java_cnv,
            gatk = gatk_cnv,
            ref = REF_MALE
    conda: CONDA_GATK_CNV
    resources:
            n = 1,
            mem_mb = 8500 # if tool doesn't have enough RAM it ends without error and create corrupted file. Therefore here pretty big RAM (mean usage much lower)
    log: pj(LOG,"gCNV_gatk", '{sample}.{cohort}.collectreadcounts_gatkcnv.log')
    benchmark: pj(BENCH, '{sample}.{cohort}.collectreadcounts_gatkcnv.txt')
    shell: """
            {params.java} -jar {params.gatk} CollectReadCounts -L {input.intervals} -R {params.ref} -imr OVERLAPPING_ONLY -I {input.bam} -O {output.ReadCounts} 2> {log}
            """

rule annotateintervals:
    input: interval = rules.preprocess_wgs_gcnv_intervals.output.intervals
    output: annotated_tsv = GCNV_WGS_ANNOTATED_INTERVALS
    conda: CONDA_GATK_CNV
    params: java= java_cnv,
            gatk=gatk_cnv,
            ref= REF_MALE,
            hard_mappability_track = HARD_MAPPABILITY_TRACK,
    shell:
            """
            {params.java} -jar {params.gatk} AnnotateIntervals -L {input.interval} -R {params.ref} -imr OVERLAPPING_ONLY -O {output.annotated_tsv} --mappability-track {params.hard_mappability_track}
            """

rule filterintervals:
    input: samples = input_func,
            annotation = rules.annotateintervals.output.annotated_tsv,
            intervals = rules.preprocess_wgs_gcnv_intervals.output.intervals,
    output: filtered_intervals = pj(GATK_gCNV, 'filtered_intervals', '{cohort}_filtered.interval_list')
    params: inputs = sample_list_per_cohort,
            java= java_cnv,
            gatk=gatk_cnv,
            PAR_and_CENTROMERIC = PAR_and_CENTROMERIC,
            cut_off_samples = '80',
            cut_off_count_threshold = '10',
    conda: CONDA_GATK_CNV
    log: pj(LOG,"gCNV_gatk", '{cohort}.filterintervals_gatkcnv.log')
    shell: """
            {params.java} -jar {params.gatk} FilterIntervals --annotated-intervals {input.annotation} -L {input.intervals} {params.inputs} -imr OVERLAPPING_ONLY -O {output.filtered_intervals} -XL {params.PAR_and_CENTROMERIC} --low-count-filter-percentage-of-samples {params.cut_off_samples} --low-count-filter-count-threshold {params.cut_off_count_threshold}  2> {log}
            """

rule make_scatters:
    input: rules.filterintervals.output.filtered_intervals
    output: scatterd = expand(pj(GATK_gCNV, 'filtered_intervals', 'cohort-{cohort}', 'temp_{scatter}', 'scattered.interval_list'), scatter = gcnv_scatter_shards, allow_missing = True)
    conda: CONDA_GATK_CNV
    log: pj(LOG,"gCNV_gatk",'{cohort}.make_scatters.log')
    params:
            java = java_cnv,
            gatk = gatk_cnv,
            scatter_count = GCNV_WGS_SCATTER_COUNT,
            output_dir = pj(GATK_gCNV, 'filtered_intervals', 'cohort-{cohort}')
    shell:
        """
        mkdir -p {params.output_dir}
        {params.java} -jar {params.gatk} IntervalListTools -I {input} --OUTPUT {params.output_dir} --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_COUNT {params.scatter_count}
        """


rule DetermineGCP:
    input: samples = input_func,
            intervals = pj(GATK_gCNV, 'filtered_intervals', '{cohort}_filtered.interval_list')
    output: CPC = pj(GATK_gCNV,  '{cohort}-model', 'contig_ploidy_prior.tsv'),
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            contig_ploydi_priors = PL_PR_TABLE,
            output_dir = GATK_gCNV
    conda: CONDA_GATK_CNV
    resources:
            n = 1
    shell: """
            export OMP_NUM_THREADS={resources.n}
            {params.java} -jar {params.gatk} DetermineGermlineContigPloidy  --output-prefix {wildcards.cohort} --output {params.output_dir} {params.inputs} -L {input.intervals} -imr OVERLAPPING_ONLY --contig-ploidy-priors {params.contig_ploydi_priors}
            """

rule GermlineCNVCaller:
    input:  samples = input_func,
            scatters = pj(GATK_gCNV, 'filtered_intervals', 'cohort-{cohort}', 'temp_{scatter}', 'scattered.interval_list'),
            contig_ploudi_calls = (pj(GATK_gCNV,  '{cohort}-model', 'contig_ploidy_prior.tsv')),
    output: models = pj(GATK_gCNV,'{cohort}_scatter_{scatter}','scatterd_{cohort}_{scatter}-model', 'calling_config.json'),
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            contig_ploydi_calls = (pj(GATK_gCNV,  '{cohort}-calls')),
            theano_complie_dir = pj(TMPDIR, '.theano-{cohort}-{scatter}'),
            output_dir = pj(GATK_gCNV, '{cohort}_scatter_{scatter}')
    conda: CONDA_GATK_CNV
    resources:
            mem_mb = lambda wildcards, attempt: attempt* 45000,
            n = 17
    log: pj(LOG,"gCNV_gatk", '{cohort}.{scatter}.germlinecnvcalling.log')
    benchmark: pj(BENCH, '{cohort}.{scatter}.germlinecnvcalling.txt')
    shell:
        """
            mkdir -p {params.theano_complie_dir}
            export OMP_NUM_THREADS={resources.n} 
            THEANO_FLAGS="base_compiledir={params.theano_complie_dir}"  {params.java} -jar {params.gatk} GermlineCNVCaller {params.inputs} -L {input.scatters} --contig-ploidy-calls  {params.contig_ploydi_calls} --interval-merging-rule OVERLAPPING_ONLY --run-mode COHORT --output {params.output_dir} --output-prefix scatterd_{wildcards.cohort}_{wildcards.scatter} 2> {log}
        """

rule PostprocessGermlineCNVCalls:
    input:
        samples = input_func,
        models= expand(pj(GATK_gCNV,'{cohort}_scatter_{scatter}','scatterd_{cohort}_{scatter}-model','calling_config.json'),scatter = gcnv_scatter_shards, allow_missing = True),
         # sample_index = expand(pj(GATK_gCNV, '{cohort}_scatter_{scatter}', 'scatterd_{cohort}_{scatter}-calls', 'SAMPLE_{index}', 'sample_name.txt'), scatter = gcnv_scatter_shards, allow_missing = True),
    output:
        genotyped_intervals = pj(GATK_gCNV, 'GENOTYPED_CALLS_intervals_{cohort}', 'COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz'),
        genotyped_segments = pj(GATK_gCNV, 'GENOTYPED_CALLS_segments_{cohort}', 'COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz'),
        denoised_copy_ratio= pj(GATK_gCNV,'GENOTYPED_CALLS_denoised_copy_ratio_{cohort}','COHORT_{cohort}_SAMPLE_{sample}_{index}.tsv')
    conda: CONDA_GATK_CNV
    params:
        java= java_cnv,
        gatk= gatk_cnv,
        model_shards = expand(' --model-shard-path {gatk_gcnv}/{cohort}_scatter_{scatter}/scatterd_{cohort}_{scatter}-model ', gatk_gcnv=GATK_gCNV, scatter = gcnv_scatter_shards, allow_missing = True),
        calls_shards = expand(' --calls-shard-path {gatk_gcnv}/{cohort}_scatter_{scatter}/scatterd_{cohort}_{scatter}-calls ', gatk_gcnv=GATK_gCNV, scatter = gcnv_scatter_shards, allow_missing = True),
        CPC = (pj(GATK_gCNV,  '{cohort}-calls')),
        SD = REF_MALE_DICT,
        ref = REF_MALE
    benchmark: pj(BENCH, '{cohort}.{sample}.{index}.PostprocessGermlineCNVcalls.txt')
    log: pj(LOG,"gCNV_gatk", '{cohort}.{sample}.{index}.PostprocessGermlineCNVcalls')
    resources:
            mem_mb = lambda wildcards, attempt: (attempt - 1) * 0.5 * 5000 + 9000
    shell:
            """
            {params.java} -jar {params.gatk} PostprocessGermlineCNVCalls -R {params.ref} --sequence-dictionary {params.SD}  {params.model_shards} {params.calls_shards} --contig-ploidy-calls {params.CPC} --sample-index {wildcards.index} --allosomal-contig chrX --allosomal-contig chrY --output-genotyped-intervals {output.genotyped_intervals} --output-genotyped-segments {output.genotyped_segments} --output-denoised-copy-ratios {output.denoised_copy_ratio} 2> {log}
            """
