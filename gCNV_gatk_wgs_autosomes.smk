from common import *
import csv
import os

GCNV_AUTOSOMES_DIR = pj(GATK_gCNV, "wgs_autosomes_all_samples")

onsuccess: shell(f"if [ -d {quote(pj(LOG, 'gCNV_gatk_wgs_autosomes'))} ]; then mkdir -p {quote(pj(GCNV_AUTOSOMES_DIR, 'debug_archives'))}; tar -C {quote(LOG)} -czf {quote(pj(GCNV_AUTOSOMES_DIR, 'debug_archives', 'gcnv_wgs_autosomes_logs.tar.gz'))} gCNV_gatk_wgs_autosomes; find {quote(pj(LOG, 'gCNV_gatk_wgs_autosomes'))} -type f -delete; find {quote(pj(LOG, 'gCNV_gatk_wgs_autosomes'))} -type d -empty -delete; fi")
wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    cohort=r"\d+",
    index=r"\d+",
    scatter=r"\d{4}_of_\d+",

module DownloadAndExtractBam:
    snakefile: 'Download_and_extract_bam.smk'
    config: config
use rule download_encrypted_cram_from_tape from DownloadAndExtractBam
use rule decrypt_cram from DownloadAndExtractBam
use rule cram_to_markdup_bam from DownloadAndExtractBam

GCNV_WGS_BIN_LENGTH = int(config.get("gcnv_wgs_bin_length", 1000))
GCNV_WGS_SCATTER_COUNT = int(config.get("gcnv_wgs_scatter_count", 300))
GCNV_WGS_INTERVAL_SOURCE = AUTO_ONLY_BED
GCNV_INTERVALS_DIR = pj(GCNV_AUTOSOMES_DIR, "intervals")
GCNV_WGS_INTERVALS = pj(GCNV_INTERVALS_DIR, f"wgs_{GCNV_WGS_BIN_LENGTH}bp.preprocessed.interval_list")
GCNV_WGS_ANNOTATED_INTERVALS = pj(GCNV_INTERVALS_DIR, f"wgs_{GCNV_WGS_BIN_LENGTH}bp.annotated.tsv")
GCNV_COHORT_REPORT_DIR = pj(GCNV_AUTOSOMES_DIR, "cohort_report")
GCNV_COHORT_SAMPLE_LIST = pj(GCNV_COHORT_REPORT_DIR, "all_samples_autosomes.samples_per_cohort.tsv")
GCNV_COHORT_REPORT_OUTPUTS = [
    GCNV_COHORT_SAMPLE_LIST,
]
GCNV_AUTOSOMES_COHORT = 0

def samples_for_samplefile(samplefile):
    return sorted(SAMPLEFILE_TO_SAMPLES[samplefile].keys())


def get_gcnv_cohort_assignments():
    samples = sorted(SAMPLEINFO.keys())
    if not samples:
        raise ValueError("No samples were selected for the autosomal gCNV workflow.")
    cohort_to_samples = {GCNV_AUTOSOMES_COHORT: samples}
    sample_to_cohort = {sample: GCNV_AUTOSOMES_COHORT for sample in samples}
    return cohort_to_samples, sample_to_cohort


def sample_to_samplefile_map():
    return {
        sample: samplefile
        for samplefile in SAMPLE_FILES
        for sample in samples_for_samplefile(samplefile)
    }


cohort_to_samples, sample_to_cohort = get_gcnv_cohort_assignments()
groups = set(cohort_to_samples)


def get_samples_in_group(cohort):
    samples = cohort_to_samples[int(cohort)]
    idx = list(range(len(samples)))
    return samples, idx

paths_output = []
for cohort in groups:
    sample_names, idxs = get_samples_in_group(cohort=cohort)
    for sample_name, idx in zip(sample_names, idxs):
        path = pj(GCNV_AUTOSOMES_DIR, f'GENOTYPED_CALLS_intervals_{cohort}', f'COHORT_{cohort}_SAMPLE_{sample_name}_{idx}.vcf.gz')
        paths_output.append(path)



def input_func(wildcards):
    cohort = wildcards['cohort']
    samples_in_group, _ = get_samples_in_group(cohort=cohort)
    return expand(
        '{gatk_gcnv}/Read_counts_hdf5/{cohort}/{sample}_readcounts.hdf5',
        gatk_gcnv=GCNV_AUTOSOMES_DIR,
        cohort=cohort,
        sample=samples_in_group,
        allow_missing=True
    )


def sample_list_per_cohort(wildcards):
    cohort = wildcards['cohort']
    samples_in_group, _ = get_samples_in_group(cohort)
    return expand(
        ' -I {gatk_gcnv}/Read_counts_hdf5/{cohort}/{sample}_readcounts.hdf5 ',
        gatk_gcnv=GCNV_AUTOSOMES_DIR,
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
gcnv_debug_archives = [
    pj(GCNV_AUTOSOMES_DIR, 'debug_archives', f'cohort_{cohort}.gcnv_debug.tar.gz')
    for cohort in groups
]

def gcnv_shard_dir(cohort, scatter):
    return pj(GCNV_AUTOSOMES_DIR, f'{cohort}_scatter_{scatter}')

def gcnv_model_shard_args(wildcards):
    return " ".join(
        f"--model-shard-path {quote(pj(gcnv_shard_dir(wildcards.cohort, scatter), f'scatterd_{wildcards.cohort}_{scatter}-model'))}"
        for scatter in gcnv_scatter_shards
    )

def gcnv_calls_shard_args(wildcards):
    return " ".join(
        f"--calls-shard-path {quote(pj(gcnv_shard_dir(wildcards.cohort, scatter), f'scatterd_{wildcards.cohort}_{scatter}-calls'))}"
        for scatter in gcnv_scatter_shards
    )

def gcnv_debug_archive_items(wildcards):
    ploidy_dirs = [
        quote(f"{wildcards.cohort}-model"),
        quote(f"{wildcards.cohort}-calls"),
    ]
    shard_dirs = [
        quote(f"{wildcards.cohort}_scatter_{scatter}")
        for scatter in gcnv_scatter_shards
    ]
    return " ".join(ploidy_dirs + shard_dirs)

rule gCNV_gatk_wgs_autosomes_all:
    input:
        paths_output + GCNV_COHORT_REPORT_OUTPUTS + gcnv_debug_archives
    default_target: True

rule preprocess_wgs_gcnv_intervals:
    input:
        intervals = GCNV_WGS_INTERVAL_SOURCE
    output:
        intervals = temp(GCNV_WGS_INTERVALS)
    params:
        java = java_cnv,
        gatk = gatk_cnv,
        ref = REF_MALE,
        bin_length = GCNV_WGS_BIN_LENGTH,
        output_dir = GCNV_INTERVALS_DIR
    conda: CONDA_GATK_CNV
    log: pj(LOG, "gCNV_gatk_wgs_autosomes", f"preprocess_wgs_{GCNV_WGS_BIN_LENGTH}bp_intervals.log")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.java} -jar {params.gatk} PreprocessIntervals -R {params.ref} -L {input.intervals} -imr OVERLAPPING_ONLY --bin-length {params.bin_length} -O {output.intervals} 2> {log}
        """

rule collect_read_counts:
    input: bam=markdup_bam_input,
            bai=markdup_bai_input,
            intervals = rules.preprocess_wgs_gcnv_intervals.output.intervals
    output: ReadCounts = ensure(pj(GCNV_AUTOSOMES_DIR, 'Read_counts_hdf5', '{cohort}', '{sample}_readcounts.hdf5'), non_empty=True)
    params: java = java_cnv,
            gatk = gatk_cnv,
            ref = REF_MALE
    conda: CONDA_GATK_CNV
    resources:
            n = 1,
            mem_mb = 8500 # if tool doesn't have enough RAM it ends without error and create corrupted file. Therefore here pretty big RAM (mean usage much lower)
    log: pj(LOG,"gCNV_gatk_wgs_autosomes", '{sample}.{cohort}.collectreadcounts_gatkcnv.log')
    benchmark: pj(BENCH, '{sample}.{cohort}.collectreadcounts_wgs_autosomes_gatkcnv.txt')
    shell: """
            {params.java} -jar {params.gatk} CollectReadCounts -L {input.intervals} -R {params.ref} -imr OVERLAPPING_ONLY -I {input.bam} -O {output.ReadCounts} 2> {log}
            """

rule annotateintervals:
    input: interval = rules.preprocess_wgs_gcnv_intervals.output.intervals
    output: annotated_tsv = temp(GCNV_WGS_ANNOTATED_INTERVALS)
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
    output: filtered_intervals = temp(pj(GCNV_AUTOSOMES_DIR, 'filtered_intervals', '{cohort}_filtered.interval_list'))
    params: inputs = sample_list_per_cohort,
            java= java_cnv,
            gatk=gatk_cnv,
            PAR_and_CENTROMERIC = PAR_and_CENTROMERIC,
            cut_off_samples = '80',
            cut_off_count_threshold = '10',
    conda: CONDA_GATK_CNV
    log: pj(LOG,"gCNV_gatk_wgs_autosomes", '{cohort}.filterintervals_gatkcnv.log')
    shell: """
            {params.java} -jar {params.gatk} FilterIntervals --annotated-intervals {input.annotation} -L {input.intervals} {params.inputs} -imr OVERLAPPING_ONLY -O {output.filtered_intervals} -XL {params.PAR_and_CENTROMERIC} --low-count-filter-percentage-of-samples {params.cut_off_samples} --low-count-filter-count-threshold {params.cut_off_count_threshold}  2> {log}
            """

rule make_scatters:
    input: rules.filterintervals.output.filtered_intervals
    output: scatter_dir = temp(directory(pj(GCNV_AUTOSOMES_DIR, 'filtered_intervals', 'cohort-{cohort}')))
    conda: CONDA_GATK_CNV
    log: pj(LOG,"gCNV_gatk_wgs_autosomes",'{cohort}.make_scatters.log')
    params:
            java = java_cnv,
            gatk = gatk_cnv,
            scatter_count = GCNV_WGS_SCATTER_COUNT,
            output_dir = pj(GCNV_AUTOSOMES_DIR, 'filtered_intervals', 'cohort-{cohort}')
    shell:
        """
        mkdir -p {params.output_dir}
        {params.java} -jar {params.gatk} IntervalListTools -I {input} --OUTPUT {params.output_dir} --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_COUNT {params.scatter_count}
        """


rule DetermineGCP:
    input: samples = input_func,
            intervals = pj(GCNV_AUTOSOMES_DIR, 'filtered_intervals', '{cohort}_filtered.interval_list')
    output:
            model_dir = temp(directory(pj(GCNV_AUTOSOMES_DIR,  '{cohort}-model'))),
            calls_dir = temp(directory(pj(GCNV_AUTOSOMES_DIR,  '{cohort}-calls'))),
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            contig_ploydi_priors = PL_PR_TABLE,
            output_dir = GCNV_AUTOSOMES_DIR
    conda: CONDA_GATK_CNV
    resources:
            n = 1
    shell: """
            export OMP_NUM_THREADS={resources.n}
            {params.java} -jar {params.gatk} DetermineGermlineContigPloidy  --output-prefix {wildcards.cohort} --output {params.output_dir} {params.inputs} -L {input.intervals} -imr OVERLAPPING_ONLY --contig-ploidy-priors {params.contig_ploydi_priors}
            """

rule GermlineCNVCaller:
    input:  samples = input_func,
            scatter_dir = rules.make_scatters.output.scatter_dir,
            contig_ploudi_calls = rules.DetermineGCP.output.calls_dir,
    output: shard_dir = temp(directory(pj(GCNV_AUTOSOMES_DIR,'{cohort}_scatter_{scatter}'))),
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            scatter_interval = pj(GCNV_AUTOSOMES_DIR, 'filtered_intervals', 'cohort-{cohort}', 'temp_{scatter}', 'scattered.interval_list'),
            theano_complie_dir = pj(TMPDIR, '.theano-{cohort}-{scatter}'),
            output_dir = pj(GCNV_AUTOSOMES_DIR, '{cohort}_scatter_{scatter}')
    conda: CONDA_GATK_CNV
    resources:
            mem_mb = lambda wildcards, attempt: attempt* 45000,
            n = 17
    log: pj(LOG,"gCNV_gatk_wgs_autosomes", '{cohort}.{scatter}.germlinecnvcalling.log')
    benchmark: pj(BENCH, '{cohort}.{scatter}.wgs_autosomes.germlinecnvcalling.txt')
    shell:
        """
            mkdir -p {params.theano_complie_dir}
            export OMP_NUM_THREADS={resources.n} 
            THEANO_FLAGS="base_compiledir={params.theano_complie_dir}"  {params.java} -jar {params.gatk} GermlineCNVCaller {params.inputs} -L {params.scatter_interval} --contig-ploidy-calls  {input.contig_ploudi_calls} --interval-merging-rule OVERLAPPING_ONLY --run-mode COHORT --output {params.output_dir} --output-prefix scatterd_{wildcards.cohort}_{wildcards.scatter} 2> {log}
            rm -rf {params.theano_complie_dir}
        """

rule PostprocessGermlineCNVCalls:
    input:
        samples = input_func,
        shards= expand(pj(GCNV_AUTOSOMES_DIR,'{cohort}_scatter_{scatter}'), scatter = gcnv_scatter_shards, allow_missing = True),
        cpc = rules.DetermineGCP.output.calls_dir,
    output:
        genotyped_intervals = pj(GCNV_AUTOSOMES_DIR, 'GENOTYPED_CALLS_intervals_{cohort}', 'COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz'),
        genotyped_segments = pj(GCNV_AUTOSOMES_DIR, 'GENOTYPED_CALLS_segments_{cohort}', 'COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz'),
        denoised_copy_ratio= pj(GCNV_AUTOSOMES_DIR,'GENOTYPED_CALLS_denoised_copy_ratio_{cohort}','COHORT_{cohort}_SAMPLE_{sample}_{index}.tsv')
    conda: CONDA_GATK_CNV
    params:
        java= java_cnv,
        gatk= gatk_cnv,
        model_shards = gcnv_model_shard_args,
        calls_shards = gcnv_calls_shard_args,
        SD = REF_MALE_DICT,
        ref = REF_MALE
    benchmark: pj(BENCH, '{cohort}.{sample}.{index}.wgs_autosomes.PostprocessGermlineCNVcalls.txt')
    log: pj(LOG,"gCNV_gatk_wgs_autosomes", '{cohort}.{sample}.{index}.PostprocessGermlineCNVcalls')
    resources:
            mem_mb = lambda wildcards, attempt: (attempt - 1) * 0.5 * 5000 + 9000
    shell:
            """
            {params.java} -jar {params.gatk} PostprocessGermlineCNVCalls -R {params.ref} --sequence-dictionary {params.SD}  {params.model_shards} {params.calls_shards} --contig-ploidy-calls {input.cpc} --sample-index {wildcards.index} --output-genotyped-intervals {output.genotyped_intervals} --output-genotyped-segments {output.genotyped_segments} --output-denoised-copy-ratios {output.denoised_copy_ratio} 2> {log}
            """

rule archive_gcnv_debug:
    input:
        ploidy_model = rules.DetermineGCP.output.model_dir,
        ploidy_calls = rules.DetermineGCP.output.calls_dir,
        shards = expand(pj(GCNV_AUTOSOMES_DIR,'{cohort}_scatter_{scatter}'), scatter = gcnv_scatter_shards, allow_missing = True)
    output:
        archive = pj(GCNV_AUTOSOMES_DIR, 'debug_archives', 'cohort_{cohort}.gcnv_debug.tar.gz')
    params:
        gatk_gcnv = GCNV_AUTOSOMES_DIR,
        archive_dir = pj(GCNV_AUTOSOMES_DIR, 'debug_archives'),
        archive_items = gcnv_debug_archive_items
    shell:
        """
        mkdir -p {params.archive_dir}
        tar -C {params.gatk_gcnv} -czf {output.archive} {params.archive_items}
        """


rule write_gcnv_cohort_report:
    input:
        calls = paths_output
    output:
        samples = GCNV_COHORT_SAMPLE_LIST
    resources:
        n = 1,
        mem_mb = 100
    run:
        os.makedirs(GCNV_COHORT_REPORT_DIR, exist_ok=True)

        sample_to_samplefile = sample_to_samplefile_map()

        with open(output.samples, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["cohort", "samplefile", "sample"])
            for cohort in sorted(cohort_to_samples):
                for sample in cohort_to_samples[cohort]:
                    writer.writerow([cohort, sample_to_samplefile.get(sample, ""), sample])
