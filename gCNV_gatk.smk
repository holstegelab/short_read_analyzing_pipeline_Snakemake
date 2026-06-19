from common import *
import csv
import h5py
import numpy as np
import os

onsuccess: shell(f"rm -fr {pj(LOG, 'gCNV_gatk')}/*")
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
GCNV_INTERVALS_DIR = pj(GATK_gCNV, "intervals")
GCNV_WGS_INTERVALS = pj(GCNV_INTERVALS_DIR, f"wgs_{GCNV_WGS_BIN_LENGTH}bp.preprocessed.interval_list")
GCNV_WGS_ANNOTATED_INTERVALS = pj(GCNV_INTERVALS_DIR, f"wgs_{GCNV_WGS_BIN_LENGTH}bp.annotated.tsv")
GCNV_CLUSTER_DATASET = config.get("gcnv_cluster_dataset", "cluster_3_cor_cov")
GCNV_MISSING_CLUSTER_MODE = config.get("gcnv_missing_cluster_mode", "single")
GCNV_COHORT_PLOT_BATCH_SIZE = max(2, int(config.get("gcnv_cohort_plot_batch_size", 1000)))
GCNV_COHORT_PLOT_MAX_DELTA = float(config.get("gcnv_cohort_plot_max_delta", 1e-5))
GCNV_COHORT_REPORT_DIR = pj(GATK_gCNV, "cohort_report")
GCNV_COHORT_SAMPLE_LIST = pj(GCNV_COHORT_REPORT_DIR, f"{GCNV_CLUSTER_DATASET}.samples_per_cohort.tsv")
GCNV_UNASSIGNED_SAMPLE_LIST = pj(GCNV_COHORT_REPORT_DIR, f"{GCNV_CLUSTER_DATASET}.unassigned_samples.tsv")
GCNV_COHORT_PCA_COORDINATES = pj(GCNV_COHORT_REPORT_DIR, f"{GCNV_CLUSTER_DATASET}.pca_2d_coordinates.tsv")
GCNV_COHORT_DISTANCE_PLOT = pj(GCNV_COHORT_REPORT_DIR, f"{GCNV_CLUSTER_DATASET}.cohort_distance_2d.png")
GCNV_COHORT_REPORT_OUTPUTS = [
    GCNV_COHORT_SAMPLE_LIST,
    GCNV_UNASSIGNED_SAMPLE_LIST,
    GCNV_COHORT_PCA_COORDINATES,
    GCNV_COHORT_DISTANCE_PLOT,
]


def samplefile_stat_path(samplefile, suffix):
    return pj(get_samplefile_folder(samplefile), f"{samplefile}.{suffix}")

def samples_for_samplefile(samplefile):
    return sorted(SAMPLEFILE_TO_SAMPLES[samplefile].keys())


def decode_hdf5_string(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def missing_cluster_cohort(samplefile_index):
    if GCNV_MISSING_CLUSTER_MODE == "single":
        return 0
    if GCNV_MISSING_CLUSTER_MODE == "samplefile":
        return samplefile_index
    if GCNV_MISSING_CLUSTER_MODE == "error":
        raise ValueError(
            f"Missing {GCNV_CLUSTER_DATASET}; run PCA.smk first or set "
            "gcnv_missing_cluster_mode to 'single' or 'samplefile'."
        )
    raise ValueError(
        "Unsupported gcnv_missing_cluster_mode "
        f"{GCNV_MISSING_CLUSTER_MODE!r}; use 'single', 'samplefile', or 'error'."
    )


def add_fallback_cohort(cohort_to_samples, samplefile, samplefile_index, reason):
    cohort = missing_cluster_cohort(samplefile_index)
    samples = samples_for_samplefile(samplefile)
    cohort_to_samples.setdefault(cohort, []).extend(samples)
    print(
        f"gCNV: {reason}; assigning {samplefile} samples to cohort {cohort} "
        f"with gcnv_missing_cluster_mode={GCNV_MISSING_CLUSTER_MODE!r}."
    )
    return cohort, samples


def get_gcnv_cohort_assignments():
    cohort_to_samples = {}
    sample_to_cohort = {}
    unassigned_samples = {}
    for samplefile_index, samplefile in enumerate(SAMPLE_FILES):
        selected_samples = samples_for_samplefile(samplefile)
        selected_sample_set = set(selected_samples)
        hdf5_file = samplefile_stat_path(samplefile, "coverage.hdf5")
        if not os.path.exists(hdf5_file):
            cohort, samples = add_fallback_cohort(
                cohort_to_samples,
                samplefile,
                samplefile_index,
                f"{hdf5_file} does not exist",
            )
            sample_to_cohort.update({sample: cohort for sample in samples})
            continue

        with h5py.File(hdf5_file, "r") as f:
            data_group = f["coverage"]
            if GCNV_CLUSTER_DATASET not in data_group:
                cohort, samples = add_fallback_cohort(
                    cohort_to_samples,
                    samplefile,
                    samplefile_index,
                    f"{GCNV_CLUSTER_DATASET} is not present in {hdf5_file}",
                )
                sample_to_cohort.update({sample: cohort for sample in samples})
                continue

            hdf5_samples = [decode_hdf5_string(sample) for sample in data_group["samples"][()]]
            labels = np.asarray(data_group[GCNV_CLUSTER_DATASET][()])

        if len(hdf5_samples) != len(labels):
            raise ValueError(
                f"{hdf5_file}:{GCNV_CLUSTER_DATASET} has {len(labels)} labels "
                f"for {len(hdf5_samples)} samples."
            )

        seen_samples = set()
        for sample, label in zip(hdf5_samples, labels):
            if sample not in selected_sample_set:
                continue
            seen_samples.add(sample)
            label = int(label)
            if label == -1:
                unassigned_samples[sample] = (samplefile, f"{GCNV_CLUSTER_DATASET}=-1")
                continue
            cohort_to_samples.setdefault(label, []).append(sample)
            sample_to_cohort[sample] = label

        for sample in sorted(selected_sample_set - seen_samples):
            unassigned_samples[sample] = (samplefile, "missing_from_coverage_hdf5")

    for sample in sample_to_cohort:
        unassigned_samples.pop(sample, None)

    if not cohort_to_samples:
        raise ValueError(
            "No samples were assigned to gCNV cohorts. Check the selected "
            f"{GCNV_CLUSTER_DATASET} labels or use gcnv_missing_cluster_mode=single."
        )

    cohort_to_samples = {cohort: sorted(set(samples)) for cohort, samples in cohort_to_samples.items()}
    return cohort_to_samples, sample_to_cohort, unassigned_samples


def get_cohort_to_samples():
    cohort_to_samples, _, _ = get_gcnv_cohort_assignments()
    return cohort_to_samples


def sample_to_samplefile_map():
    return {
        sample: samplefile
        for samplefile in SAMPLE_FILES
        for sample in samples_for_samplefile(samplefile)
    }


cohort_to_samples = get_cohort_to_samples()
groups = set(cohort_to_samples)


def get_samples_in_group(cohort):
    samples = cohort_to_samples[int(cohort)]
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
        paths_output + GCNV_COHORT_REPORT_OUTPUTS
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


rule write_gcnv_cohort_report:
    input:
        calls = paths_output
    output:
        samples = GCNV_COHORT_SAMPLE_LIST,
        unassigned = GCNV_UNASSIGNED_SAMPLE_LIST,
        coordinates = GCNV_COHORT_PCA_COORDINATES,
        plot = GCNV_COHORT_DISTANCE_PLOT
    conda: CONDA_PCA
    resources:
        n = 1,
        mem_mb = 32000
    run:
        import pandas as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sklearn.decomposition import IncrementalPCA

        os.makedirs(GCNV_COHORT_REPORT_DIR, exist_ok=True)

        cohort_to_samples, sample_to_cohort, unassigned_samples = get_gcnv_cohort_assignments()
        sample_to_samplefile = sample_to_samplefile_map()

        with open(output.samples, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["cohort", "samplefile", "sample"])
            for cohort in sorted(cohort_to_samples):
                for sample in cohort_to_samples[cohort]:
                    writer.writerow([cohort, sample_to_samplefile.get(sample, ""), sample])

        with open(output.unassigned, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["samplefile", "sample", "reason"])
            for sample, (samplefile, reason) in sorted(
                unassigned_samples.items(),
                key=lambda item: (item[1][0], item[0]),
            ):
                writer.writerow([samplefile, sample, reason])

        def compute_pca_coordinates():
            selected_samples = set(sample_to_samplefile)
            hdf5_files = [
                samplefile_stat_path(samplefile, "coverage.hdf5")
                for samplefile in SAMPLE_FILES
                if os.path.exists(samplefile_stat_path(samplefile, "coverage.hdf5"))
            ]
            if len(selected_samples) < 2 or not hdf5_files:
                return {}

            total_intervals = 0
            for hdf5_file in hdf5_files:
                with h5py.File(hdf5_file, "r") as f:
                    total_intervals = max(total_intervals, f["coverage"]["coverage"].shape[1])
            if total_intervals < 2:
                return {}

            main_chr_prefixes = tuple([f"chr{i}:" for i in range(1, 23)] + ["chrX:", "chrY:"])
            ipca = IncrementalPCA(n_components=2)
            previous_components = None
            max_delta = np.inf
            sample_names = None
            fitted = False
            start_interval = 0

            while start_interval < total_intervals and max_delta > GCNV_COHORT_PLOT_MAX_DELTA:
                batch_df = pd.DataFrame()
                end_interval = start_interval + GCNV_COHORT_PLOT_BATCH_SIZE

                for hdf5_file in hdf5_files:
                    with h5py.File(hdf5_file, "r") as f:
                        data_group = f["coverage"]
                        coverage = data_group["coverage"]
                        if start_interval >= coverage.shape[1]:
                            continue
                        hdf5_end = min(end_interval, coverage.shape[1])
                        samples = [decode_hdf5_string(sample) for sample in data_group["samples"][()]]
                        bases_mapped_raw = np.asarray(data_group["bases_mapped"][:], dtype=np.float64)
                        selected_mask = np.asarray([sample in selected_samples for sample in samples], dtype=bool)
                        valid_mask = selected_mask & np.isfinite(bases_mapped_raw) & (bases_mapped_raw > 0)
                        if not np.any(valid_mask):
                            continue

                        valid_samples = [
                            sample
                            for sample, valid in zip(samples, valid_mask)
                            if valid
                        ]
                        chroms = [decode_hdf5_string(chrom) for chrom in data_group["chrom"][start_interval:hdf5_end]]
                        starts = data_group["start"][start_interval:hdf5_end]
                        ends = data_group["end"][start_interval:hdf5_end]
                        interval_index = [
                            f"{chrom}:{start}-{end}"
                            for chrom, start, end in zip(chroms, starts, ends)
                        ]
                        coverage_values = np.asarray(
                            coverage[valid_mask, start_interval:hdf5_end],
                            dtype=np.float64,
                        )
                        bases_mapped = bases_mapped_raw[valid_mask][:, np.newaxis] / 1000000000.0
                        corrected_coverage = coverage_values / bases_mapped
                        hdf5_df = pd.DataFrame(
                            columns=valid_samples,
                            data=corrected_coverage.T,
                            index=interval_index,
                        )
                        hdf5_df = hdf5_df[hdf5_df.index.str.startswith(main_chr_prefixes)]
                        batch_df = pd.concat([batch_df, hdf5_df], axis=1)

                batch_df = batch_df.replace([np.inf, -np.inf], np.nan)
                if sample_names is None:
                    batch_df = batch_df.dropna(axis=0, how="any")
                    if batch_df.empty:
                        start_interval = end_interval
                        continue
                    sample_names = list(batch_df.columns)
                    if len(sample_names) < 2:
                        return {}
                else:
                    batch_df = batch_df.reindex(columns=sample_names)
                    batch_df = batch_df.dropna(axis=0, how="any")
                    if batch_df.empty:
                        start_interval = end_interval
                        continue

                if batch_df.shape[0] < 2:
                    start_interval = end_interval
                    continue

                ipca.partial_fit(batch_df.to_numpy())
                fitted = True
                components = ipca.components_.copy()
                if previous_components is not None:
                    max_delta = np.max(np.abs(components - previous_components))
                previous_components = components
                start_interval = end_interval

            if not fitted:
                return {}

            return {
                sample: (float(pc1), float(pc2))
                for sample, (pc1, pc2) in zip(sample_names, ipca.components_.T)
            }

        pca_coordinates = compute_pca_coordinates()

        with open(output.coordinates, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["samplefile", "sample", "assignment", "cohort", "pc1", "pc2"])
            for sample in sorted(sample_to_samplefile, key=lambda s: (sample_to_samplefile[s], s)):
                cohort = sample_to_cohort.get(sample, "")
                assignment = "assigned" if sample in sample_to_cohort else "unassigned"
                pc1, pc2 = pca_coordinates.get(sample, ("", ""))
                writer.writerow([sample_to_samplefile[sample], sample, assignment, cohort, pc1, pc2])

        fig, ax = plt.subplots(figsize=(9, 6.5))
        if len(pca_coordinates) < 2:
            ax.text(
                0.5,
                0.5,
                "2D PCA unavailable\ncoverage HDF5 data had fewer than two usable samples",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_axis_off()
        else:
            color_map = plt.get_cmap("tab20")
            cohorts = sorted(cohort_to_samples)
            cohort_to_color = {
                cohort: color_map(index % color_map.N)
                for index, cohort in enumerate(cohorts)
            }

            for cohort in cohorts:
                cohort_samples = [
                    sample
                    for sample in cohort_to_samples[cohort]
                    if sample in pca_coordinates
                ]
                if not cohort_samples:
                    continue
                xs = [pca_coordinates[sample][0] for sample in cohort_samples]
                ys = [pca_coordinates[sample][1] for sample in cohort_samples]
                ax.scatter(
                    xs,
                    ys,
                    s=24,
                    alpha=0.72,
                    color=cohort_to_color[cohort],
                    label=f"Cohort {cohort} (n={len(cohort_samples)})",
                )
                centroid_x = float(np.mean(xs))
                centroid_y = float(np.mean(ys))
                ax.scatter(
                    [centroid_x],
                    [centroid_y],
                    s=150,
                    marker="X",
                    color=cohort_to_color[cohort],
                    edgecolor="black",
                    linewidth=0.8,
                )
                ax.text(
                    centroid_x,
                    centroid_y,
                    f" {cohort}",
                    fontsize=9,
                    fontweight="bold",
                    va="center",
                )

            unassigned_with_coordinates = [
                sample
                for sample in unassigned_samples
                if sample in pca_coordinates
            ]
            if unassigned_with_coordinates:
                ax.scatter(
                    [pca_coordinates[sample][0] for sample in unassigned_with_coordinates],
                    [pca_coordinates[sample][1] for sample in unassigned_with_coordinates],
                    s=32,
                    marker="x",
                    color="#666666",
                    label=f"Unassigned (n={len(unassigned_with_coordinates)})",
                )

            ax.axhline(0, color="#d0d0d0", linewidth=0.6)
            ax.axvline(0, color="#d0d0d0", linewidth=0.6)
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_title(f"gCNV cohort distances: {GCNV_CLUSTER_DATASET}")
            ax.legend(loc="best", fontsize=8, frameon=False)

        fig.tight_layout()
        fig.savefig(output.plot, dpi=200)
        plt.close(fig)
