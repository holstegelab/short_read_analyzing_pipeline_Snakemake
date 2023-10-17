from common import *
import h5py
import numpy as np

wildcard_constraints:
    sample="[\w\d_\-@]+",
    cohorts = "[\d]",
    index = "(\d+)"
    # readgroup="[\w\d_\-@]+"

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

module Stat:
    snakefile: 'Stat.smk'
    config: config
use rule * from Stat


module Tools:
    snakefile: 'Tools.smk'
    config: config
# use rule * from Tools

def get_sample_name(sample_index, cohort):
    # Get the sample name from the "sample_name.txt" file
    sample_name_file = os.path.join(GATK_gCNV, f'{cohort}_scatter_0001_of_148', f'scatterd_{cohort}_0001_of_148-calls', f'SAMPLE_{sample_index}', 'sample_name.txt')
    with open(sample_name_file, 'r') as f:
        sample_name = f.read().strip()
    return sample_name


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




def get_n_samples_in_group(cohort, hdf5_files = hdf5_files):
    samples_per_cohort = 0
    for hdf5_file in hdf5_files:
        with h5py.File(hdf5_file, 'r') as f:
            data_group = f['coverage']
            sample_list_per_hdf5_file = [s.decode('utf-8') for s in data_group['samples'][()]]
            cohort_mask = np.equal(data_group['cluster_3_cor_cov'], int(cohort))
            cohort_samples = [sample_list_per_hdf5_file[i] for i, mask in enumerate(cohort_mask) if mask]
            samples_per_cohort += len(cohort_samples)
    return samples_per_cohort


paths_output = []
for cohort in groups:
    sample_names, idxs = get_samples_in_group(cohort=cohort)
    for sample_name, idx in zip(sample_names, idxs):
        path =  f'GATK_gCNV/GENOTYPED_CALLS_intervals_{cohort}/COHORT_{cohort}_SAMPLE_{sample_name}_{idx}.vcf.gz'
        paths_output.append(path)



def input_func(wildcards):
    cohort = wildcards['cohort']
    samples_in_group, index = get_samples_in_group(cohort = cohort)
    return expand(
        '{gatk_gcnv}/Read_counts_hdf5/{cohort}/{sample}_readcounts.hdf5',
        gatk_gcnv=GATK_gCNV,
        cohort=cohort,
        sample=samples_in_group,
        index = index,
        allow_missing=True
    )

# def input_bams(wildcards):
#     cohort = wildcards['cohort']
#     samples_in_group, index = get_samples_in_group(cohort = cohort)
#     return expand(
#         '{bam}/{sample}.markdup.bam',
#         bam=BAM,
#         cohort=cohort,
#         sample=samples_in_group,
#         allow_missing=True
#     )

def sample_list_per_cohort(wildcards):
    cohort = wildcards['cohort']
    samples_in_group, index = get_samples_in_group(cohort)
    return expand(
        ' -I {gatk_gcnv}/Read_counts_hdf5/{cohort}/{sample}_readcounts.hdf5 ',
        gatk_gcnv=GATK_gCNV,
        cohort=cohort,
        sample=samples_in_group,
        index = index,
        allow_missing=True
    )

def generate_scatter_list(start, end):
    string_list = [f"{i:04d}_of_{end}" for i in range(int(start), int(end) + 1)]
    return string_list

scatter_merged_cature_kit = generate_scatter_list('0001', '25')

rule gCNV_gatk_all:
    input:
        # expand('GATK_gCNV/{cohort}_scatter_{scatter}/scatterd_{cohort}_{scatter}-model/calling_config.json',scatter=scatter_merged_cature_kit,cohort=groups),
        # expand('GATK_gCNV/filtred_intervals/{cohort}_filtred.interval_list', cohort=groups),
        # expand('GATK_gCNV/{cohort}-model/contig_ploidy_prior.tsv', cohort = groups),
        paths_output
        # expand('GATK_gCNV/GENOTYPED_CALLS_intervals_{cohort}/COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz')
        # rules.Stat_all.input
    default_target: True

rule collect_read_counts:
    input: bam=rules.markdup.output.mdbams,
            # Merged capture kit for test
            capture_kit = MERGED_CAPTURE_KIT_IVL_CNV
    output: ReadCounts = ensure(pj(GATK_gCNV, 'Read_counts_hdf5', '{cohort}', '{sample}_readcounts.hdf5'), non_empty=True)
    params: java = java_cnv,
            gatk = gatk_cnv,
            ref = REF_MALE
    conda: CONDA_GATK_CNV
    resources:
            n = 1,
            mem_mb = 5500 # if tool doesn't have enough RAM it ends without error and create corrupted file. Therefore here pretty big RAM (mean usage much lower)
    log: pj(LOG, '{sample}.{cohort}.collectreadcounts_gatkcnv.log')
    benchmark: pj(BENCH, '{sample}.{cohort}.collectreadcounts_gatkcnv.txt')
    shell: """
            {params.java} -jar {params.gatk} CollectReadCounts -L {input.capture_kit} -R {params.ref} -imr OVERLAPPING_ONLY -I {input.bam} -O {output.ReadCounts} 2> {log}
            """

rule annotateintervals:
    input: interval = MERGED_CAPTURE_KIT_IVL_CNV
    output: annotated_tsv = pj(INTERVALS_DIR, 'preprocessed_intervals_for_GATK_CNV', 'merged_capture_kits_cds.annotated.tsv')
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
    output: filtered_intervals = pj(GATK_gCNV, 'filtred_intervals', '{cohort}_filtred.interval_list')
    params: inputs = sample_list_per_cohort,
            java= java_cnv,
            gatk=gatk_cnv,
            capture_kit = MERGED_CAPTURE_KIT_IVL_CNV,
            PAR_and_CENTROMERIC = PAR_and_CENTROMERIC,
            cut_off_samples = '80',
            cut_off_count_threshold = '10',
    conda: CONDA_GATK_CNV
    log: pj(LOG, '{cohort}.filterintervals_gatkcnv.log')
    benchmark: pj(BENCH, '{cohort}.filterintervals_gatkcnv.txt')
    shell: """
            {params.java} -jar {params.gatk} FilterIntervals --annotated-intervals {input.annotation} -L {params.capture_kit} {params.inputs} -imr OVERLAPPING_ONLY -O {output.filtered_intervals} -XL {params.PAR_and_CENTROMERIC} --low-count-filter-percentage-of-samples {params.cut_off_samples} --low-count-filter-count-threshold {params.cut_off_count_threshold}  2> {log}
            """

rule make_scatters:
    input: rules.filterintervals.output.filtered_intervals
    output: scatterd = expand(pj(GATK_gCNV, 'filtred_intervals', 'cohort-{cohort}', 'temp_{scatter}', 'scattered.interval_list'), scatter = scatter_merged_cature_kit, allow_missing = True)
    conda: CONDA_GATK_CNV
    log: pj(LOG,'{cohort}.make_scatters.log')
    benchmark: pj(BENCH,'{cohort}.make_scatters.txt')
    params:
            java = java_cnv,
            gatk = gatk_cnv,
            scatter_count = 25
    shell:
        """
        mkdir -p GATK_gCNV/filtred_intervals/cohort-{wildcards.cohort}
        {params.java} -jar {params.gatk} IntervalListTools -I {input} --OUTPUT GATK_gCNV/filtred_intervals/cohort-{wildcards.cohort} --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_COUNT {params.scatter_count}
        """


rule DetermineGCP:
    input: samples = input_func,
            intervals = pj(GATK_gCNV, 'filtred_intervals', '{cohort}_filtred.interval_list')
    output: CPC = pj(GATK_gCNV,  '{cohort}-model', 'contig_ploidy_prior.tsv'),
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            contig_ploydi_priors = PL_PR_TABLE
    conda: CONDA_GATK_CNV
    # log: pj(LOG, '{cohort}.determinecontigploydi.log')
    benchmark: pj(BENCH, '{cohort}.determinecontigploydi.txt')
    resources:
            n = 1
    shell: """
            export OMP_NUM_THREADS={resources.n}
            {params.java} -jar {params.gatk} DetermineGermlineContigPloidy  --output-prefix {wildcards.cohort} --output GATK_gCNV/ {params.inputs} -L {input.intervals} -imr OVERLAPPING_ONLY --contig-ploidy-priors {params.contig_ploydi_priors}
            """


rule GermlineCNVCaller:
    input:  samples = input_func,
            scatters = pj(GATK_gCNV, 'filtred_intervals', 'cohort-{cohort}', 'temp_{scatter}', 'scattered.interval_list'),
            contig_ploudi_calls = (pj(GATK_gCNV,  '{cohort}-model', 'contig_ploidy_prior.tsv')),
    output: models = pj(GATK_gCNV,'{cohort}_scatter_{scatter}','scatterd_{cohort}_{scatter}-model', 'calling_config.json'),
            # sample_files = pj(GATK_gCNV, '{cohort}_scatter_{scatter}', 'scatterd_{cohort}_{scatter}-calls', 'SAMPLE_{index}', 'sample_name.txt')
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            contig_ploydi_calls = (pj(GATK_gCNV,  '{cohort}-calls')),
            theano_complie_dir = pj(TMPDIR, '.theano-{cohort}-{scatter}')
    conda: CONDA_GATK_CNV
    resources:
            mem_mb = 25000,
            n = 4
    log: pj(LOG, '{cohort}.{scatter}.germlinecnvcalling.log')
    benchmark: pj(BENCH, '{cohort}.{scatter}.germlinecnvcalling.txt')
    shell:
        """
            mkdir -p {params.theano_complie_dir}
            export OMP_NUM_THREADS={resources.n} 
            THEANO_FLAGS="base_compiledir={params.theano_complie_dir}"  {params.java} -jar {params.gatk} GermlineCNVCaller {params.inputs} -L {input.scatters} --contig-ploidy-calls  {params.contig_ploydi_calls} --interval-merging-rule OVERLAPPING_ONLY --run-mode COHORT --output GATK_gCNV/{wildcards.cohort}_scatter_{wildcards.scatter} --output-prefix scatterd_{wildcards.cohort}_{wildcards.scatter} 2> {log}
        """

rule PostprocessGermlineCNVCalls:
    input:
        samples = input_func,
        models= expand(pj(GATK_gCNV,'{cohort}_scatter_{scatter}','scatterd_{cohort}_{scatter}-model','calling_config.json'),scatter = scatter_merged_cature_kit, allow_missing = True),
         # sample_index = expand(pj(GATK_gCNV, '{cohort}_scatter_{scatter}', 'scatterd_{cohort}_{scatter}-calls', 'SAMPLE_{index}', 'sample_name.txt'), scatter = scatter_merged_cature_kit, allow_missing = True),
    output:
        genotyped_intervals = pj(GATK_gCNV, 'GENOTYPED_CALLS_intervals_{cohort}', 'COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz'),
        genotyped_segments = pj(GATK_gCNV, 'GENOTYPED_CALLS_segments_{cohort}', 'COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz'),
        denoised_copy_ratio= pj(GATK_gCNV,'GENOTYPED_CALLS_denoised_copy_ratio_{cohort}','COHORT_{cohort}_SAMPLE_{sample}_{index}.tsv')
    conda: CONDA_GATK_CNV
    params:
        java= java_cnv,
        gatk= gatk_cnv,
        model_shrads = expand(' --model-shard-path GATK_gCNV/{cohort}_scatter_{scatter}/scatterd_{cohort}_{scatter}-model ',  scatter = scatter_merged_cature_kit, allow_missing = True),
        calls_shrads = expand(' --calls-shard-path GATK_gCNV/{cohort}_scatter_{scatter}/scatterd_{cohort}_{scatter}-calls ',  scatter = scatter_merged_cature_kit, allow_missing = True),
        CPC = (pj(GATK_gCNV,  '{cohort}-calls')),
        SD = REF_MALE_DICT,
    benchmark: pj(BENCH, '{cohort}.{sample}.{index}.PostprocessGermlineCNVcalls.txt')
    resources:
            mem_mb = 6000
    shell:
            """
            {params.java} -jar {params.gatk} PostprocessGermlineCNVCalls --sequence-dictionary {params.SD}  {params.model_shrads} {params.calls_shrads} --contig-ploidy-calls {params.CPC} --sample-index {wildcards.index} --allosomal-contig chrX --allosomal-contig chrY --output-genotyped-intervals {output.genotyped_intervals} --output-genotyped-segments {output.genotyped_segments} --output-denoised-copy-ratios {output.denoised_copy_ratio}
            """