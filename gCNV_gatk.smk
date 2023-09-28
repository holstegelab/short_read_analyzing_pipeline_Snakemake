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

scatter_merged_cature_kit = generate_scatter_list('0001', '148')

rule gCNV_gatk_all:
    input:
        expand('GATK_gCNV/{cohort}_scatter_{scatter}/scatterd_{cohort}_{scatter}-model/calling_config.json',scatter=scatter_merged_cature_kit,cohort=groups),
        expand('GATK_gCNV/filtred_intervals/{cohort}_filtred.interval_list', cohort=groups),
        expand('GATK_gCNV/{cohort}-calls/contig_ploidy_prior.tsv', cohort = groups),
        paths_output
        # expand('GATK_gCNV/GENOTYPED_CALLS_intervals_{cohort}/COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz')
        # rules.Stat_all.input
    default_target: True


rule extract_bams:
    input: bam=rules.markdup.output.mdbams,
    output: ex_bam = temp(pj(GATK_gCNV, 'extracted_bams','{sample}_extract.bam')),
            ex_bai = temp(ensure(pj(GATK_gCNV, 'extracted_bams','{sample}_extract.bai'), non_empty=True))
    conda: CONDA_MAIN
    params: reg = MAIN_CHRS_BED
    resources:
        n = 4
    shell:
        """
        samtools view -@ 3 -L {params.reg} -o {output.ex_bam} {input} && samtools index -@ 3 {output.ex_bam}
        """



rule collect_read_counts:
    input: bam = rules.extract_bams.output.ex_bam,
            # Merged capture kit for test
            capture_kit = MERGED_CAPTURE_KIT_IVL
    output: ReadCounts = ensure(pj(GATK_gCNV, 'Read_counts_hdf5', '{cohort}', '{sample}_readcounts.hdf5'), non_empty=True)
    params: java = java_cnv,
            gatk = gatk_cnv,
            ref = REF_MALE
    # wildcard_constraints:
    #     sample = "|".join(SAMPLE_TO_INDEX.keys()),
    #     index = "|".join(str(idx) for idx in SAMPLE_TO_INDEX.values())
    conda: CONDA_GATK_CNV
    resources:
            n = 2,
            mem_mb = 5500
    log: pj(LOG, '{sample}.{cohort}.collectreadcounts_gatkcnv.log')
    benchmark: pj(BENCH, '{sample}.{cohort}.collectreadcounts_gatkcnv.txt')
    shell: """
            {params.java} -jar {params.gatk} CollectReadCounts -L {input.capture_kit} -R {params.ref} -imr OVERLAPPING_ONLY -I {input.bam} -O {output.ReadCounts} 2> {log}
            """

rule filterintervals:
    input: input_func,
    output: filtered_intervals = pj(GATK_gCNV, 'filtred_intervals', '{cohort}_filtred.interval_list')
    params: inputs = sample_list_per_cohort,
            java= java_cnv,
            gatk=gatk_cnv,
            capture_kit = pj(INTERVALS_DIR, 'preprocessed_intervals_for_GATK_CNV', 'merged_capture_kits_cds.interval_list' )
    conda: CONDA_GATK_CNV
    log: pj(LOG, '{cohort}.filterintervals_gatkcnv.log')
    benchmark: pj(BENCH, '{cohort}.filterintervals_gatkcnv.txt')
    shell: """
            {params.java} -jar {params.gatk} FilterIntervals -L {params.capture_kit} {params.inputs} -imr OVERLAPPING_ONLY -O {output.filtered_intervals} 2> {log}
            """

rule DetermineGCP:
    input: samples = input_func,
            intervals = pj(GATK_gCNV, 'filtred_intervals', '{cohort}_filtred.interval_list')
    output: CPC = pj(GATK_gCNV,  '{cohort}-calls', 'contig_ploidy_prior.tsv'),
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            capture_kit= pj(INTERVALS_DIR,'preprocessed_intervals_for_GATK_CNV','merged_capture_kits_cds.interval_list'),
            contig_ploydi_priors = PL_PR_TABLE
    conda: CONDA_GATK_CNV
    # log: pj(LOG, '{cohort}.determinecontigploydi.log')
    # benchmark: pj(BENCH, '{cohort}.determinecontigploydi.txt')
    shell: """
            {params.java} -jar {params.gatk} DetermineGermlineContigPloidy  --output-prefix {wildcards.cohort} --output GATK_gCNV/ {params.inputs} -L {input.intervals} -imr OVERLAPPING_ONLY --contig-ploidy-priors {params.contig_ploydi_priors}
            """


rule GermlineCNVCaller:
    input:  samples = input_func,
            scatters = pj(INTERVALS_DIR, 'scatter_merged_capture_kits_cds', 'temp_{scatter}', 'scattered.interval_list'),
            contig_ploudi_calls = (pj(GATK_gCNV,  '{cohort}-calls', 'contig_ploidy_prior.tsv')),
    output: models = pj(GATK_gCNV,'{cohort}_scatter_{scatter}','scatterd_{cohort}_{scatter}-model', 'calling_config.json'),
            # sample_files = pj(GATK_gCNV, '{cohort}_scatter_{scatter}', 'scatterd_{cohort}_{scatter}-calls', 'SAMPLE_{index}', 'sample_name.txt')
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            capture_kit= pj(INTERVALS_DIR,'preprocessed_intervals_for_GATK_CNV','merged_capture_kits_cds.interval_list'),
    conda: CONDA_GATK_CNV
    # log: pj(LOG, '{cohort}.{scatter}.germlinecnvcalling.log')
    # benchmark: pj(BENCH, '{cohort}.{scatter}.germlinecnvcalling.txt')
    shell: """
           {params.java} -jar {params.gatk} GermlineCNVCaller {params.inputs} -L {input.scatters} --contig-ploidy-calls  {input.contig_ploudi_calls} --interval-merging-rule OVERLAPPING_ONLY --run-mode COHORT --output GATK_gCNV/{wildcards.cohort}_scatter_{wildcards.scatter} --output-prefix scatterd_{wildcards.cohort}_{wildcards.scatter}
    """

rule PostprocessGermlineCNVCalls:
    input:
        samples = input_func,
        models= expand(pj(GATK_gCNV,'{cohort}_scatter_{scatter}','scatterd_{cohort}_{scatter}-model','calling_config.json'),scatter = scatter_merged_cature_kit, allow_missing = True),
         # sample_index = expand(pj(GATK_gCNV, '{cohort}_scatter_{scatter}', 'scatterd_{cohort}_{scatter}-calls', 'SAMPLE_{index}', 'sample_name.txt'), scatter = scatter_merged_cature_kit, allow_missing = True),
    output:
        genotyped_intervals = pj(GATK_gCNV, 'GENOTYPED_CALLS_intervals_{cohort}', 'COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz'),
        genotyped_segments = pj(GATK_gCNV, 'GENOTYPED_CALLS_segments_{cohort}', 'COHORT_{cohort}_SAMPLE_{sample}_{index}.vcf.gz'),
    params:
        java= java_cnv,
        gatk= gatk_cnv,
        model_shrads = expand(' --model-shard-path GATK_gCNV/{cohort}_scatter_{scatter}/scatterd_{cohort}_{scatter}-model ',  scatter = scatter_merged_cature_kit, allow_missing = True),
        calls_shrads = expand(' --calls-shard-path GATK_gCNV/{cohort}_scatter_{scatter}/scatterd_{cohort}_{scatter}-calls ',  scatter = scatter_merged_cature_kit, allow_missing = True),
        CPC = (pj(GATK_gCNV,  '{cohort}-calls')),
        SD = REF_MALE_DICT
    shell:
            """
            {params.java} -jar {params.gatk} PostprocessGermlineCNVCalls --sequence-dictionary {params.SD}  {params.model_shrads} {params.calls_shrads} --contig-ploidy-calls {params.CPC} --sample-index {wildcards.index} --allosomal-contig chrX --allosomal-contig chrY --output-genotyped-intervals {output.genotyped_intervals} --output-genotyped-segments {output.genotyped_segments} 
            """