from common import *
import h5py


wildcard_constraints:
    sample="[\w\d_\-@]+",
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
use rule * from Tools

def get_sample_name(sample_index, cohort):
    # Get the sample name from the "sample_name.txt" file
    sample_name_file = os.path.join(GATK_gCNV, f'{cohort}_scatter_0001', f'scatterd_{cohort}_0001-calls', f'SAMPLE_{sample_index}', 'sample_name.txt')
    with open(sample_name_file, 'r') as f:
        sample_name = f.read().strip()
    return sample_name


def get_groups_from_hdf5(hdf5_file):
    with h5py.File(hdf5_file, 'r') as f:
        groups = list(f['coverage']['cluster_3_cor_cov'])
    return groups

hdf5_files = [f"{SF}.coverage.hdf5" for SF in SAMPLE_FILES]
groups = set()
for hdf5_file in hdf5_files:
    groups.update(get_groups_from_hdf5(hdf5_file))
# Filter out the '-1' group
filtered_groups = sorted([group for group in groups if group != '-1'])

def get_samples_in_group(hdf5_file, group):
    with h5py.File(hdf5_file, 'r') as f:
        samples = list(f['coverage']['samples'][f['coverage']['cluster_3_cor_cov'] == group])
    return samples

def get_samples_in_all_groups(hdf5_file):
    with h5py.File(hdf5_file, 'r') as f:
        samples = list(f['coverage']['samples'][f['coverage']['cluster_3_cor_cov'] != -1])
    return samples
grouped_samples = []
for hdf5_file in hdf5_files:
    grouped_samples.append(get_samples_in_all_groups(hdf5_file))


def input_func(wildcards):
    group = wildcards.group
    hdf5_file = wildcards.hdf5_file
    samples_in_group = get_samples_in_group(hdf5_file, group)
    return expand(
        '{gatk_gcnv}/Read_counts_hdf5/{cohort}/{sample}_readcounts.hdf5',
        gatk_gcnv=GATK_gCNV,
        cohort=hdf5_file,
        sample=samples_in_group,
        allow_missing=True
    )

def sample_list_per_cohort(wildcards):
    group = wildcards.group
    hdf5_file = wildcards.hdf5_file
    samples_in_group = get_samples_in_group(hdf5_file, group)
    return expand(
        ' -I {gatk_gcnv}/Read_counts_hdf5/{cohort}/{sample}_readcounts.hdf5 ',
        gatk_gcnv=GATK_gCNV,
        cohort=hdf5_file,
        sample=samples_in_group,
        allow_missing=True
    )

def generate_scatter_list(start, end):
    string_list = [f"{i:04d}_of_{end}" for i in range(int(start), int(end) + 1)]
    return string_list

scatter_merged_cature_kit = generate_scatter_list('0001', '0148')

rule gCNV_gatk_all:
    input:
        expand(pj(GATK_gCNV,'{cohort}_scatter_{scatter}','scatterd_{cohort}_{scatter}-model'),scatter=scatter_merged_cature_kit,cohort=filtered_groups),
        rules.Stat_all.input
    default_target: True


# expand("scatter_{part}", part = parts)

def get_capture_kit_path(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    if SAMPLEINFO[wildcards['sample']]['sample_type'].startswith('illumina_wgs'):
        capture_kit_path = None
    else:
        capture_kit_path = pj(INTERVALS_DIR, capture_kit + '.interval_list')
    return capture_kit_path

def get_capture_kit(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    return capture_kit

def get_preprocessed_capture_kit(wildcards):
    capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit']
    return pj(INTERVALS_DIR, 'preprocessed_intervals_for_GATK_CNV', capture_kit + '.interval_list')

rule collect_read_counts:
    input: bam = rules.markdup.output.mdbams,
            # Merged capture kit for test
            capture_kit = MERGED_CAPTURE_KIT_IVL
            # capture_kit= ancient(get_preprocessed_capture_kit),
    output: ReadCounts = pj(GATK_gCNV, 'Read_counts_hdf5', '{cohort}', '{sample}_readcounts.hdf5')
    params: java = java_cnv,
            gatk = gatk_cnv,
            ref = REF_MALE
    conda: CONDA_GATK_CNV
    log: pj(LOG, '{sample}.{cohort}.collectreadcounts_gatkcnv.log')
    benchmark: pj(BENCH, '{sample}.{cohort}.collectreadcounts_gatkcnv.txt')
    shell: """
            {params.java} -jar {params.gatk} CollectReadCounts -L {input.capture_kit} -R {params.ref} -imr OVERLAPPING_ONLY -I {input.bam} -O {output.ReadCounts}
            """

rule filterintervals:
    input: input_func,

    output: filtered_intervals = pj(GATK_gCNV, 'filtred_intervals', '{cohort}_filtred.interval_list')
    params: inputs = sample_list_per_cohort,
            java= java_cnv,
            gatk=gatk_cnv,
            # function to choose capture kit based on cohort
            # merged CDS for test only!
            capture_kit = pj(INTERVALS_DIR, 'preprocessed_intervals_for_GATK_CNV', 'merged_capture_kits_cds.interval_list' )
    conda: CONDA_GATK_CNV
    log: pj(LOG, '{cohort}.filterintervals_gatkcnv.log')
    benchmark: pj(BENCH, '{cohort}.filterintervals_gatkcnv.txt')
    shell: """
            {params.java} -jar {params.gatk} FilterIntervals -L {params.capture_kit} {params.inputs} -imr OVERLAPPING_ONLY -O {output.filtered_intervals}
            """

rule DetermineGCP:
    input: samples = input_func,
            intervals = rules.filterintervals.output.filtered_intervals
    output: CPC = dir(pj(GATK_gCNV, '/{cohort}-calls')),
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            capture_kit= pj(INTERVALS_DIR,'preprocessed_intervals_for_GATK_CNV','merged_capture_kits_cds.interval_list' ),
            contig_ploydi_priors = PL_PR_TABLE
    conda: CONDA_GATK_CNV
    log: pj(LOG, '{cohort}.determinecontigploydi.log')
    benchmark: pj(BENCH, '{cohort}.determinecontigploydi.txt')
    shell: """
            {params.java} -jar {params.gatk} DetermineGermlineContigPloidy --output gCNV/ --output-prefix {cohort} {params.inputs} -L {input.intervals} -imr OVERLAPPING_ONLY --contig-ploidy-priors {params.contig_ploydi_priors}
            """



rule GermlineCNVCaller:
    input: samples = input_func,
            scatters = pj(INTERVALS_DIR, 'scatter_merged_capture_kits_cds', 'temp_{scatter}'),
            contig_ploudi_calls = rules.DetermineGCP.output.CPC,
    output: OD = dir(pj(GATK_gCNV, '{cohort}_scatter_{scatter}')),
            calls = dir(pj(GATK_gCNV, '{cohort}_scatter_{scatter}', 'scatterd_{cohort}_{scatter}-calls')),
            models = dir(pj(GATK_gCNV, '{cohort}_scatter_{scatter}', 'scatterd_{cohort}_{scatter}-model'))
    params: inputs = sample_list_per_cohort,
            java = java_cnv,
            gatk = gatk_cnv,
            capture_kit= pj(INTERVALS_DIR,'preprocessed_intervals_for_GATK_CNV','merged_capture_kits_cds.interval_list' ),
    conda: CONDA_GATK_CNV
    log: pj(LOG, '{cohort}.{scatter}.germlinecnvcalling.log')
    benchmark: pj(BENCH, '{cohort}.{scatter}.germlinecnvcalling.txt')
    shell: """
           {params.java} -jar {params.gatk} GermlineCNVCaller {params.inputs} -L {input.scatters} --contig-ploidy-calls  {input.contig_ploudi_calls} --interval-merging-rule OVERLAPPING_ONLY --run-mode COHORT --output {output.OD} --output-prefix scatterd_{cohort}_{scatter}
    """

# rule PostprocessGermlineCNVCalls:
#     input:
#         models = expand(pj(GATK_gCNV, '{cohort}_scatter_{scatter}', 'scatterd_{cohort}_{scatter}-model'), scatter = scatter_merged_cature_kit, allow_missing = True),
#         calls = expand(pj(GATK_gCNV, '{cohort}_scatter_{scatter}', 'scatterd_{cohort}_{scatter}-calls'), scatter = scatter_merged_cature_kit, allow_missing = True),
#         sample_index = pj(GATK_gCNV, '{cohort}_scatter_0001', 'scatterd_{cohort}_0001-calls', 'SAMPLE_{index}'),
#     output:
#         genotyped_intervals = pj(GATK_gCNV, 'GENOTYPED_CALLS_{cohort}', 'COHORT_{cohort}_SAMPLE_{sample}')
#     params:
#         sample = get_sample_name(str(wildcards.index).zfill(4), wildcards.cohort)  # Add leading zeros to index
#     wildcard_constraints:
#         cohort = '|'.join(groups),
#         index = '|'.join([str(idx).zfill(4) for idx in get_samples_in_group(hdf5_files[0], groups)])
#     run:

#     how to come from INDEX to SAMPLE NAME?
#     function to read a tsv file and get a sample name from it
#     how to iterate over indexes?




# ~/data/hg38_res/software/java/jdk-17.0.7/bin/java -jar ~/data/hg38_res/software/gatk_4.4/build/bundle-files-collected/gatk-package-4.4.0.0-27-gabe8148-SNAPSHOT-local.jar GermlineCNVCaller
#
# java = '~/data/hg38_res/software/java/jdk-17.0.7/bin/java'
# gatk = '~/data/hg38_res/software/gatk_4.4/build/bundle-files-collected/gatk-package-4.4.0.0-27-gabe8148-SNAPSHOT-local.jar'
#
# # java = os.path.join(config['RES'], config['SOFTWARE'], config['JAVA_for_GATK'])
# # gatk = os.path.join(config['RES'], config['SOFTWARE'], config['GATK'])
#
# rule gCNV_scatter:
#     input: scatter = "scatter/temp_{part}_of_90/scattered.interval_list"
#     output: touch("scatter_{part}")
#     shell: "{java} -jar {gatk} GermlineCNVCaller --input tsvs/NL_VUMC_KG-001369.markdup.tsv --input tsvs/NL_VUMC_KG-001374.markdup.tsv --input tsvs/NL_VUMC_KG-001436.markdup.tsv --input tsvs/NL_VUMC_KG-001474.markdup.tsv --input tsvs/NL_VUMC_KG-001492.markdup.tsv --input tsvs/NL_VUMC_KG-001515.markdup.tsv --input tsvs/NL_VUMC_KG-001532.markdup.tsv --input tsvs/NL_VUMC_KG-001584.markdup.tsv --input tsvs/NL_VUMC_KG-001661.markdup.tsv --input tsvs/NL_VUMC_KG-001769.markdup.tsv --input tsvs/NL_VUMC_KG-001781.markdup.tsv --input tsvs/NL_VUMC_KG-001883.markdup.tsv --input tsvs/NL_VUMC_KG-001884.markdup.tsv --input tsvs/NL_VUMC_KG-002070.markdup.tsv --input tsvs/NL_VUMC_KG-002165.markdup.tsv --input tsvs/NL_VUMC_KG-002169.markdup.tsv --input tsvs/NL_VUMC_KG-002246.markdup.tsv --input tsvs/NL_VUMC_KG-002271.markdup.tsv --input tsvs/NL_VUMC_KG-002273.markdup.tsv --input tsvs/NL_VUMC_KG-002320.markdup.tsv --input tsvs/NL_VUMC_KG-002335.markdup.tsv --input tsvs/NL_VUMC_KG-002347.markdup.tsv --input tsvs/NL_VUMC_KG-002348.markdup.tsv --input tsvs/NL_VUMC_KG-002482.markdup.tsv --input tsvs/NL_VUMC_KG-002510.markdup.tsv --input tsvs/NL_VUMC_KG-002511.markdup.tsv --input tsvs/NL_VUMC_KG-002522.markdup.tsv --input tsvs/NL_VUMC_KG-002564.markdup.tsv --input tsvs/NL_VUMC_KG-002592.markdup.tsv --input tsvs/NL_VUMC_KG-002682.markdup.tsv --input tsvs/NL_VUMC_KG-002694.markdup.tsv --input tsvs/NL_VUMC_KG-002757.markdup.tsv --input tsvs/NL_VUMC_KG-002768.markdup.tsv --input tsvs/NL_VUMC_KG-002851.markdup.tsv --input tsvs/NL_VUMC_KG-002852.markdup.tsv --input tsvs/NL_VUMC_KG-002863.markdup.tsv --input tsvs/NL_VUMC_KG-002868.markdup.tsv --input tsvs/NL_VUMC_KG-002880.markdup.tsv --input tsvs/NL_VUMC_KG-002900.markdup.tsv --input tsvs/NL_VUMC_KG-002917.markdup.tsv --input tsvs/NL_VUMC_KG-002957.markdup.tsv --input tsvs/NL_VUMC_KG-002958.markdup.tsv --input tsvs/NL_VUMC_KG-002962.markdup.tsv --input tsvs/NL_VUMC_KG-003940.markdup.tsv --input tsvs/NL_VUMC_KG-003947.markdup.tsv --input tsvs/NL_VUMC_KG-003950.markdup.tsv --input tsvs/NL_VUMC_KG-003951.markdup.tsv --input tsvs/NL_VUMC_KG-003974.markdup.tsv --input tsvs/NL_VUMC_KG-004019.markdup.tsv --input tsvs/NL_VUMC_KG-004025.markdup.tsv --input tsvs/NL_VUMC_KG-004047.markdup.tsv --input tsvs/NL_VUMC_KG-004050.markdup.tsv --input tsvs/NL_VUMC_KG-004060.markdup.tsv --input tsvs/NL_VUMC_KG-004086.markdup.tsv --input tsvs/NL_VUMC_KG-004140.markdup.tsv --input tsvs/NL_VUMC_KG-004172.markdup.tsv --input tsvs/NL_VUMC_KG-004214.markdup.tsv --input tsvs/NL_VUMC_KG-004244.markdup.tsv --input tsvs/NL_VUMC_KG-004251.markdup.tsv --input tsvs/NL_VUMC_KG-004257.markdup.tsv --input tsvs/NL_VUMC_KG-004615.markdup.tsv --input tsvs/NL_VUMC_KG-004626.markdup.tsv --input tsvs/NL_VUMC_KG-004785.markdup.tsv --input tsvs/NL_VUMC_KG-004787.markdup.tsv --input tsvs/NL_VUMC_KG-004792.markdup.tsv --input tsvs/NL_VUMC_KG-004810.markdup.tsv --input tsvs/NL_VUMC_KG-004900.markdup.tsv --input tsvs/NL_VUMC_KG-004965.markdup.tsv --input tsvs/NL_VUMC_KG-004973.markdup.tsv --input tsvs/NL_VUMC_KG-004984.markdup.tsv --input tsvs/NL_VUMC_KG-004986.markdup.tsv --input tsvs/NL_VUMC_KG-004988.markdup.tsv --input tsvs/NL_VUMC_KG-005037.markdup.tsv --input tsvs/NL_VUMC_KG-005130.markdup.tsv --input tsvs/NL_VUMC_KG-005181.markdup.tsv --input tsvs/NL_VUMC_KG-005240.markdup.tsv --input tsvs/NL_VUMC_KG-005277.markdup.tsv --input tsvs/NL_VUMC_KG-005289.markdup.tsv --input tsvs/NL_VUMC_KG-005299.markdup.tsv --input tsvs/NL_VUMC_KG-005316.markdup.tsv --input tsvs/NL_VUMC_KG-005420.markdup.tsv --input tsvs/NL_VUMC_KG-005454.markdup.tsv --input tsvs/NL_VUMC_KG-005469.markdup.tsv --input tsvs/NL_VUMC_KG-005571.markdup.tsv --input tsvs/NL_VUMC_KG-005594.markdup.tsv --input tsvs/NL_VUMC_KG-005626.markdup.tsv --input tsvs/NL_VUMC_KG-005713.markdup.tsv --input tsvs/NL_VUMC_KG-005731.markdup.tsv --input tsvs/NL_VUMC_KG-005741.markdup.tsv --input tsvs/NL_VUMC_KG-005773.markdup.tsv --input tsvs/NL_VUMC_KG-005803.markdup.tsv --input tsvs/NL_VUMC_KG-005830.markdup.tsv --input tsvs/NL_VUMC_KG-005831.markdup.tsv --input tsvs/NL_VUMC_KG-005890.markdup.tsv --input tsvs/NL_VUMC_KG-005974.markdup.tsv --input tsvs/NL_VUMC_KG-006031.markdup.tsv --input tsvs/NL_VUMC_KG-006060.markdup.tsv --input tsvs/NL_VUMC_KG-006086.markdup.tsv --input tsvs/NL_VUMC_KG-006087.markdup.tsv --input tsvs/NL_VUMC_KG-006117.markdup.tsv --input tsvs/NL_VUMC_KG-006129.markdup.tsv --input tsvs/NL_VUMC_KG-006136.markdup.tsv --input tsvs/NL_VUMC_KG-006157.markdup.tsv --input tsvs/NL_VUMC_KG-006302.markdup.tsv --input tsvs/NL_VUMC_KG-006340.markdup.tsv --input tsvs/NL_VUMC_KG-006378.markdup.tsv --input tsvs/NL_VUMC_KG-006383.markdup.tsv --input tsvs/NL_VUMC_KG-006551.markdup.tsv --input tsvs/NL_VUMC_KG-006670.markdup.tsv --input tsvs/NL_VUMC_KG-006684.markdup.tsv --input tsvs/NL_VUMC_KG-006692.markdup.tsv --input tsvs/NL_VUMC_KG-006736.markdup.tsv --input tsvs/NL_VUMC_KG-006752.markdup.tsv --input tsvs/NL_VUMC_KG-006778.markdup.tsv --input tsvs/NL_VUMC_KG-006794.markdup.tsv --input tsvs/NL_VUMC_KG-006835.markdup.tsv --input tsvs/NL_VUMC_KG-006970.markdup.tsv --input tsvs/NL_VUMC_KG-007036.markdup.tsv --input tsvs/NL_VUMC_KG-007082.markdup.tsv --input tsvs/NL_VUMC_KG-007167.markdup.tsv --input tsvs/NL_VUMC_KG-007210.markdup.tsv --input tsvs/NL_VUMC_KG-007226.markdup.tsv --input tsvs/NL_VUMC_KG-007241.markdup.tsv --input tsvs/NL_VUMC_KG-007245.markdup.tsv --input tsvs/NL_VUMC_KG-007255.markdup.tsv --input tsvs/NL_VUMC_KG-007265.markdup.tsv --input tsvs/NL_VUMC_KG-007327.markdup.tsv --input tsvs/NL_VUMC_KG-007338.markdup.tsv --input tsvs/NL_VUMC_KG-007419.markdup.tsv --input tsvs/NL_VUMC_KG-007483.markdup.tsv --input tsvs/NL_VUMC_KG-007510.markdup.tsv --input tsvs/NL_VUMC_KG-007542.markdup.tsv --input tsvs/NL_VUMC_KG-007570.markdup.tsv --input tsvs/NL_VUMC_KG-007573.markdup.tsv --input tsvs/NL_VUMC_KG-007620.markdup.tsv --input tsvs/NL_VUMC_KG-007621.markdup.tsv --input tsvs/NL_VUMC_KG-007685.markdup.tsv --input tsvs/NL_VUMC_KG-007692.markdup.tsv --input tsvs/NL_VUMC_KG-007696.markdup.tsv --input tsvs/NL_VUMC_KG-007725.markdup.tsv --input tsvs/NL_VUMC_KG-007745.markdup.tsv --input tsvs/NL_VUMC_KG-007779.markdup.tsv --input tsvs/NL_VUMC_KG-007782.markdup.tsv --input tsvs/NL_VUMC_KG-007815.markdup.tsv --input tsvs/NL_VUMC_KG-007873.markdup.tsv --input tsvs/NL_VUMC_KG-007991.markdup.tsv --input tsvs/NL_VUMC_KG-008119.markdup.tsv --input tsvs/NL_VUMC_KG-008131.markdup.tsv --input tsvs/NL_VUMC_KG-008141.markdup.tsv --input tsvs/NL_VUMC_KG-008266.markdup.tsv --input tsvs/NL_VUMC_KG-008285.markdup.tsv --input tsvs/NL_VUMC_KG-008294.markdup.tsv --input tsvs/NL_VUMC_KG-008527.markdup.tsv --input tsvs/NL_VUMC_KG-008872.markdup.tsv --input tsvs/NL_VUMC_KG-008888.markdup.tsv --input tsvs/NL_VUMC_KG-008970.markdup.tsv --input tsvs/NL_VUMC_KG-008983.markdup.tsv --input tsvs/NL_VUMC_KG-009026.markdup.tsv --input tsvs/NL_VUMC_KG-009112.markdup.tsv --input tsvs/NL_VUMC_KG-009114.markdup.tsv --input tsvs/NL_VUMC_KG-009154.markdup.tsv --input tsvs/NL_VUMC_KG-009171.markdup.tsv --input tsvs/NL_VUMC_KG-009199.markdup.tsv --input tsvs/NL_VUMC_KG-009299.markdup.tsv --input tsvs/NL_VUMC_KG-009330.markdup.tsv --input tsvs/NL_VUMC_KG-009331.markdup.tsv --input tsvs/NL_VUMC_KG-009335.markdup.tsv --input tsvs/NL_VUMC_KG-009384.markdup.tsv --input tsvs/NL_VUMC_KG-009470.markdup.tsv --input tsvs/NL_VUMC_KG-009497.markdup.tsv --input tsvs/NL_VUMC_KG-009733.markdup.tsv --input tsvs/NL_VUMC_KG-009746.markdup.tsv --input tsvs/NL_VUMC_KG-009748.markdup.tsv --input tsvs/NL_VUMC_KG-009769.markdup.tsv --input tsvs/NL_VUMC_KG-009814.markdup.tsv --input tsvs/NL_VUMC_KG-009884.markdup.tsv --input tsvs/NL_VUMC_KG-010015.markdup.tsv --input tsvs/NL_VUMC_KG-010024.markdup.tsv --input tsvs/NL_VUMC_KG-010048.markdup.tsv --input tsvs/NL_VUMC_KG-010063.markdup.tsv --input tsvs/NL_VUMC_KG-010072.markdup.tsv --input tsvs/NL_VUMC_KG-010151.markdup.tsv --input tsvs/NL_VUMC_KG-010161.markdup.tsv --input tsvs/NL_VUMC_KG-010175.markdup.tsv --input tsvs/NL_VUMC_KG-010227.markdup.tsv --input tsvs/NL_VUMC_KG-010254.markdup.tsv --input tsvs/NL_VUMC_KG-010302.markdup.tsv --input tsvs/NL_VUMC_KG-010532.markdup.tsv --input tsvs/NL_VUMC_KG-010580.markdup.tsv --input tsvs/NL_VUMC_KG-010581.markdup.tsv --input tsvs/NL_VUMC_KG-010583.markdup.tsv --input tsvs/NL_VUMC_KG-010604.markdup.tsv --input tsvs/NL_VUMC_KG-010619.markdup.tsv --input tsvs/NL_VUMC_KG-010645.markdup.tsv --input tsvs/NL_VUMC_KG-010648.markdup.tsv --input tsvs/NL_VUMC_KG-010664.markdup.tsv --input tsvs/NL_VUMC_KG-010669.markdup.tsv --input tsvs/NL_VUMC_KG-010706.markdup.tsv --input tsvs/NL_VUMC_KG-010827.markdup.tsv -L {input} --contig-ploidy-calls test/test_ploidy-calls/ --interval-merging-rule OVERLAPPING_ONLY --verbosity DEBUG --output scatter_{wildcards.part} --output-prefix scatter_{wildcards.part}_90_cohort --run-mode COHORT"