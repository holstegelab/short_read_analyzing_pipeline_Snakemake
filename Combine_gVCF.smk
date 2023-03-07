configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")
gatk = config['gatk']
samtools = config['samtools']
bcftools = config['bcftools']
dragmap = config['dragmap']
verifybamid2 = config['verifybamid2']

ref = config['RES'] + config['ref']

wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+",

from read_samples import *
from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config

module gVCF:
    snakefile: 'gVCF.smk'
    config: config
use rule * from gVCF

# rule all:
#     expand('{samplefile}.done', samplefile=SAMPLE_FILES)
mode = config.get("computing_mode", "WES")

rule Combine_gVCF_all:
    input:
        expand(["{gvcf}/MERGED/bin_level/{chr}_{chr_p}.{mode}.g.vcf.gz"], zip, chr = main_chrs_db, chr_p = chr_p, gvcf = [config['gVCF']]*853, mode = mode*853),
    default_target: True


def get_mem_mb_combine_gvcf(wildcrads, attempt):
    return attempt*(config['combinegvcfs']['mem'])

rule combinegvcfs:
    input:
        gvcf = expand("{gvcf}/reblock/{chr}/{sample}.{chr}.{mode}.g.vcf.gz", gvcf = config['gVCF'], sample = sample_names, allow_missing=True),
        intervals = os.path.join(config['RES'], config['kit_folder'], 'BINS', 'interval_list', '{chr}_{chr_p}.interval_list')
    output: chr_gvcfs = config['gVCF'] + "/MERGED/bin_level/{chr}_{chr_p}.{mode}.g.vcf.gz"
    log: combine =config['LOG'] + "/{chr}_{chr_p}_{mode}.combine.log"
    benchmark:
        config['BENCH'] + "/{chr}_{chr_p}_{mode}.combinegvcf_NV.txt"
    conda: "envs/preprocess.yaml"
    priority: 30
    resources: mem_mb = get_mem_mb_combine_gvcf
    params:
        inputs = expand("--variant {gvcf}/reblock/{chr}/{sample}.{chr}.{mode}.g.vcf.gz", gvcf = config['gVCF'], sample = sample_names, allow_missing=True),
    shell:
        """{gatk} CombineGVCFs --java-options "-Xmx{resources.mem_mb}M"  -G StandardAnnotation -G AS_StandardAnnotation {params.inputs} -O {output} -R {ref} -L {input.intervals} 2> {log}"""
