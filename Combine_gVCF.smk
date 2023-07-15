configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")

wildcard_constraints:
    sample="[\w\d_\-@]+",
    readgroup="[\w\d_\-@]+",

from common import *


module Aligner:
    snakefile: 'Aligner.smk'
    config: config

module gVCF:
    snakefile: 'gVCF.smk'
    config: config
use rule * from gVCF


mode = config.get("computing_mode", "WES")

rule Combine_gVCF_all:
    input:
        expand(os.path.join(GVCF, "MERGED/{samplefile}.{region}.wg.vcf.gz"), samplefile=SAMPLE_FILES, region=level1_regions),
    default_target: True


def get_input_combine_files(wildcards):
    samples = SAMPLEFILE_TO_SAMPLES[wildcards.samplefile]
    res = []
    for sample in samples:
        
        if 'wgs' in SAMPLEINFO[sample]['sample_type']:
            res.append(pj(GVCF, wildcards.region,  sample + '.' + wildcards.region + '.wg.vcf.gz'))
        else:
            region = convert_to_level0(wildcards.region)
            res.append(pj(GVCF, region,  sample + '.' + region + '.wg.vcf.gz'))

    return res


rule combinegvcfs:
    """Create a combined gVCF file for each chromosome per batch of samples"""
    input:
        gvcf = get_input_combine_files,
        interval =lambda wildcards: region_to_file(wildcards.region, wgs=True, extension='interval_list')
    output: gvcfs = pj(GVCF, "MERGED/{samplefile}.{region}.wg.vcf.gz")
    log: combine = pj(LOG, "{samplefile}.{region}.combine.log")
    benchmark: pj(BENCH, "{samplefile}.{region}.combinegvcf_NV.txt")
    conda: CONDA_VCF
    priority: 30
    resources: 
        mem_mb =lambda wildcards, attempt: attempt*7000,
        n = 1
    params:
        inputs = lambda wildcards, input: [f"--variant {file}" for file in input.gvcf]
        
    shell:
        """{gatk} CombineGVCFs --java-options "-Xmx{resources.mem_mb}M"  -G StandardAnnotation -G AS_StandardAnnotation {params.inputs} -O {output} -R {REF_MALE} -L {input.interval} 2> {log}"""


