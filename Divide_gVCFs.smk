from common import *


def get_deepvariant_files(wildcards):#{{{
    sample = wildcards['sample']
    if 'wgs' in SAMPLEINFO[sample]['sample_type']:
        return [pj(DEEPVARIANT,  'gVCF', 'exome_extract', region, f'{sample}.{region}.wg.vcf.gz') for region in level1_regions]
    else:
        return [pj(DEEPVARIANT,  'gVCF', region, f'{sample}.{region}.wg.vcf.gz') for region in level0_regions]#}}}

rule Divide_gVCFs_all:
    input: expand(pj(DEEPVARIANT, 'gVCF', 'DIVIDED', '{region}', '{sample}.{region}.wg.vcf.gz'), region=level2_regions, sample=sample_names)
    default_target: True


rule divide_deepvariant:
    input: get_deepvariant_files
    output: expand(pj(DEEPVARIANT, 'gVCF', 'DIVIDED', '{region}', '{sample}.{region}.wg.vcf.gz'), region=level2_regions, allow_missing = True)
    conda: CONDA_MAIN
    run:
        for region in level2_regions:
            bed_file = region_to_file(region, wgs=True, extension='bed')
            shell(f'bcftools view -r {bed_file} {input} -O z -o deepvariant/gVCF/DIVIDED/{region}/{wildcards.sample}.{region}.wg.vcf.gz')