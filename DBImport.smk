import os
from common import *
import yaml

wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    mode = r"WES|WGS"


module Tools:
    snakefile: 'Tools.smk'
    config: config
use rule * from Tools

gvcf_caller = config.get("caller", "HaplotypeCaller") #not implemented yet

genotype_mode = config.get("genotype_mode", "WES") #or WGS
genotype_level = int(config.get("genotype_level", 3))
DBIpath = config.get("DBIpath", "genomicsdb_")

if genotype_level == 2:
    parts = get_regions(level2_range_so)
elif genotype_level == 3:
    parts = get_regions(level3_range_so)
elif genotype_level == 4:
    parts = get_regions(level4_range_so)
else:
    raise RuntimeError(f'Unknown level {genotype_level}')



print(f"Sample type: {genotype_mode}")
print(f"Genotype level: {genotype_level}")

def region_to_IL_file(wildcards):#{{{
    """Converts a region to a interval_list file location (see common.py and Tools.smk)"""
    region = wildcards['region']
    # WGS files have fewer regions so DBI works faster and could use multiple cores
    return region_to_file(region, wgs=True, classic=True, extension='interval_list')#}}}


rule DBImport_all:
    input:
        expand(['labels/done_p{region}.txt'],region = parts),
        # expand("{chr}_gvcfs.list", chr = main_chrs)
    default_target: True

def generate_gvcf_input(wildcards):
    region = wildcards['region']
    res = []
    for samplefile in SAMPLE_FILES:
        sample_names = SAMPLEFILE_TO_SAMPLES[samplefile]
        samplefile_folder = get_samplefile_folder(samplefile)
        sample_sex = read_sexchrom(pj(samplefile_folder, samplefile + '.sex_chrom.tab'))
        gvcf_input = []
        for sample in sample_names:
            if sample_sex[sample] == 'F' and (region.startswith('Y') or region.endswith('H')):
                continue
            # Determine if it is WGS or WES
            if 'wgs' in SAMPLEINFO[sample]["sample_type"] or "WGS" in SAMPLEINFO[sample]["sample_type"] :
                chunk = convert_to_level1(region)
                if genotype_mode == 'WES':
                    filenames = expand("{cd}/{GVCF}/exome_extract/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=GVCF,region = chunk, sample=sample,allow_missing=True)
                else:
                    filenames = expand("{cd}/{GVCF}/reblock/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=GVCF,region = chunk, sample=sample,allow_missing=True)
            else:  # WES
                chunk = convert_to_level0(region)
                filenames = expand("{cd}/{GVCF}/reblock/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=GVCF,region = chunk, sample=sample,allow_missing=True)
            gvcf_input.extend(filenames)
        res.extend(gvcf_input)
    return res

rule GenomicDBImport:
    input:
        g=generate_gvcf_input,
        intervals = region_to_IL_file
    conda: CONDA_VCF
    output:
        ready=touch('labels/done_p{region}.txt'),
        dbi=directory(DBIpath + "p{region}")
    threads: 3
    params:
        inputs=lambda wildcards,input: ' '.join([f'-V {gvcf}' for gvcf in input.g]),
        batches='75',
        ref = REF_MALE,
        merge_contigs=lambda wildcards: ' --merge-contigs-into-num-partitions 1 '  if 'O' in wildcards['region'] else ''
    priority: 30
    resources: 
        n="10",
        mem_mb = lambda wildcards, attempt: attempt*20500,
        mem_mb_reduced = lambda wildcards, attempt: attempt * 16500, #tile db is not included in java memory
        tmpdir= TMPDIR
    shell:
        """
            {gatk} GenomicsDBImport --java-options "-Xmx{resources.mem_mb_reduced}M"  --reader-threads {threads} {params.inputs}  --consolidate True --max-num-intervals-to-import-in-parallel {threads} \
            --intervals {input.intervals} {params.merge_contigs} -R {params.ref} --genomicsdb-workspace-path {output.dbi}/ --batch-size {params.batches} --tmp-dir {resources.tmpdir} --merge-input-intervals \
         --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader"""



