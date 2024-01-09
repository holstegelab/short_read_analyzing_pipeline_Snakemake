import os
from common import *

wildcard_constraints:
    sample="[\w\d_\-@]+",
    mode = "WES|WGS"

gvcf_caller = config.get("caller", "HaplotypeCaller")
sample_types = config.get("sample_types","WGS")
DBImethod = config.get("DBI_method", "new")
DBIpath = config.get("DBIpath", "genomicsdb_")

parts = level2_regions



print(f"Caller: {gvcf_caller}")
print(f"Sample type: {sample_types}")
print(f"DBImethod: {DBImethod}, DBIpath: {DBIpath}")

regions = []
for part in parts:
    regions.append(convert_to_level0(part))

print(parts, regions)

if DBImethod == "new":
    # if want to
    DBI_method_params = "--genomicsdb-workspace-path "
    path_to_dbi = DBIpath
    labels = []

elif DBImethod == "update" and len(DBIpath) != 1:
    
    DBI_method_params = "--genomicsdb-update-workspace-path "
    path_to_dbi = DBIpath
    number_of_splits = len(regions)
    labels = expand(["labels/done_backup_{samplefile}_{region}.p{part}"],zip, part = parts,samplefile=SAMPLE_FILES * number_of_splits)

elif DBImethod == "update" and len(DBIpath) == 1:
    raise ValueError(
        "If you want to update existing DB please provide path to this DB in format 'DBIpath=/path/to/directory_with_DB-s/genomicsdb_' \n"
        "Don't provide {chr}.p{chr_p} part of path!"
    )

else:
    raise ValueError(
        "invalid option provided to 'DBImethod'; please choose either 'new' or 'update'."
    )


def region_to_IL_file(wildcards):#{{{
    """Converts a region to a interval_list file location (see common.py and Tools.smk)"""
    region = wildcards['part']
    # WGS files have fewer regions so DBI works faster and could use multiple cores
    return region_to_file(region, wgs=sample_types=='WGS', extension='interval_list')#}}}


rule DBImport_all:
    input:
        expand(['labels/done_p{part}.txt'],part = parts),
        # expand("{chr}_gvcfs.list", chr = main_chrs)
    default_target: True

rule backup_gdbi:
    input: gdbi = path_to_dbi + 'p{part}'
    output: label = touch(temp('labels/done_backup_{samplefile}_p{part}'))
    params: tar = "{samplefile}_gdbi_p{part}.tar.gz"
    shell: """
            mkdir -p BACKUPS/previous &&
            find . -maxdepth 2 -name '*_gdbi_p{part}.tar.gz' -type f -print0 | xargs -0r mv -t BACKUPS/previous/ && 
            tar -czv -f BACKUPS/{params.tar} {input}
            """
def generate_gvcf_input(gvcf_folder, part):
    res = []
    for samplefile in SAMPLE_FILES:
        sample_names = SAMPLEFILE_TO_SAMPLES[samplefile]
        samplefile_folder = get_samplefile_folder(samplefile)
        gvcf_input = []
        for sample in sample_names:
            #determine if it is wgs or wes
            if SAMPLEINFO[sample]["sample_type"] == "WGS":
                sex = SAMPLEINFO[sample]["sex"]
                if sex == "M" or not part.startswith('Y'):
                    region = convert_to_level1(part)
                else:
                    region = []
            else: #wes
                sex = SAMPLEINFO[sample]["sex"]
                if sex == "M" or not part.startswith('Y'):
                    region = convert_to_level0(part)
                else:
                    region = []
            filename = expand("{cd}/{GVCF}/{region}/{sample}.{region}.wg.vcf.gz",cd=samplefile_folder,GVCF=gvcf_folder,region = region, sample=sample_names,allow_missing=True)
            gvcf_input.append(filename)
        res.extend(gvcf_input)
    return res

# gvcf_input = generate_gvcf_input(GVCF + "/exome_gatk", wildcards.part)

rule GenomicDBImport:
    input:
        g=lambda wildcards: generate_gvcf_input(GVCF + "/exome_gatk", wildcards.part),
        labels = labels
    conda: CONDA_VCF
    output:
        ready=touch(temp('labels/done_p{part}.txt'))
    threads: 3
    params:
        inputs=lambda wildcards,input: ' '.join([f'-V {gvcf}' for gvcf in input.g]),
        dbi=os.path.join(path_to_dbi + "p{part}"),
        method=DBI_method_params,
        batches='75',
        intervals = region_to_IL_file,
        ref = REF_MALE
    priority: 30
    resources: 
        n="3",
        mem_mb = lambda wildcards, attempt: attempt*8500,
        mem_mb_reduced = lambda wildcards, attempt: attempt * 6500, #tile db is not included in java memory
        tmpdir= TMPDIR
    shell:
        """{gatk} GenomicsDBImport --java-options "-Xmx{resources.mem_mb_reduced}M"  --reader-threads {threads} {params.inputs}  --consolidate True --max-num-intervals-to-import-in-parallel {threads} \
            --intervals {params.intervals}  -R {params.ref} {params.method} {params.dbi}/ --batch-size {params.batches} --tmp-dir {resources.tmpdir} --merge-input-intervals \
         --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader"""
