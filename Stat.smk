from common import *
import read_stats
import utils
import os
import tarfile
import csv
import glob
from shlex import quote
import math
import read_samples
import datetime

onsuccess: shell("rm -fr logs/Stats/*")

wildcard_constraints:
    sample=r"[\w\d_\-@]+"

module Aligner:
    snakefile: 'Aligner.smk'
    config: config

use rule * from Aligner

module Tools:
    snakefile: 'Tools.smk'
    config: config

use rule BedToIntervalList from Tools

##THIS FUNCION IS COPIED ALSO in Aligner.smk and Kraken.smk
def sampleinfo(SAMPLEINFO, sample, checkpoint=False):  #{{{
    """If samples are on tape, we do not have sample readgroup info.
    That is, the 'readgroups' field is empty.

    This function first checks if the readgroup info is available on disk,
    in the file SAMPLEINFODIR/<sample>.dat. 

    Alternatively, the function injects a checkpoint rule to load this readgroup info.
    """

    sinfo = SAMPLEINFO[sample]
    if not 'readgroups' in sinfo:
        rgpath = pj(SAMPLEINFODIR,sample + ".dat")
        if os.path.exists(rgpath):
            xsample = utils.load(rgpath)
        elif checkpoint:
            #no readgroup info yet
            filename = Aligner.checkpoints.get_readgroups.get(sample=sample).output[0]
            xsample = utils.load(filename)
        sinfo = sinfo.copy()
        sinfo['readgroups'] = xsample['readgroups']
        sinfo['alternative_names'] = sinfo.get('alternative_names',set()).union(xsample['alternative_names'])
        SAMPLEINFO[sample] = sinfo
    return sinfo  #}}}

rule Stat_all:
    input:
        expand("{samplefile}.oxo_quality.tab",samplefile=SAMPLE_FILES),
        expand("{samplefile}.bam_quality.tab",samplefile=SAMPLE_FILES),
        expand("{samplefile}.bam_rg_quality.tab",samplefile=SAMPLE_FILES),
        expand("{samplefile}.sex_chrom.tab",samplefile=SAMPLE_FILES),
        expand("{samplefile}.coverage.hdf5",samplefile=SAMPLE_FILES),
        expand("{samplefile}.kraken.tab",samplefile=SAMPLE_FILES),
        expand("{samplefile}.deepvariant_bcftools.tab",samplefile=SAMPLE_FILES),
        expand("{samplefile}.phase_quality.tab",samplefile=SAMPLE_FILES),
        # expand(pj(STAT,"{sample}.capture_kit_stats.tsv"), sample=sample_names),

def get_rg_files(wildcards):
    """Get sorted bam index files for all readgroups for a given sample."""
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []
    for readgroup in readgroups_b:
        files.append(pj(STAT, f"{wildcards['sample']}.{readgroup['info']['ID']}.dragmap.log"))
        files.append(pj(STAT, f"{wildcards['sample']}.{readgroup['info']['ID']}.adapter_removal.log"))
        files.append(pj(STAT, f"{wildcards['sample']}.{readgroup['info']['ID']}.fastq.adapters"))
        files.append(pj(STAT, f"{wildcards['sample']}.{readgroup['info']['ID']}.merge_stats.tsv"))
        files.append(pj(STAT, f"{wildcards['sample']}.{readgroup['info']['ID']}.dechimer_stats.tsv"))
       
    return files

rule tar_stats_per_sample:
    input:
        hs=pj(STAT,"{sample}.hs_metrics"),
        samtools=pj(STAT,"{sample}.samtools.stat"),
        exome=pj(STAT,'{sample}.samtools.exome.stat'),
        contam=pj(STAT,'contam','{sample}.verifybamid.pca2.selfSM'),
        bam_all=pj(STAT,'{sample}.bam_all.tsv'),
        bam_exome=pj(STAT,'{sample}.bam_exome.tsv'),
        pread=pj(STAT,'{sample}.pre_adapter_summary_metrics'),
        biat_bias=pj(STAT,'{sample}.bait_bias_summary_metrics'),
        pread_det=pj(STAT,'{sample}.pre_adapter_detail_metrics'),
        biat_bias_det=pj(STAT,'{sample}.bait_bias_detail_metrics'),
        cov=pj(STAT,'cov','{sample}.regions.bed.gz'),
        markdup=pj(STAT,"{sample}.markdup.stat"),
        sex_y=pj(KMER,"{sample}.result.yaml"),
        rg_logs=get_rg_files,
        error_summary=pj(STAT,"{sample}.error_summary_metrics"),
        mosdepth_dist=pj(STAT,'cov','{sample}.mosdepth.region.dist.txt'),
        ancestry=pj(STAT,"contam","{sample}.verifybamid.pca2.Ancestry"),
        chrM=pj(STAT,'{sample}.chrM_read_stats.tsv'),
        numt=pj(STAT,'{sample}.numt_read_stats.tsv'),
        phase=pj(STAT,'{sample}.phase_stats.tsv'),
    output:
        tar=temp(pj(STAT,"{sample}.stats.tar.gz"))
    resources:
        n="1",
        mem_mb=100
    shell: """
            tar -czvf {output.tar} {input.error_summary} {input.mosdepth_dist} {input.ancestry} {input.markdup} {input.hs} {input.samtools} {input.exome} {input.contam} {input.bam_all} {input.bam_exome} {input.pread} {input.biat_bias} {input.pread_det} {input.biat_bias_det} {input.cov} {input.sex_y} {input.rg_logs} {input.chrM} {input.numt} {input.phase}
            """

rule coverage:
    """Estimates coverage using the mosdepth tool"""
    input:
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
    output:
        pj(STAT,'cov','{sample}.regions.bed.gz'),
        pj(STAT,'cov','{sample}.regions.bed.gz.csi'),
        pj(STAT,'cov','{sample}.mosdepth.global.dist.txt'),
        pj(STAT,'cov','{sample}.mosdepth.summary.txt'),
        pj(STAT,'cov','{sample}.mosdepth.region.dist.txt')
    priority: 27
    params:
        bed=WINDOWS,
        prefix=pj(STAT,'cov','{sample}')
    resources:
        mem_mb=2200,
        n="1.5"
    conda: CONDA_MOSDEPTH
    shell:
        """
            mkdir -p `dirname {output[0]}`
            mosdepth  --threads 2 -b {params.bed} --no-per-base {params.prefix} {input.bam}
        """

rule chrM_and_numt_read_stats:
    input:
        chrm_bam=pj(chrM, '{sample}_chrM_orig.reads.bam'),
        chrm_bai=pj(chrM,'{sample}_chrM_orig.reads.bai'),
        numt_bam=pj(chrM, 'NUMTs', '{sample}_NUMTs.realign.bam'),
        numt_bai=pj(chrM, 'NUMTs', '{sample}_NUMTs.realign.bai')
    output:
        chrM=ensure(temp(pj(STAT, '{sample}.chrM_read_stats.tsv')), non_empty=True),
        numt=ensure(temp(pj(STAT, '{sample}.numt_read_stats.tsv')), non_empty=True)
    params:
        script=srcdir('scripts/region_read_stats.py')
    resources:
        n=1,
        mem_mb=200
    conda: CONDA_MAIN
    run:
        script = str(params.script)
        n = resources.n
        shell(f"python {quote(script)} --bam {quote(str(input.chrm_bam))} --threads {n} --exclude-flags 1024 --output {quote(str(output.chrM))}")
        shell(f"python {quote(script)} --bam {quote(str(input.numt_bam))} --threads {n} --exclude-flags 1024 --output {quote(str(output.numt))}")

def get_whatsHap_stats_inputs(wildcards):  #{{{
    sample = wildcards['sample']
    # For WGS, stats are produced for level0 regions (F, X, Y)
    # For WES, stats are produced for level1 regions (autosplit1 + X/Y)
    if 'wgs' in SAMPLEINFO[sample]['sample_type'] or 'WGS' in SAMPLEINFO[sample]['sample_type']:
        regions = level1_regions
    else:
        regions = level0_regions
    return [pj(STAT, 'whatshap_dvphasing', f"{sample}.{region}.stats") for region in regions]

#}}}

rule whatsHap_phase_stats:
    input:
        validated_sex=pj(KMER,"{sample}.result.yaml"),
        wstats=get_whatsHap_stats_inputs
    output:
        ensure(temp(pj(STAT,'{sample}.phase_stats.tsv')), non_empty=True)
    params:
        script=srcdir('scripts/whatsHap_phase_stats.py')
    resources:
        n=1,
        mem_mb=200
    conda: CONDA_MAIN
    run:
        files = [str(p) for p in input.wstats]
        inputs_q = ' '.join(quote(p) for p in files)
        sex_yaml = str(input.validated_sex)
        cmd = f"python {params.script} --sample {quote(wildcards.sample)} --sex-yaml {quote(sex_yaml)} --inputs {inputs_q} --output {quote(str(output[0]))}"
        shell(cmd)

rule copy_stats_tar_to_dcache:
    input:
        tar=pj(STAT, "{samplefile}.stats_bundle.tar.gz")
    output:
        copied=temp(pj(STAT, "{samplefile}.stats_bundle.tar.copied")),
        checksum=temp(pj(STAT, "{samplefile}.stats_bundle.tar.ADLER32"))
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1"
        ,dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        remote_dir = os.path.join(remote_base_for_samplefile(wildcards.samplefile), "stat")
        remote_name = os.path.basename(input.tar)

        copy_with_checksum(str(input.tar), remote_dir, remote_name, str(output.checksum), AGH_DCACHE_CONFIG, params.ada_script)
        shell(f"touch {quote(str(output.copied))}")

rule copy_excluded_to_dcache:
    input:
        samplefile=lambda wildcards: _samplefile_local_path(wildcards.samplefile),
        stats_tabs_copied=pj(STAT, "{samplefile}.stats_tabs.copied") #not used, but needed to wait for copy_samplefile_stats_to_dcache
    output:
        copied=temp(pj(STAT, "{samplefile}.excluded.copied")),
        checksums=temp(pj(STAT, "{samplefile}.excluded.ADLER32"))
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1"
        ,dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        basename = wildcards.samplefile
        excluded_map = SAMPLEFILE_TO_EXCLUDED_SAMPLES.get(basename, {})
        excluded_samples = [s for s, d in excluded_map.items() if d.get('info') is not None]

        if not excluded_samples:
            print(f"No excluded samples with info for {basename}, skipping uploads")
            # still write empty checksum file and touch copied marker
            with open(output.checksums, 'w') as sum_out:
                sum_out.write("")
            shell(f"touch {quote(str(output.copied))}")
            return

        try:
            base_target = remote_base_for_samplefile(basename)
        except Exception:
            all_map = SAMPLEFILE_TO_ALL_SAMPLES.get(basename, {})
            if not all_map:
                raise
            first = next(iter(all_map.keys()))
            sinfo_fb = all_map[first]
            target = sinfo_fb.get('target')
            samplefile_bn = os.path.basename(sinfo_fb.get('samplefile', basename))
            if not target:
                base_target = os.path.join(sinfo_fb['study'], samplefile_bn)
            else:
                base_target = target[:-1] if target.endswith('/') else target

        checksum_entries = []

        for sample in sorted(excluded_samples):
            sinfo = excluded_map.get(sample, {}).get('info')
            if not sinfo:
                continue

            remote_dir = os.path.join(base_target, "excluded", sample)

            def _add_file(path):
                if not path:
                    return
                lp = str(path)
                # If original was .bz2 and was converted to .gz on active, prefer .gz
                if lp.endswith('.bz2') and not os.path.exists(lp):
                    alt = lp[:-4] + '.gz'
                    if os.path.exists(alt):
                        lp = alt
                if not os.path.exists(lp):
                    print(f"* WARNING: Missing file for excluded sample {sample}: {lp}")
                    return
                remote_name = os.path.basename(lp)
                checksum_tmp = f"{output.checksums}.{sample}.{remote_name}.tmp"
                copy_with_checksum(lp, remote_dir, remote_name, checksum_tmp, AGH_DCACHE_CONFIG, params.ada_script)
                with open(checksum_tmp, 'r') as handle:
                    checksum_entries.append(f"{sample}/{remote_name}\t{handle.readline().strip()}\n")

            prefix = sinfo.get('prefix', '')
            ftype = sinfo.get('file_type')
            from_external = bool(sinfo.get('from_external'))
            data_base = pj(SOURCEDIR, f"{sample}.data")

            def _localize(f):
                if not f:
                    return None
                if os.path.isabs(f):
                    return f
                if from_external:
                    return pj(data_base, f)
                return read_samples.append_prefix(prefix, f)

            rgs = sinfo.get('readgroups', [])
            if rgs:
                for rg in rgs:
                    r_ftype = rg.get('file_type')
                    if r_ftype in ('fastq_paired', 'fastq'):
                        _add_file(_localize(rg.get('file1')))
                        _add_file(_localize(rg.get('file2')))
                    elif r_ftype == 'gvcf':
                        p = _localize(rg.get('file'))
                        _add_file(p)
                        if p:
                            idx = p + '.tbi'
                            if os.path.exists(idx):
                                _add_file(idx)
                    else:
                        p = _localize(rg.get('file'))
                        _add_file(p)
                        if p:
                            for idx in (p + '.bai', p + '.crai'):
                                if os.path.exists(idx):
                                    _add_file(idx)
            else:
                if ftype in ('fastq_paired', 'fastq'):
                    files1 = sinfo.get('file1', [])
                    files2 = sinfo.get('file2', [])
                    for f1, f2 in zip(files1, files2):
                        _add_file(_localize(f1))
                        _add_file(_localize(f2))
                elif ftype == 'gvcf':
                    files = sinfo.get('file1', [])
                    for f in files[:1]:
                        p = _localize(f)
                        _add_file(p)
                        if p:
                            idx = p + '.tbi'
                            if os.path.exists(idx):
                                _add_file(idx)
                else:
                    for f in sinfo.get('file1', []):
                        p = _localize(f)
                        _add_file(p)
                        if p:
                            for idx in (p + '.bai', p + '.crai'):
                                if os.path.exists(idx):
                                    _add_file(idx)

        with open(output.checksums, 'w') as sum_out:
            sum_out.writelines(checksum_entries)

        shell(f"touch {quote(str(output.copied))}")

def _samplefile_local_path(samplefile):
    return os.path.realpath(samplefile + '.tsv')

def _exclude_filename(samplefile_path):
    base = os.path.realpath(samplefile_path)
    if base.endswith('.tsv'):
        base = base[:-4]
    return base + '.exclude'

rule copy_samplefile_stats_to_dcache:
    input:
        bam=pj("{samplefile}.bam_quality.tab"),
        bam_rg=pj("{samplefile}.bam_rg_quality.tab"),
        oxo=pj("{samplefile}.oxo_quality.tab"),
        sex=pj("{samplefile}.sex_chrom.tab"),
        cov=pj("{samplefile}.coverage.hdf5"),
        kraken=pj("{samplefile}.kraken.tab"),
        phase=pj("{samplefile}.phase_quality.tab"),
        deepvariant=pj("{samplefile}.deepvariant_bcftools.tab"),
        samplefile=lambda wildcards: _samplefile_local_path(wildcards.samplefile)
    output:
        copied=temp(pj(STAT, "{samplefile}.stats_tabs.copied")),
        checksums=temp(pj(STAT, "{samplefile}.stats_tabs.ADLER32"))
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1"
        ,dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        target_dir = remote_base_for_samplefile(wildcards.samplefile)

        files = [
            (str(input.bam), os.path.basename(str(input.bam))),
            (str(input.bam_rg), os.path.basename(str(input.bam_rg))),
            (str(input.oxo), os.path.basename(str(input.oxo))),
            (str(input.sex), os.path.basename(str(input.sex))),
            (str(input.cov), os.path.basename(str(input.cov))),
            (str(input.kraken), os.path.basename(str(input.kraken))),
            (str(input.phase), os.path.basename(str(input.phase))),
            (str(input.deepvariant), os.path.basename(str(input.deepvariant))),
            (str(input.samplefile), os.path.basename(str(input.samplefile)))
        ]

        exclude_path = _exclude_filename(str(input.samplefile))
        if os.path.isfile(exclude_path):
            files.append((exclude_path, os.path.basename(exclude_path)))
        else:
            print(f"No exclusion file found for {wildcards.samplefile}, skipping upload")

        checksum_entries = []
        for local_file, remote_name in files:
            checksum_path = f"{output.checksums}.{remote_name}.tmp"
            copy_with_checksum(local_file, target_dir, remote_name, checksum_path, AGH_DCACHE_CONFIG, params.ada_script)
            with open(checksum_path, 'r') as handle:
                checksum_entries.append(f"{remote_name}\t{handle.readline().strip()}\n")

        with open(output.checksums, 'w') as sum_out:
            sum_out.writelines(checksum_entries)

        shell(f"touch {quote(str(output.copied))}")

rule tar_badmap_fastqs:
    input:
        fastq1=pj(FQ_BADMAP, "{sample}.badmap.R1.fastq.gz"),
        fastq2=pj(FQ_BADMAP, "{sample}.badmap.R2.fastq.gz")
    output:
        tar=temp(pj(FQ_BADMAP, "{sample}.badmap.fastqs.tar.gz"))
    resources:
        mem_mb=500,
        n="0.5"
    run:
        out_path = str(output.tar)
        out_dir = os.path.dirname(out_path)
        os.makedirs(out_dir, exist_ok=True)

        staging_path = f"{out_path}.tmp"
        file_map = [
            (str(input.fastq1), os.path.basename(str(input.fastq1))),
            (str(input.fastq2), os.path.basename(str(input.fastq2))),
        ]

        with tarfile.open(staging_path, "w:gz") as tar_handle:
            for local_path, arcname in file_map:
                tar_handle.add(local_path, arcname=arcname)

        os.replace(staging_path, out_path)

rule copy_badmap_to_dcache:
    input:
        tar=pj(FQ_BADMAP, "{sample}.badmap.fastqs.tar.gz")
    output:
        copied=temp(pj(FQ_BADMAP, "{sample}.badmap.tar.copied")),
        checksum=temp(pj(FQ_BADMAP, "{sample}.badmap.tar.ADLER32"))
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1"
        ,dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        base_target = remote_base_for_sample(wildcards.sample)
        remote_dir = os.path.join(base_target, "stat", "fq_badmap")

        local_tar = str(input.tar)
        remote_name = os.path.basename(local_tar)
        checksum_tmp = f"{output.checksum}.tmp"

        copy_with_checksum(local_tar, remote_dir, remote_name, checksum_tmp, AGH_DCACHE_CONFIG, params.ada_script)

        with open(checksum_tmp, 'r') as handle:
            checksum_value = handle.readline().strip()

        with open(output.checksum, 'w') as sum_out:
            sum_out.write(f"{remote_name}\t{checksum_value}\n")

        os.remove(checksum_tmp)

        shell(f"touch {quote(str(output.copied))}")

rule badmap_tar_all:
    input:
        expand(pj(FQ_BADMAP, "{sample}.badmap.tar.copied"), sample=sorted(sample_names))
    output:
        done=touch(pj(FQ_BADMAP, "badmap_uploads.done"))

rule stats_to_dcache_all:
    input:
        expand(pj(STAT, "{samplefile}.stats_bundle.tar.copied"), samplefile=SAMPLE_FILES),
        expand(pj(STAT, "{samplefile}.stats_tabs.copied"), samplefile=SAMPLE_FILES),
        expand(pj(FQ_BADMAP, "{sample}.badmap.tar.copied"), sample=sorted(sample_names))
    output:
        done=touch(pj(STAT, "stats_uploads.done"))

rule excluded_to_dcache_all:
    input:
        expand(pj(STAT, "{samplefile}.excluded.copied"), samplefile=SAMPLE_FILES)
    output:
        done=touch(pj(STAT, "excluded_uploads.done"))

def get_regions(wildcards):  #{{{
    samples = list(SAMPLEFILE_TO_SAMPLES[wildcards['samplefile']])
    samples.sort()
    return [pj(STAT,'cov','{sample}.regions.bed.gz'.format(sample=sample)) for sample in samples]  #}}}

def get_stats_samtools(wildcards):  #{{{
    samples = list(SAMPLEFILE_TO_SAMPLES[wildcards['samplefile']])
    samples.sort()
    return [pj(STAT,'{sample}.samtools.stat'.format(sample=sample)) for sample in samples]  #}}}

def get_samplefile_stat_tars(wildcards):  #{{{
    samples = list(SAMPLEFILE_TO_SAMPLES[wildcards['samplefile']])
    samples.sort()
    return [pj(STAT, f"{sample}.stats.tar.gz") for sample in samples]  #}}}

rule bundle_stats_tar_per_samplefile:
    input:
        get_samplefile_stat_tars
    output:
        tar=temp(pj(STAT, "{samplefile}.stats_bundle.tar.gz"))
    resources:
        n="0.5",
        mem_mb=200
    run:
        os.makedirs(os.path.dirname(output.tar), exist_ok=True)
        with tarfile.open(output.tar, "w:gz") as tar_handle:
            for path in input:
                print(path)
                tar_handle.add(path, arcname=os.path.basename(path))

rule write_samplefile_coverage_hdf5:
    input:
        bam=pj("{samplefile}.bam_quality.tab"),
        mapped=get_stats_samtools,
        cov=get_regions
    output:
        hdf5=pj('{samplefile}.coverage.hdf5')
    resources:
        n="1.0",
        mem_mb=14700
    run:
        samples = list(SAMPLEFILE_TO_SAMPLES[wildcards['samplefile']])
        samples.sort()
        annotation = WINDOWS_ANNOTATED

        read_stats.write_coverage_to_hdf5(annotation,samples,list(input.mapped),list(input.cov),output.hdf5)

def get_svd(wildcards):  #{{{
    """Returns the VerifyBamID SVD file for the sample type of the sample"""
    sinfo = SAMPLEINFO[wildcards['sample']]
    return VERIFYBAMID_WGS if 'wgs' in sinfo['sample_type'] else VERIFYBAMID_EXOME  #}}}

rule verifybamid:
    """Estimates contamination in a sample using the verifybamid2 tool"""
    input:
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        validated_sex=pj(KMER,"{sample}.result.yaml"),
    output:
        VBID_stat=temp(pj(STAT, 'contam/{sample}.verifybamid.pca2.selfSM')),
        VBID_ancestry=temp(pj(STAT, 'contam/{sample}.verifybamid.pca2.Ancestry'))
    # end of this file hardcoded in Haplotypecaller and read_contam_w
    priority: 27
    params:
        ref=get_ref_by_validated_sex,
        #minimum Base Quality for a base to contribute cov (default=20)
        VBID_prefix=pj(STAT,'contam/{sample}.verifybamid.pca2'),
        SVD=get_svd
    resources:
        mem_mb=300,
        n="1.4"
    conda: CONDA_VERIFYBAMID
    shell:
        """verifybamid2 --BamFile {input.bam} --SVDPrefix {params.SVD} --Reference {params.ref} --DisableSanityCheck --NumThread 2 --Output {params.VBID_prefix}"""


def get_capture_kit_interval_list(wildcards):  #{{{
    """Returns the capture kit interval list file for the sample type of the sample"""
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        capture_kit = ancient(MERGED_CAPTURE_KIT_IVL)
    else:
        if SAMPLEINFO[wildcards['sample']]['capture_kit']  == '':
            capture_kit = ancient(INTERSECT_CAPTURE_KIT_IVL)
        else: 
            capture_kit = pj(INTERVALS_DIR,SAMPLEINFO[wildcards['sample']]['capture_kit'] + '.interval_list')
    return capture_kit  #}}}

rule hs_stats:
    """Collects HS metrics for a sample using the gatk CollectHsMetrics tool"""
    input:
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        interval=ancient(MERGED_CAPTURE_KIT_IVL),
        validated_sex=pj(KMER,"{sample}.result.yaml"),
        targets=ancient(TARGETS_IVL),
    output:
        HS_metrics=(pj(STAT,"{sample}.hs_metrics"))
    priority: 99
    params:
        ref=get_ref_by_validated_sex,
        #minimum Base Quality for a base to contribute cov (default=20)
        Q=10,
        #minimum Mapping Quality for a read to contribute cov(default=20)
        MQ=10
    resources: mem_mb=lambda wildcards, attempt: attempt * 2000,
        tmpdir=tmpdir,
        n="1.0"
    conda: CONDA_VCF
    shell:
        """
            TMP_SSD="/scratch-node/${{USER}}.${{SLURM_JOB_ID}}"
            if [ ! -d "$TMP_SSD" ] || [ ! -w "$TMP_SSD" ]; then CAND=$(ls -1dt /scratch-node/${{USER}}.* 2>/dev/null | head -n1); if [ -n "$CAND" ] && [ -d "$CAND" ] && [ -w "$CAND" ]; then TMP_SSD="$CAND"; fi; fi
            if [ -d "$TMP_SSD" ] && [ -w "$TMP_SSD" ]; then TMPDIR_USE="$TMP_SSD"; elif [ -n "$SLURM_TMPDIR" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then TMPDIR_USE="$SLURM_TMPDIR"; else TMPDIR_USE="{resources.tmpdir}"; fi
            gatk  --java-options "-Xmx{resources.mem_mb}M  {DEFAULT_JAVA_OPTIONS}" CollectHsMetrics  --TMP_DIR "$TMPDIR_USE" \
                -I {input.bam} -R {params.ref} -BI {input.interval} -TI {input.targets} \
                -Q {params.Q} -MQ {params.MQ} \
                -O stats/{wildcards.sample}.hs_metrics
        """


rule artifacts_and_oxog_metrics:
    input:
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        interval=get_capture_kit_interval_list,
        validated_sex=pj(KMER,"{sample}.result.yaml")
    output:
        Bait_bias = temp(pj(STAT, '{sample}.bait_bias_summary_metrics')),
        Pre_adapter = temp(ensure(pj(STAT, '{sample}.pre_adapter_summary_metrics'),non_empty=True)),
        Bait_bias_det = temp(ensure(pj(STAT,'{sample}.bait_bias_detail_metrics'),non_empty=True)),
        Pre_adapter_det = temp(ensure(pj(STAT, '{sample}.pre_adapter_detail_metrics'),non_empty=True)),
        Error_summary = temp(ensure(pj(STAT, '{sample}.error_summary_metrics'),non_empty=True)),
        OXOG=temp(pj(STAT,"{sample}.OXOG"))
    priority: 99
    log: pj(LOG,"Stats","Artifact_OXOG_stats_{sample}.log")
    params:
        ref=get_ref_by_validated_sex,
        out=pj(STAT,"{sample}")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2600,
        tmpdir=tmpdir,
        n="1.0"
    conda: CONDA_VCF
    shell:
        """
            TMP_SSD="/scratch-node/${{USER}}.${{SLURM_JOB_ID}}"
            if [ ! -d "$TMP_SSD" ] || [ ! -w "$TMP_SSD" ]; then CAND=$(ls -1dt /scratch-node/${{USER}}.* 2>/dev/null | head -n1); if [ -n "$CAND" ] && [ -d "$CAND" ] && [ -w "$CAND" ]; then TMP_SSD="$CAND"; fi; fi
            if [ -d "$TMP_SSD" ] && [ -w "$TMP_SSD" ]; then TMPDIR_USE="$TMP_SSD"; elif [ -n "$SLURM_TMPDIR" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then TMPDIR_USE="$SLURM_TMPDIR"; else TMPDIR_USE="{resources.tmpdir}"; fi
            gatk --java-options "-Xmx{resources.mem_mb}M {DEFAULT_JAVA_OPTIONS}" CollectSequencingArtifactMetrics  --TMP_DIR "$TMPDIR_USE" -I {input.bam} -O {params.out} \
                    -R {params.ref} --DB_SNP {DBSNP} --INTERVALS {input.interval} 2> {log}
            gatk  --java-options "-Xmx{resources.mem_mb}M {DEFAULT_JAVA_OPTIONS}" CollectOxoGMetrics -I {input.bam} -O {output.OXOG} -R {params.ref} \
              --INTERVALS {input.interval} 2>> {log}
        """

# extract info about capture kit from SAMPLEFILE
# assume that all kits bed and interval_list files are existing and download to res folder
def get_capture_kit_bed(wildcards):  #{{{
    if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type']:
        capture_kit = MERGED_CAPTURE_KIT_BED
    elif SAMPLEINFO[wildcards['sample']]['capture_kit'] == '':
        capture_kit = MERGED_CAPTURE_KIT_BED
    else:
        capture_kit = SAMPLEINFO[wildcards['sample']]['capture_kit'] + '.bed'

    return pj(INTERVALS_DIR,capture_kit)  #}}}

rule samtools_stats:
    input:
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        validated_sex=pj(KMER,"{sample}.result.yaml")
    output:
        genome=temp(ensure(pj(STAT,"{sample}.samtools.stat"),non_empty=True)),
        exome=temp(ensure(pj(STAT,"{sample}.samtools.exome.stat"),non_empty=True))
    priority: 99
    log: pj(LOG,"Stats","samtools_{sample}.log")
    resources:
        mem_mb=130,
        n=1
    conda: CONDA_MAIN
    params:
        ref=get_ref_by_validated_sex,
        bed_interval=get_capture_kit_bed
    shell:
        """
        samtools stat -@ {resources.n} -r {params.ref} -d -p {input.bam} > {output.genome}
        samtools stat -@ {resources.n} -t {params.bed_interval} -d -p -r {params.ref} {input.bam} > {output.exome}
        """


    

rule bamstats_all_and_exome:
    input:
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
    output:
        all=temp(ensure(pj(STAT,'{sample}.bam_all.tsv'),non_empty=True)),
        exome=temp(ensure(pj(STAT,'{sample}.bam_exome.tsv'),non_empty=True))
    resources:
        mem_mb=250,
        n=1
    params:
        py_stats=srcdir(BAMSTATS),
        bed_interval=get_capture_kit_bed,
    conda: CONDA_PYPY
    shell:
        """
        samtools view -s 0.05 -h {input.bam} --threads {resources.n}  | pypy {params.py_stats} stats > {output.all}
        samtools view -s 0.05 -h {input.bam} --threads {resources.n} -L {params.bed_interval} | pypy {params.py_stats} stats > {output.exome}
        """

def get_quality_stats(wildcards):  #{{{
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    samples = sorted(sampleinfo.keys())
    required = []
    for sample in samples:
        required.extend(
            [
                pj(STAT, f"{sample}.samtools.stat"),
                pj(STAT, f"{sample}.samtools.exome.stat"),
                pj(STAT, 'contam', f"{sample}.verifybamid.pca2.selfSM"),
                pj(STAT, f"{sample}.bam_all.tsv"),
                pj(STAT, f"{sample}.bam_exome.tsv"),
                pj(STAT, f"{sample}.pre_adapter_summary_metrics"),
                pj(STAT, f"{sample}.bait_bias_summary_metrics"),
                pj(STAT, f"{sample}.hs_metrics"),
                pj(STAT, f"{sample}.chrM_read_stats.tsv"),
                pj(STAT, f"{sample}.numt_read_stats.tsv"),
            ]
        )

    return required

#}}}

rule gatherstats:
    # keep in mind, that samtools_stat create file even if it it's finished with error or you force to stop it
    # if you force to stop samtools_stat delete all output to prevent errors
    # rm -r stats/*samtools*
    input:
        get_quality_stats
    output:
        '{samplefile}.bam_quality.tab'
    resources:
        n="1.0",
        mem_mb=10000
    run:
        samples = sorted(SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])].keys())

        def paths_for(pattern):
            return [pj(STAT, pattern.format(sample=sample)) for sample in samples]

        stats = paths_for("{sample}.samtools.stat")
        exome_stats = paths_for("{sample}.samtools.exome.stat")
        vpca2 = [pj(STAT, 'contam', f"{sample}.verifybamid.pca2.selfSM") for sample in samples]
        bam_extra_all = paths_for("{sample}.bam_all.tsv")
        bam_extra_exome = paths_for("{sample}.bam_exome.tsv")
        pre_adapter = paths_for("{sample}.pre_adapter_summary_metrics")
        bait_bias = paths_for("{sample}.bait_bias_summary_metrics")
        hs_stats = paths_for("{sample}.hs_metrics")
        chrM_stats = paths_for("{sample}.chrM_read_stats.tsv")
        numt_stats = paths_for("{sample}.numt_read_stats.tsv")

        header, data = read_stats.combine_quality_stats(
            samples,
            stats,
            exome_stats,
            vpca2,
            bam_extra_all,
            bam_extra_exome,
            pre_adapter,
            bait_bias,
            hs_stats,
            chrM_stats,
            numt_stats,
        )
        read_stats.write_tsv(str(output),header,data)

def get_rg_quality_stats(wildcards):  #{{{
    smsinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    sample_readgroups = []
    samples = sorted(smsinfo.keys())    
    for i, sample in enumerate(samples):        
        x = sampleinfo(SAMPLEINFO, sample, checkpoint=True)
        rgs = x['readgroups']        
        for readgroup in rgs:
            sample_readgroups.append((sample, readgroup['info']['ID']))

    sample_readgroups.sort()

    required = []
    for sample, rg in sample_readgroups:
        required.append(pj(STAT, f"{sample}.{rg}.adapter_removal.log"))
        required.append(pj(STAT, f"{sample}.{rg}.fastq.adapters"))
        required.append(pj(STAT, f"{sample}.{rg}.merge_stats.tsv"))
        required.append(pj(STAT, f"{sample}.{rg}.dragmap.log"))
        required.append(pj(STAT, f"{sample}.{rg}.dechimer_stats.tsv"))
    
    return required

#}}}


def get_phase_stats(wildcards):  #{{{
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    return [pj(STAT, f"{sample}.phase_stats.tsv") for sample in sorted(sampleinfo.keys())]


#}}}


rule gather_phase_stats:
    input:
        get_phase_stats
    output:
        ensure(temp(pj("{samplefile}.phase_quality.tab")), non_empty=True)
    run:
        files = list(input)
        header = None
        rows = []
        for path in files:
            if not os.path.exists(path):
                continue
            with open(path, 'r') as handle:
                reader = csv.reader(handle, delimiter='\t')
                hdr = next(reader, None)
                rec = next(reader, None)
                if hdr is None or rec is None:
                    continue
                if header is None:
                    header = hdr
                rows.append(rec)
        if header is None:
            header = [
                'sample','sex',
                'phase_auto_variants','phase_auto_het','phase_auto_phased','phase_auto_het_phased_pct','phase_auto_blocks','phase_auto_block_sizes_sum_variants','phase_auto_block_lengths_sum_bp','phase_auto_block_longest_bp','phase_auto_block_shortest_bp',
                'phase_chrx_variants','phase_chrx_het','phase_chrx_phased','phase_chrx_het_phased_pct','phase_chrx_blocks','phase_chrx_block_sizes_sum_variants','phase_chrx_block_lengths_sum_bp','phase_chrx_block_longest_bp','phase_chrx_block_shortest_bp',
                'phase_chrm_variants','phase_chrm_het','phase_chrm_phased','phase_chrm_het_phased_pct','phase_chrm_blocks','phase_chrm_block_sizes_sum_variants','phase_chrm_block_lengths_sum_bp','phase_chrm_block_longest_bp','phase_chrm_block_shortest_bp',
                'phase_all_variants','phase_all_het','phase_all_phased','phase_all_het_phased_pct','phase_all_blocks','phase_all_block_sizes_sum_variants','phase_all_block_lengths_sum_bp','phase_all_block_longest_bp','phase_all_block_shortest_bp',
            ]
        with open(output[0], 'w') as out:
            w = csv.writer(out, delimiter='\t')
            w.writerow(header)
            for r in rows:
                w.writerow(r)

rule gather_rg_stats:
    input:
        get_rg_quality_stats
    output:
        '{samplefile}.bam_rg_quality.tab'
    resources:
        n="1",
        mem_mb=1000
    run:
        smsinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
        sample_readgroups = []
        samples = sorted(smsinfo.keys())
        for i, sample in enumerate(samples):
            x = sampleinfo(SAMPLEINFO, sample, checkpoint=True)
            rgs = x['readgroups']
            for readgroup in rgs:
                sample_readgroups.append((sample, readgroup['info']['ID']))

        sample_readgroups.sort()

        aremoval = [pj(STAT, f"{sample}.{rg}.adapter_removal.log") for sample, rg in sample_readgroups]
        aidentify = [pj(STAT, f"{sample}.{rg}.fastq.adapters") for sample, rg in sample_readgroups]
        mergestats = [pj(STAT, f"{sample}.{rg}.merge_stats.tsv") for sample, rg in sample_readgroups]
        dragmap_stats = [pj(STAT, f"{sample}.{rg}.dragmap.log") for sample, rg in sample_readgroups]
        dechimer_stats = [pj(STAT, f"{sample}.{rg}.dechimer_stats.tsv") for sample, rg in sample_readgroups]
        
        missing_inputs = [p for p in list(input) if not os.path.exists(str(p))]
        if missing_inputs:
            raise FileNotFoundError(f"gather_rg_stats missing {len(missing_inputs)} inputs; first_missing={missing_inputs[:5]}")
        
        header, data = read_stats.combine_rg_quality_stats(sample_readgroups,aremoval,aidentify,mergestats,dragmap_stats,dechimer_stats)        
        out_path = str(output[0])        
        read_stats.write_tsv(out_path,header,data)        

def get_oxo_stats(wildcards):  #{{{
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    samples = sorted(sampleinfo.keys())
    required = []
    for sample in samples:
        required.append(pj(STAT, f"{sample}.pre_adapter_detail_metrics"))
        required.append(pj(STAT, f"{sample}.bait_bias_detail_metrics"))
    return required  #}}}


def get_deepvariant_bcftools_summaries(wildcards):  #{{{
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    samples = sorted(sampleinfo.keys())
    outputs = []
    for sample in samples:
        if 'wgs' in SAMPLEINFO[sample]['sample_type']:
            regions = level1_regions
        else:
            regions = level0_regions
        for region in regions:
            outputs.append(pj(STAT, "deepvariant_bcftools", f"{sample}.{region}.summary.tsv"))
    return outputs  #}}}

def get_deepvariant_gvcf_outputs(wildcards):  #{{{
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    samples = sorted(sampleinfo.keys())
    outputs = []
    for sample in samples:
        if 'wgs' in SAMPLEINFO[sample]['sample_type']:
            for region in level1_regions:
                outputs.append(pj(DEEPVARIANT, 'gVCF', 'exome_extract', region, f'{sample}.{region}.wg.vcf.gz'))
        else:
            for region in level0_regions:
                outputs.append(pj(DEEPVARIANT, 'gVCF', region, f'{sample}.{region}.wg.vcf.gz'))
    return outputs  #}}}


def get_kraken_stats(wildcards):  #{{{
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    return [pj(KRAKEN, f"{sample}.kraken_summary.tsv") for sample in sampleinfo.keys()]


#}}}


rule gather_deepvariant_bcftools_stats:
    input:
        summaries=get_deepvariant_bcftools_summaries,
        gvcfs=get_deepvariant_gvcf_outputs
    output:
        ensure(pj("{samplefile}.deepvariant_bcftools.tab"), non_empty=True)
    run:
        summaries = sorted([str(p) for p in input.summaries])
        aggregated = {}
        field_order = []

        def remember_field(field):
            if field not in field_order:
                field_order.append(field)

        def parse_number(value):
            if value is None:
                return None
            text = str(value).strip()
            if text == '' or text.lower() == 'nan':
                return None
            try:
                if any(ch in text for ch in ['.', 'e', 'E']):
                    return float(text)
                return int(text)
            except ValueError:
                return None

        def region_to_slice(region):
            if not region:
                return ''
            if region.startswith('A'):
                return 'A'
            if region.startswith('F'):
                return 'F'
            if region.startswith('X'):
                return 'XH' if region.endswith('H') else 'X'
            if region.startswith('Y'):
                return 'YH' if region.endswith('H') else 'Y'
            return region

        for path in summaries:
            if not os.path.exists(path):
                continue
            basename = os.path.basename(path)
            parts = basename.split('.')
            sample_name = parts[0] if parts else ''
            region_name = parts[1] if len(parts) > 1 else ''
            with open(path, "r") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fieldnames = reader.fieldnames or []
                for field in fieldnames:
                    if field in {"sample", "region"}:
                        continue
                    remember_field(field)
                for record in reader:
                    sample = record.get('sample', sample_name)
                    # Determine slice from filename (record doesn't have region field)
                    region_tok = record.get('region', region_name)
                    slice_name = region_to_slice(region_tok)
                    agg = aggregated.setdefault((sample, slice_name), {"sums": {}, "strings": {}})
                    for field, value in record.items():
                        if field in {"sample", "region"}:
                            continue
                        remember_field(field)
                        numeric = parse_number(value)
                        if numeric is None:
                            if value not in ('', None):
                                agg["strings"].setdefault(field, value)
                            continue
                        if field == 'number_of_samples':
                            prev = agg["sums"].get(field)
                            agg["sums"][field] = numeric if prev is None else max(prev, numeric)
                        else:
                            agg["sums"][field] = agg["sums"].get(field, 0.0) + numeric

        derived_fields = [
            "tstv_ratio",
            "tstv_ratio_first_alt",
            "indel_repeat_ratio",
            "indel_repeat_consistent_frac",
            "indel_repeat_inconsistent_frac",
            "avg_length_consistent_del",
            "avg_length_inconsistent_del",
            "avg_length_consistent_ins",
            "avg_length_inconsistent_ins",
        ]

        def format_number(value):
            if isinstance(value, float):
                if math.isfinite(value) and abs(value - round(value)) < 1e-9:
                    return str(int(round(value)))
                return f"{value:.6f}"
            return str(value)

        def safe_ratio(numerator, denominator):
            if denominator in (None, 0) or numerator is None:
                return "NaN"
            return f"{numerator / denominator:.6f}"

        def safe_avg(total, count):
            if count in (None, 0) or total is None:
                return "NaN"
            return f"{total / count:.6f}"

        header = ["sample", "slice"] + field_order + derived_fields

        if not aggregated:
            with open(output[0], 'w') as out:
                writer = csv.writer(out, delimiter='\t')
                writer.writerow(header)
            return

        with open(output[0], 'w') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(header)
            for sample, slice_name in sorted(aggregated):
                sums = aggregated[(sample, slice_name)]["sums"]
                strings = aggregated[(sample, slice_name)]["strings"]

                total_indels = sums.get("number_of_indels")
                consistent = sums.get("indel_repeat_consistent")
                inconsistent = sums.get("indel_repeat_inconsistent")

                derived = {
                    "tstv_ratio": safe_ratio(sums.get("tstv_ts"), sums.get("tstv_tv")),
                    "tstv_ratio_first_alt": safe_ratio(sums.get("tstv_ts_first_alt"), sums.get("tstv_tv_first_alt")),
                    "indel_repeat_ratio": safe_ratio(consistent, (consistent or 0) + (inconsistent or 0)),
                    "indel_repeat_consistent_frac": safe_ratio(consistent, total_indels),
                    "indel_repeat_inconsistent_frac": safe_ratio(inconsistent, total_indels),
                }

                for key, avg_key in [
                    ("sum_length_consistent_del", "avg_length_consistent_del"),
                    ("sum_length_inconsistent_del", "avg_length_inconsistent_del"),
                    ("sum_length_consistent_ins", "avg_length_consistent_ins"),
                    ("sum_length_inconsistent_ins", "avg_length_inconsistent_ins"),
                ]:
                    total = sums.get(key)
                    count = sums.get(f"{key}_count")
                    derived[avg_key] = safe_avg(total, count)

                row = [sample, slice_name]
                for field in field_order:
                    if field in sums:
                        row.append(format_number(sums[field]))
                    elif field in strings:
                        row.append(strings[field])
                    else:
                        row.append('')
                for field in derived_fields:
                    row.append(derived.get(field, "NaN"))
                writer.writerow(row)


rule gatherkraken:
    input:
        get_kraken_stats
    output:
        '{samplefile}.kraken.tab'
    resources:
        n="1",
        mem_mb=500
    params:
        merge_script=srcdir("scripts/kraken_merge.py"),
    run:
        files = [str(p) for p in sorted(list(input))]
        manifest = f"{output[0]}.summaries.txt"
        with open(manifest, 'w') as handle:
            handle.write("\n".join(files) + "\n")
        shell(f"python {params.merge_script} --manifest {quote(manifest)} --output {quote(str(output[0]))}")
        os.remove(manifest)


rule gatheroxostats:
    input:
        get_oxo_stats
    output:
        '{samplefile}.oxo_quality.tab'
    resources:
        n="1",
        mem_mb=1000
    run:
        samples = sorted(SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])].keys())
        pre_adapter = [pj(STAT, f"{sample}.pre_adapter_detail_metrics") for sample in samples]
        bait_bias = [pj(STAT, f"{sample}.bait_bias_detail_metrics") for sample in samples]

        header, data = read_stats.combine_oxo_stats(samples, pre_adapter, bait_bias)
        read_stats.write_tsv(str(output),header,data)


def get_sex_stats(wildcards):  #{{{
    sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
    return [pj(KMER,f'{sample}.result.yaml') for sample in sampleinfo.keys()]


#}}}

rule gathersexstats:
    input:
        get_sex_stats
    output:
        '{samplefile}.sex_chrom.tab'
    resources:
        n="1",
        mem_mb=1000
    run:
        sampleinfo = SAMPLEFILE_TO_SAMPLES[os.path.basename(wildcards['samplefile'])]
        samples = list(sampleinfo.keys())
        samples.sort()

        kmer_stats = [pj(KMER,f"{sample}.result.yaml") for sample in samples]
        sex_reported = [sampleinfo[sample]['sex'] for sample in samples]

        header, data = read_stats.combine_sex_stats(samples,kmer_stats,sex_reported)
        read_stats.write_tsv(str(output),header,data)



rule mospeth_mergedCK:
    input:
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        interval=MERGED_CAPTURE_KIT_BED
    output:         pj(STAT,'cov','{sample}_MERGED_CK.regions.bed.gz'),
        pj(STAT,'cov','{sample}_MERGED_CK.regions.bed.gz.csi'),
        temp(pj(STAT,'cov','{sample}_MERGED_CK.mosdepth.global.dist.txt')),
        temp(pj(STAT,'cov','{sample}_MERGED_CK.mosdepth.summary.txt'))
    params:
        prefix=pj(STAT,'cov','{sample}_MERGED_CK')
    resources:
        mem_mb=2200,
        n="1.8"
    conda: CONDA_MOSDEPTH
    shell:
        """
            mkdir -p `dirname {output[0]}`
            mosdepth  --threads 2 -b {input.interval} --no-per-base {params.prefix} {input.bam}
        """


rule get_capture_kit:
    input: 
        bam = pj(BAM,"{sample}.markdup.bam"),
        target = MERGED_CAPTURE_KIT_BED,
    output: 
        capture_kit_stats = pj(STAT,"{sample}.capture_kit_stats.tsv"),
        cov_decompressed = temp(pj(STAT,"cov","{sample}.regions.bed")),
    params: 
        expected_kit = lambda wildcards: SAMPLEINFO[wildcards['sample']]['capture_kit'],
        CAPTURE_KIT_CHECKER = srcdir(CAPTURE_KIT_CHECKER)
    conda: CONDA_CK_FINDER
    resources: 
        n = "16",
        mem_mb = 8000
    shell:
        """
        set -ex
        mosdepth --threads {resources.n} -n --by {input.target} {wildcards.sample} {input.bam}
        mkdir -p $(dirname {output.cov_decompressed})
        gunzip -c {wildcards.sample}.regions.bed.gz > {output.cov_decompressed}
        python {params.CAPTURE_KIT_CHECKER} \
            --coverage {output.cov_decompressed} \
            --metadata_capture {params.expected_kit} \
            --output {output.capture_kit_stats}
        rm {wildcards.sample}.regions.bed.gz {wildcards.sample}.mosdepth.summary.txt
        """
