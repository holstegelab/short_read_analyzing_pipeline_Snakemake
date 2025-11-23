wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    region = r"[\w\d]+",
    # readgroup="[\w\d_\-@]+"
onsuccess: shell("rm -fr logs/Deepvariant/*")
from common import *
from shlex import quote
import tarfile
import tempfile
import shutil
from pathlib import Path


module Tools:
    snakefile: 'Tools.smk'
    config: config
use rule * from Tools
mode = config.get("computing_mode", "WES")


def level2_parent_level1(region):
    """Return the parent level1 region for a given level2 region."""
    return convert_to_level1(region)


def level2_interval(region, wgs):
    """Return the padded interval list for the level2 region, using WGS or WES bins."""
    return region_to_file(region, wgs=wgs, extension="bed")


def deepvariant_level2_inputs_wgs(wildcards):
    """Collect per-sample WGS gVCF paths for a given level2 region and samplefile."""
    region = wildcards.region
    samplefile = wildcards.samplefile
    parent_level1 = level2_parent_level1(region)
    files = []
    for sample, sinfo in SAMPLEFILE_TO_SAMPLES[samplefile].items():
        if 'wgs' not in sinfo['sample_type']:
            continue
        base_path = pj(DEEPVARIANT, "gVCF", parent_level1, f"{sample}.{parent_level1}.wg.vcf.gz")
        files.append(base_path)
        files.append(base_path + ".tbi")
    if not files:
        raise ValueError(f"No WGS DeepVariant gVCFs found for samplefile {samplefile} in region {region}")
    return files


def deepvariant_level2_inputs_wes(wildcards):
    """Collect per-sample exome gVCF paths for a given level2 region and samplefile."""
    region = wildcards.region
    samplefile = wildcards.samplefile
    parent_level1 = level2_parent_level1(region)
    files = []
    for sample, sinfo in SAMPLEFILE_TO_SAMPLES[samplefile].items():
        if 'wgs' in sinfo['sample_type']:
            base_path = pj(DEEPVARIANT, "gVCF", "exome_extract", parent_level1, f"{sample}.{parent_level1}.wg.vcf.gz")
        else:
            base0 = convert_to_level0(parent_level1)
            base_path = pj(DEEPVARIANT, "gVCF", base0, f"{sample}.{base0}.wg.vcf.gz")
        files.append(base_path)
        files.append(base_path + ".tbi")
    return files


def deepvariant_level2_wgs_samples(wildcards):
    return [
        sample
        for sample, sinfo in SAMPLEFILE_TO_SAMPLES[wildcards.samplefile].items()
        if 'wgs' in sinfo['sample_type']
    ]


def deepvariant_level2_samples(wildcards):
    return list(SAMPLEFILE_TO_SAMPLES[wildcards.samplefile].keys())


rule extract_and_tar_deepvariant_level2_wgs:
    input:
        gvcfs=deepvariant_level2_inputs_wgs
    output:
        tar=temp(pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar"))
    params:
        samples=deepvariant_level2_wgs_samples,
        region=lambda wc: wc.region,
        samplefile=lambda wc: wc.samplefile
    resources:
        mem_mb=4000,
        n="1.5"
    run:
        region = wildcards.region
        samplefile = wildcards.samplefile
        samples = params.samples
        tmpdir = Path(tempfile.mkdtemp(prefix=f"dv_l2_{samplefile}_{region}_"))
        try:
            Path(output.tar).parent.mkdir(parents=True, exist_ok=True)
            staged_paths = []
            parent_level1 = level2_parent_level1(region)
            interval = level2_interval(region, wgs=True)
            for sample in samples:
                source_path = Path(pj(DEEPVARIANT, "gVCF", parent_level1, f"{sample}.{parent_level1}.wg.vcf.gz"))
                if not source_path.exists():
                    raise FileNotFoundError(f"Source gVCF missing: {source_path}")
                staged_gvcf = tmpdir / f"{sample}.{region}.dv.wgs.g.vcf.gz"
                staged_tbi = tmpdir / f"{sample}.{region}.dv.wgs.g.vcf.gz.tbi"
                shell(
                    "bcftools view -R {interval:q} {source_path:q} -O z -o {staged_gvcf:q}"
                )
                shell(
                    "bcftools index -f -t {staged_gvcf:q}"
                )
                if not staged_gvcf.exists() or not staged_tbi.exists():
                    raise FileNotFoundError(f"Extraction failed for sample {sample} region {region}")
                staged_paths.extend([staged_gvcf, staged_tbi])
            if not staged_paths:
                raise ValueError(f"No staged gVCFs for {samplefile} {region}")
            with tarfile.open(output.tar, "w") as tar_handle:
                for path in staged_paths:
                    tar_handle.add(path, arcname=path.name)
            with tarfile.open(output.tar, "r") as tar_handle:
                members = tar_handle.getnames()
            expected = [f"{sample}.{region}.dv.wgs.g.vcf.gz" for sample in samples]
            expected += [f"{sample}.{region}.dv.wgs.g.vcf.gz.tbi" for sample in samples]
            missing = sorted(set(expected) - set(members))
            if missing:
                raise ValueError(f"Tarball missing entries: {missing}")
        finally:
            shutil.rmtree(tmpdir)


rule extract_and_tar_deepvariant_level2_wes:
    input:
        gvcfs=deepvariant_level2_inputs_wes
    output:
        tar=temp(pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar"))
    params:
        samples=deepvariant_level2_samples,
        region=lambda wc: wc.region,
        samplefile=lambda wc: wc.samplefile
    resources:
        mem_mb=4000,
        n="1.5"
    run:
        region = wildcards.region
        samplefile = wildcards.samplefile
        samples = params.samples
        tmpdir = Path(tempfile.mkdtemp(prefix=f"dv_l2_wes_{samplefile}_{region}_"))
        try:
            Path(output.tar).parent.mkdir(parents=True, exist_ok=True)
            staged_paths = []
            parent_level1 = level2_parent_level1(region)
            interval = level2_interval(region, wgs=False)
            for sample in samples:
                sinfo = SAMPLEFILE_TO_SAMPLES[samplefile][sample]
                if 'wgs' in sinfo['sample_type']:
                    source_path = Path(pj(DEEPVARIANT, "gVCF", "exome_extract", parent_level1, f"{sample}.{parent_level1}.wg.vcf.gz"))
                else:
                    source_region = convert_to_level0(parent_level1)
                    source_path = Path(pj(DEEPVARIANT, "gVCF", source_region, f"{sample}.{source_region}.wg.vcf.gz"))
                if not source_path.exists():
                    raise FileNotFoundError(f"Exome gVCF missing: {source_path}")
                staged_gvcf = tmpdir / f"{sample}.{region}.dv.wes.g.vcf.gz"
                staged_tbi = tmpdir / f"{sample}.{region}.dv.wes.g.vcf.gz.tbi"
                shell(
                    "bcftools view -R {interval:q} {source_path:q} -O z -o {staged_gvcf:q}"
                )
                shell(
                    "bcftools index -f -t {staged_gvcf:q}"
                )
                if not staged_gvcf.exists() or not staged_tbi.exists():
                    raise FileNotFoundError(f"Exome extraction failed for sample {sample} region {region}")
                staged_paths.extend([staged_gvcf, staged_tbi])
            if not staged_paths:
                raise ValueError(f"No staged exome gVCFs for {samplefile} {region}")
            with tarfile.open(output.tar, "w") as tar_handle:
                for path in staged_paths:
                    tar_handle.add(path, arcname=path.name)
            with tarfile.open(output.tar, "r") as tar_handle:
                members = tar_handle.getnames()
            expected = [f"{sample}.{region}.dv.wes.g.vcf.gz" for sample in samples]
            expected += [f"{sample}.{region}.dv.wes.g.vcf.gz.tbi" for sample in samples]
            missing = sorted(set(expected) - set(members))
            if missing:
                raise ValueError(f"WES tarball missing entries: {missing}")
        finally:
            shutil.rmtree(tmpdir)


rule DeepVariant_all:
    input:
        expand("{dv}/{sample}.done",sample=sample_names, dv = DEEPVARIANT),


WGS_SAMPLEFILES = [
    samplefile
    for samplefile in SAMPLE_FILES
    if any('wgs' in sinfo['sample_type'] for sinfo in SAMPLEFILE_TO_SAMPLES[samplefile].values())
]


def deepvariant_region_inputs(wildcards):
    region = wildcards.region
    chunk = convert_to_level1(region)
    files = []
    for sample, sinfo in SAMPLEFILE_TO_SAMPLES[wildcards.samplefile].items():
        if 'wgs' not in sinfo['sample_type']:
            continue
        base = pj(DEEPVARIANT, "gVCF", chunk, f"{sample}.{chunk}.wg.vcf.gz")
        files.append(base)
        files.append(base + ".tbi")
    if not files:
        raise ValueError(f"No DeepVariant gVCFs found for samplefile {wildcards.samplefile} in region {region}")
    return files


rule copy_deepvariant_wgs_region_to_dcache:
    input:
        tar=pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar")
    output:
        copied=temp(pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar.copied")),
        checksum=temp(pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar.ADLER32"))
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1",
        dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        remote_dir = os.path.join(remote_base_for_samplefile(wildcards.samplefile), "gvcf", "deepvariant", "level2", "wgs")
        remote_name = os.path.basename(input.tar)
        copy_with_checksum(str(input.tar), remote_dir, remote_name, str(output.checksum), AGH_DCACHE_CONFIG, params.ada_script)
        shell(f"touch {quote(str(output.copied))}")


rule copy_deepvariant_wes_region_to_dcache:
    input:
        tar=pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar")
    output:
        copied=pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar.copied"),
        checksum=pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar.ADLER32")
    params:
        ada_script=srcdir(ADA)
    resources:
        mem_mb=2000,
        n="0.1",
        dcache_use_add=config.get('dcache_use_add', 0),
        dcache_use_remove=config.get('dcache_use_remove', 0)
    run:
        remote_dir = os.path.join(remote_base_for_samplefile(wildcards.samplefile), "gvcf", "deepvariant", "level2", "wes")
        remote_name = os.path.basename(input.tar)
        copy_with_checksum(str(input.tar), remote_dir, remote_name, str(output.checksum), AGH_DCACHE_CONFIG, params.ada_script)
        shell(f"touch {quote(str(output.copied))}")



rule deepvariant_tar_wgs_all:
    input:
        expand(pj(GVCF_TAR, "deepvariant_level2_wgs", "{samplefile}.{region}.dv.wgs.gvcf.tar.copied"), samplefile=WGS_SAMPLEFILES, region=level2_regions)
    output:
        done=temp(touch(pj(GVCF_TAR, "deepvariant_gvcf_wgs_uploads.done")))


rule deepvariant_tar_wes_all:
    input:
        expand(pj(GVCF_TAR, "deepvariant_level2_wes", "{samplefile}.{region}.dv.wes.gvcf.tar.copied"), samplefile=SAMPLE_FILES, region=level2_regions)
    output:
        done=temp(touch(pj(GVCF_TAR, "deepvariant_gvcf_wes_uploads.done")))


def get_deepvariant_files(wildcards):#{{{
    sample = wildcards['sample']
    if 'wgs' in SAMPLEINFO[sample]['sample_type']:
        return [pj(DEEPVARIANT,  'gVCF', 'exome_extract', region, f'{sample}.{region}.wg.vcf.gz') for region in level1_regions]
    else:
        return [pj(DEEPVARIANT,  'gVCF', region, f'{sample}.{region}.wg.vcf.gz') for region in level0_regions]
#}}}


rule deepvariant_sample_done:
    input:
        get_deepvariant_files
    output:
        done=temp(touch(pj(DEEPVARIANT, "{sample}.done")))
    resources:
        mem_mb = 100,
        n = "1.0"

def get_sequencing_mode(wildcards):#{{{
    return "WGS" if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else "WES"#}}}

def get_mem_mb_deepvariant(wildcards, attempt):#{{{
    res = 10000
    return (attempt - 1) * 0.5 * res + res#}}}


def region_to_bed_file(wildcards):#{{{
    """Converts a region to a bed file location (see common.py and Tools.smk)"""
    sample = wildcards['sample']
    region = wildcards['region']
    return region_to_file(region, wgs='wgs' in SAMPLEINFO[sample]['sample_type'], extension='bed')#}}}

def region_to_bed_file_wgs(wildcards):#{{{
    region = wildcards['region']
    return region_to_file(region, wgs=True, extension='bed')#}}}

rule deepvariant:
    input:
        bed = region_to_bed_file,
        bed_wgs = region_to_bed_file_wgs,
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        validated_sex=pj(KMER,"{sample}.result.yaml"),
    output:
        vcf = (temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz"))),
        vcf_tbi = (temp(pj(DEEPVARIANT,'VCF', "{region}","{sample}.{region}.vcf.gz.tbi"))),
        gvcf = (temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz"))),
        gvcf_tbi = (temp(pj(DEEPVARIANT,'gVCF', "{region}","{sample}.{region}.g.vcf.gz.tbi")))
    params:
            mode=get_sequencing_mode,
            haploid_contigs=lambda wildcards: 'chrX,chrX_KI270880v1_alt,chrX_KI270881v1_alt,chrX_KI270913v1_alt,chrY,chrY_KI270740v1_random' if wildcards['region'].endswith("H") else 'chrNONE',
            skipsex = lambda wildcards, input: int(get_validated_sex_file(input) == 'female' and wildcards['region'].startswith('Y')),
            inter_dir = pj(DEEPVARIANT,'DV_intermediate'),
            # check = CHECKEMPTY
    container: 'docker://google/deepvariant:1.9.0'
    resources:
        n="7",
        nshards=8,
        mem_mb=get_mem_mb_deepvariant,
        time = 6600
    shell:
        """
        if [ {params.skipsex} -eq 0 ]
        then
            TMP_SSD="/scratch-node/${{USER}}.${{SLURM_JOB_ID}}"
            JOB_ID="${{SLURM_JOB_ID}}"
            if [ -z "$JOB_ID" ]; then JOB_ID="${{SLURM_JOBID}}"; fi
            if [ -z "$JOB_ID" ]; then JOB_ID="$$"; fi
            if [ ! -d "$TMP_SSD" ] || [ ! -w "$TMP_SSD" ]; then CAND=$(ls -1dt /scratch-node/${{USER}}.* 2>/dev/null | head -n1); if [ -n "$CAND" ] && [ -d "$CAND" ] && [ -w "$CAND" ]; then TMP_SSD="$CAND"; fi; fi
            if [ -d "$TMP_SSD" ] && [ -w "$TMP_SSD" ]; then RUNDIR_BASE="$TMP_SSD/deepvariant/$JOB_ID"; elif [ -n "$SLURM_TMPDIR" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then RUNDIR_BASE="$SLURM_TMPDIR/deepvariant/$JOB_ID"; elif [ -d "/tmp" ] && [ -w "/tmp" ]; then RUNDIR_BASE="/tmp/${{USER}}/deepvariant/$JOB_ID"; else RUNDIR_BASE="{params.inter_dir}"; fi
            RUNDIR="$RUNDIR_BASE/{wildcards.sample}.{wildcards.region}"
            echo "SSD base: $TMP_SSD" >&2
            echo "RUNDIR_BASE: $RUNDIR_BASE" >&2
            echo "JOB_ID: $JOB_ID" >&2
            echo "RUNDIR: $RUNDIR" >&2
            /bin/rm -rf "$RUNDIR" 2>/dev/null || true
            mkdir -p "$RUNDIR"
            trap '/bin/rm -rf "$RUNDIR" 2>/dev/null || true' EXIT INT TERM

            OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
            TF_NUM_INTRAOP_THREADS={resources.nshards} TF_NUM_INTEROP_THREADS={resources.nshards} \
            /opt/deepvariant/bin/run_deepvariant \
              --make_examples_extra_args "normalize_reads=true,regions={input.bed},small_model_call_multiallelics=false" \
              --call_variants_extra_args "config_string=inter_op_parallelism_threads: {resources.nshards} intra_op_parallelism_threads: {resources.nshards} device_count: {{ key: 'CPU' value: {resources.nshards} }}" \
              --num_shards={resources.nshards} \
              --model_type={params.mode} \
              --ref={REF_MALE} --reads={input.bam} \
              --output_vcf={output.vcf} --output_gvcf={output.gvcf} \
              --haploid_contigs {params.haploid_contigs} \
              --intermediate_results_dir "$RUNDIR" \
              --postprocess_cpus {resources.nshards}
        else
            printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" | bgzip -c > {output.vcf}
            tabix -f -p vcf {output.vcf}
            printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" | bgzip -c > {output.gvcf}
            tabix -f -p vcf {output.gvcf}            
        fi
        """

# python {params.check} {output.vcf}
# python {params.check} {output.gvcf}

rule DVWhatshapPhasingMerge:
    """Phase VCF with Whatshap and merge into the gVCF"""
    input:
        vcf = rules.deepvariant.output.vcf,
        vcf_tbi = rules.deepvariant.output.vcf_tbi,
        gvcf = rules.deepvariant.output.gvcf,
        gvcf_tbi = rules.deepvariant.output.gvcf_tbi,
        bams=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        validated_sex=pj(KMER,"{sample}.result.yaml"),
    output:
        vcf = temp(pj(DEEPVARIANT, "VCF/{region}/{sample}.{region}.w.vcf.gz")),
        vcf_tbi = temp(pj(DEEPVARIANT, "VCF/{region}/{sample}.{region}.w.vcf.gz.tbi")),
        wstats = pj(STAT, "whatshap_dvphasing/{sample}.{region}.stats"),
        mwstats = pj(STAT, "whatshap_dvphasing/{sample}.{region}.merge_stats"),
        bcftools_stats = temp(pj(STAT, "deepvariant_bcftools/{sample}.{region}.bcftools_stats.txt")),
        bcftools_summary = ensure(temp(pj(STAT, "deepvariant_bcftools/{sample}.{region}.summary.tsv")), non_empty=True),
        tmp_gvcf= temp(pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf")),
        gvcf= temp(pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz")),
        gvcf_tbi = temp(pj(DEEPVARIANT, "gVCF/{region}/{sample}.{region}.wg.vcf.gz.tbi")),
        gvcf_exome = ensure(temp(pj(DEEPVARIANT, "gVCF/exome_extract/{region}/{sample}.{region}.wg.vcf.gz")), non_empty = True),
        gvcf_exome_tbi = ensure(temp(pj(DEEPVARIANT, "gVCF/exome_extract/{region}/{sample}.{region}.wg.vcf.gz.tbi")), non_empty = True),
    params:
        merge_script=srcdir(MERGEPHASEDIRECT),
        stats_parser=srcdir("scripts/deepvariant_bcftools_stats_parser.py"),
        ploidy=lambda wildcards: 1 if wildcards["region"].endswith("H") else 2,
        skipsex = lambda wildcards, input: int(get_validated_sex_file(input) == 'female' and wildcards['region'].startswith('Y')),
        interval_bed=lambda wildcards: region_to_file(region=wildcards.region, extension="bed", padding=True),
        capture_intersect_bed=lambda wildcards: (
            INTERSECT_CAPTURE_KIT_AUTO_BED if wildcards.region.startswith('A') or wildcards.region.startswith('F') else (
            INTERSECT_CAPTURE_KIT_X_BED if wildcards.region.startswith('X') else (
            INTERSECT_CAPTURE_KIT_Y_BED))
        )
    log: pj(LOG, "Deepvariant", "{sample}.{region}.whatshap.log"),
    resources: 
        n="1.0",
        mem_mb = 1500
    conda: CONDA_VCF
    shell: """
        mkdir -p `dirname {output.wstats}` `dirname {output.vcf}` `dirname {output.gvcf}`
        if [ {params.ploidy} -eq 2 ] && [ {params.skipsex} -eq 0 ]
        then 
            if [ {params.skipsex} -eq 0 ]
            then 
                echo "[DEBUG] start whatshap phase: {wildcards.sample} {wildcards.region}" >&2
                echo "[DEBUG] input.vcf: {input.vcf}" >&2
                echo "[DEBUG] output.vcf: {output.vcf}" >&2
                ls -l `dirname {output.vcf}` || true
                whatshap phase  --ignore-read-groups --reference {REF} {input.vcf} {input.bams} -o {output.vcf}
                bcftools index -f -t {output.vcf}
                whatshap stats {output.vcf} > {output.wstats}
                mkdir -p `dirname {output.bcftools_stats}`
                if [ "{wildcards.region}" = "F" ]; then
                    rm -f {output.bcftools_summary}
                    echo "[DEBUG] bcftools stats AUTO -> {output.bcftools_stats}" >&2
                    bcftools stats -R {INTERSECT_CAPTURE_KIT_AUTO_BED} -F {REF} {output.vcf} > {output.bcftools_stats}
                    python {params.stats_parser} {output.bcftools_stats} {output.bcftools_summary} --sample {wildcards.sample} --region A --append
                    echo "[DEBUG] bcftools stats X -> {output.bcftools_stats}" >&2
                    bcftools stats -R {INTERSECT_CAPTURE_KIT_X_BED} -F {REF} {output.vcf} > {output.bcftools_stats}
                    python {params.stats_parser} {output.bcftools_stats} {output.bcftools_summary} --sample {wildcards.sample} --region X --append
                    echo "[DEBUG] bcftools stats Y -> {output.bcftools_stats}" >&2
                    bcftools stats -R {INTERSECT_CAPTURE_KIT_Y_BED} -F {REF} {output.vcf} > {output.bcftools_stats}
                    python {params.stats_parser} {output.bcftools_stats} {output.bcftools_summary} --sample {wildcards.sample} --region Y --append
                else
                    echo "[DEBUG] bcftools stats region {wildcards.region} -> {output.bcftools_stats}" >&2
                    bcftools stats -R {params.capture_intersect_bed} -F {REF} {output.vcf} > {output.bcftools_stats}
                    python {params.stats_parser} {output.bcftools_stats} {output.bcftools_summary}
                fi
                echo "[DEBUG] merging phased VCF into gVCF" >&2
                python {params.merge_script} {input.gvcf} {output.vcf} {output.tmp_gvcf} {output.mwstats}
                bcftools view {output.tmp_gvcf} -o {output.gvcf}
                bcftools index --tbi {output.gvcf}
            else
                touch {output.vcf}
                touch {output.vcf_tbi}
                touch {output.wstats}
                touch {output.mwstats}
                touch {output.tmp_gvcf}
                touch {output.gvcf}
                touch {output.gvcf_tbi}
                touch {output.bcftools_stats}
                touch {output.bcftools_summary}
            fi
        else
            cp {input.vcf} {output.vcf}
            cp {input.vcf_tbi} {output.vcf_tbi}
            cp {input.gvcf} {output.gvcf}
            cp {input.gvcf_tbi} {output.gvcf_tbi}

            touch {output.tmp_gvcf}
            touch {output.wstats}
            touch {output.mwstats}
            mkdir -p `dirname {output.bcftools_stats}`
            if [ "{wildcards.region}" = "F" ]; then
                rm -f {output.bcftools_summary}
                echo "[DEBUG] (copy branch) bcftools stats AUTO -> {output.bcftools_stats}" >&2
                bcftools stats -R {INTERSECT_CAPTURE_KIT_AUTO_BED} -F {REF} {output.vcf} > {output.bcftools_stats}
                python {params.stats_parser} {output.bcftools_stats} {output.bcftools_summary} --sample {wildcards.sample} --region A --append
                echo "[DEBUG] (copy branch) bcftools stats X -> {output.bcftools_stats}" >&2
                bcftools stats -R {INTERSECT_CAPTURE_KIT_X_BED} -F {REF} {output.vcf} > {output.bcftools_stats}
                python {params.stats_parser} {output.bcftools_stats} {output.bcftools_summary} --sample {wildcards.sample} --region X --append
                echo "[DEBUG] (copy branch) bcftools stats Y -> {output.bcftools_stats}" >&2
                bcftools stats -R {INTERSECT_CAPTURE_KIT_Y_BED} -F {REF} {output.vcf} > {output.bcftools_stats}
                python {params.stats_parser} {output.bcftools_stats} {output.bcftools_summary} --sample {wildcards.sample} --region Y --append
            else
                echo "[DEBUG] (copy branch) bcftools stats region {wildcards.region} -> {output.bcftools_stats}" >&2
                bcftools stats -R {params.capture_intersect_bed} -F {REF} {output.vcf} > {output.bcftools_stats}
                python {params.stats_parser} {output.bcftools_stats} {output.bcftools_summary}
            fi
        fi

        if [ {params.skipsex} -eq 0 ]
        then
            bcftools view -R {params.interval_bed} {output.gvcf} -O z -o {output.gvcf_exome}
            bcftools index -f -t {output.gvcf_exome}
        else
            printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" | bgzip -c > {output.gvcf_exome}
            tabix -f -p vcf {output.gvcf_exome}
        fi
        """
