# Optional comparison workflow; outputs never feed the production phased DAG.
# Included by Deepvariant.smk to preserve the DeepVariant_apptainer_all target.

rule deepvariant_apptainer:
    """Opt-in DeepVariant 1.9.0 fallback using Apptainer."""
    input:
        bed = region_to_bed_file,
        bed_wgs = region_to_bed_file_wgs,
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        validated_sex=pj(KMER,"{sample}.result.yaml"),
    output:
        vcf = pj(DEEPVARIANT_APPTAINER,'VCF', "{region}","{sample}.{region}.vcf.gz"),
        vcf_tbi = pj(DEEPVARIANT_APPTAINER,'VCF', "{region}","{sample}.{region}.vcf.gz.tbi"),
        gvcf = pj(DEEPVARIANT_APPTAINER,'gVCF', "{region}","{sample}.{region}.g.vcf.gz"),
        gvcf_tbi = pj(DEEPVARIANT_APPTAINER,'gVCF', "{region}","{sample}.{region}.g.vcf.gz.tbi")
    params:
            mode=get_sequencing_mode,
            haploid_contigs=lambda wildcards: 'chrX,chrX_KI270880v1_alt,chrX_KI270881v1_alt,chrX_KI270913v1_alt,chrY,chrY_KI270740v1_random' if wildcards['region'].endswith("H") else 'chrNONE',
            skipsex = lambda wildcards, input: int(get_validated_sex_file(input) == 'female' and wildcards['region'].startswith('Y')),
            inter_dir = pj(DEEPVARIANT_APPTAINER,'DV_intermediate'),
            # check = CHECKEMPTY
    container: 'docker://google/deepvariant:1.9.0'
    resources:
        n="7",
        nshards=8,
        # Limit concurrent Apptainer launches in Snakemake with
        # --resources deepvariant_container_slots=N. Without this, a retry
        # wave can exhaust the pilot node's user-namespace quota.
        deepvariant_container_slots=1,
        mem_mb=get_mem_mb_deepvariant,
        time = get_time('deepvariant_apptainer'),
        ssd_use="possible",
        # Live WGS region jobs use 0.9--1.4 GiB. Keep room for a larger
        # interval, TFRecords and transient post-processing files.
        ssd_gb=4
    shell:
        """
        if [ {params.skipsex} -eq 0 ]
        then
            TMP_SSD="/scratch-node/${{USER}}.${{SLURM_JOB_ID}}"
            JOB_ID="${{SLURM_JOB_ID}}"
            if [ -z "$JOB_ID" ]; then JOB_ID="${{SLURM_JOBID}}"; fi
            if [ -z "$JOB_ID" ]; then JOB_ID="$$"; fi
            if [ ! -d "$TMP_SSD" ] || [ ! -w "$TMP_SSD" ]; then CAND=$(ls -1dt /scratch-node/${{USER}}.* 2>/dev/null | head -n1 || true); if [ -n "${{CAND:-}}" ] && [ -d "$CAND" ] && [ -w "$CAND" ]; then TMP_SSD="$CAND"; fi; fi
            if [ -d "$TMP_SSD" ] && [ -w "$TMP_SSD" ]; then RUNDIR_BASE="$TMP_SSD/deepvariant_apptainer/$JOB_ID"; elif [ -n "${{SLURM_TMPDIR:-}}" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then RUNDIR_BASE="$SLURM_TMPDIR/deepvariant_apptainer/$JOB_ID"; else RUNDIR_BASE="{params.inter_dir}"; fi
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

def get_deepvariant_apptainer_gvcfs(wildcards):
    files = []
    for sample in sample_names:
        regions = level1_regions if 'wgs' in SAMPLEINFO[sample]['sample_type'] else level0_regions
        for region in regions:
            gvcf = pj(DEEPVARIANT_APPTAINER, 'gVCF', region, f'{sample}.{region}.g.vcf.gz')
            files.extend([gvcf, gvcf + '.tbi'])
    return files


rule DeepVariant_apptainer_all:
    input:
        get_deepvariant_apptainer_gvcfs


# python {params.check} {output.vcf}
# python {params.check} {output.gvcf}

