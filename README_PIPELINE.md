# Pipeline Overview

## Scope

This document summarizes the short read Snakemake pipeline, focusing on per-sample processing, statistics packaging, current archival behaviour, and the path toward region-level tarball uploads of gVCFs and exome subsets.

## Per-sample flow (WGS and WES)

1. **Staging and space reservation** (`Aligner.smk`)
   - Route-specific `start_sample_*` rules reserve active storage and validate or materialize inputs from active disk, `/archive`, dCache, or requester-pays S3 (`get_source_files`).
   - `get_readgroups` checkpoint populates `SAMPLEINFO` read group metadata when missing.
   - A `dcache:<remote>:/path` `.source` is staged as a stable batch directly from Snellius and downloaded with Adler-32 verification. An `s3://bucket/path` `.source` is downloaded per sample on a compute node with the configured AWS profile and `--request-payer requester`; AWS credentials never belong in the sample sheet. S3 downloads request ZSlurm's independent `s3_download_slots` pool. During a rolling manager upgrade, the executor falls back to the dCache-download pool so concurrency remains bounded.

   For an S3-backed sheet, put only the common object prefix in the sidecar,
   for example `s3://wanglab-dss-share/distribution/adsp/cram`. Keep column 7
   relative to that prefix and set the ninth-column `filesize=<GiB>` value from
   catalog metadata so active-storage accounting is resolved without querying
   S3 while building the DAG. The route is non-interactive and uses the normal
   AWS credential chain; anonymous `--no-sign-request` access is not assumed.
   Optional workflow configuration keys are `aws_cli`, `s3_max_attempts`,
   `s3_initial_backoff_seconds`, and `s3_max_backoff_seconds`.

2. **Optional readgroup split & FASTQ preparation**
   - `split_alignments_by_readgroup` splits multi-RG BAM/CRAMs.
   - `external_alignments_to_fastq` or `{sample}.{readgroup}.fastq.cut_*.fq.gz` produced via `adapter_removal` + `adapter_removal_identify` (`Aligner.smk`).

3. **Alignment and merge**
   - `align_reads_fused` runs DRAGMAP and its per-readgroup merge/check/dechimer/sort tail.
   - `markdup` consumes all validated readgroup BAMs, creates a merged BAM only on assigned node SSD when needed, and performs duplicate marking there. Only the final markdup BAM/index/statistic are atomically published to shared storage.
   - For external inputs, the staged source is reclaimed as soon as all readgroup BAMs, indexes, and validation markers exist; markdup does not pin or read the source. The long-lived active-storage reservation is 85% of the estimated peak. Markdup temporarily acquires the remaining 15% while its final BAM can overlap the readgroup BAMs, then releases that transient share together with completed intermediates.

4. **QC, contamination, duplicates**
   - `merge_rgs_badmap` for contamination FASTQs.
   - `markdup` supplies the final analysis BAM to the downstream QC and caller rules.
   - `get_validated_sex`, `verifybamid`, `hs_stats`, `Artifact_stats`, `samtools_stat*`, `bamstats_*`, `coverage`, `gatherstats` family generate QC artefacts in `stat/` and samplefile-level tabs.

5. **gVCF generation**
   - `HaplotypeCaller` produces region-split gVCFs; `reblock_gvcf` compresses blocks. For DeepVariant, equivalent pipeline in `Deepvariant.smk` culminating in `DVWhatshapPhasingMerge`.
   - `extract_exomes_gvcf` (Snakemake rule in `gVCF.smk`) produces WES subset for WGS samples via `SelectVariants`; WES samples bypass extraction (`cp`).
   - `gvcf_sample_done` marks completion once all region files exist (level 1 for WGS, level 0 for WES).

6. **Finishing and archival**
   - `stat_sample_done` creates `{sample}.done` sentinel used by `Snakefile::finished_sample` and `tar_stats_per_sample`.
   - `tar_stats_per_sample` packages per-sample QC into `stat/{sample}.stats.tar.gz`, currently using `tar --remove-files` to delete inputs.
   - `cram_encrypt_fused` creates CRAM/CRAI on assigned node SSD, encrypts and uploads both directly from SSD, verifies both Adler-32 checksums and publishes only `{sample}.mapped_hg38.cram.ADLER32` plus `{sample}.mapped_hg38.cram.copied` on active storage. Its weighted `0.05` upload-slot reservation reflects the measured 3.8% upload duty without creating a phase-time wait.
   - If the cohort has a `dcache:<remote>:/path` `.target`, all processed-data upload helpers use that remote/config directly from Snellius and retain the established `cram/`, `stat/`, `kraken/`, `gvcf/`, and `chrM/` layout below the target root.
   - `Snakefile::finished_sample` touches `{SOURCEDIR}/{sample}.finished` after CRAM copy, gVCF done, stats done, Kraken output ready.
- `{sample}.done` (from `Stat.smk::stat_sample_done`) only asserts statistics availability; `{sample}.finished` is the global completion marker consumed by `Snakefile::finished_sample` to release active-storage reservations.

## Statistics packaging dependencies

- **Per-sample QC artefacts only** `Stat.smk::tar_stats_per_sample` now takes explicit per-sample inputs (`{sample}.hs_metrics`, `{sample}.samtools.stat`, `{sample}.bam_all.tsv`, `{sample}.markdup.stat`, etc.) plus readgroup logs enumerated via `get_rg_files`. Cohort-wide tables (`{samplefile}.*.tab`/`.hdf5`) are no longer listed as inputs.
- **Temp outputs at source** The producing rules (`hs_stats`, `samtools_stat`, `Artifact_stats`, `bamstats_all`, `verifybamid`, `coverage`) now emit their main artefacts as `temp(...)`. Snakemake keeps them until `tar_stats_per_sample` (and other consumers) finish, then cleans them automatically.
- **Sample readiness** `tar_stats_per_sample` still depends on `{sample}.done` via `stat_sample_done`, and `Snakefile::finished_sample` uses the same sentinel chain to ensure CRAM, gVCF, stats, and Kraken outputs complete before marking `{sample}.finished`.

## Artefacts and persistence

Persistent by default:
- `CRAM/{sample}.mapped_hg38.cram.copied` and `.ADLER32` upload receipts.
- The encrypted CRAM and its CRAI are durable at the configured dCache target. Plaintext CRAM, encrypted CRAM and CRAI exist only in the fused job's assigned SSD directory; active storage holds only the checksum and completion receipt.
- `stat/{sample}.stats.tar.gz` (but original inputs removed; consider `temp()` refactor).
- Cohort files: `{samplefile}.bam_quality.tab`, `{samplefile}.bam_rg_quality.tab`, `{samplefile}.oxo_quality.tab`, `{samplefile}.coverage.hdf5`, `{samplefile}.sex_chrom.tab`.
- gVCF outputs: `gVCF/reblock/{region}/{sample}.{region}.wg.vcf.gz` (WES, or WGS before exome extraction), `gVCF/exome_extract/{region}/{sample}.{region}.wg.vcf.gz` (WGS exome slices), `Deepvariant/gVCF/...` equivalents when DeepVariant active.
- `SOURCEDIR/{sample}.finished` sentinel, plus Kraken deliverables.

Temporary/auto-cleaned (many declared with `temp()` or candidates):
- Intermediate FASTQs (`FQ/`, `FQ_BADMAP/`), raw HaplotypeCaller outputs (`gVCF/raw/...`), DragSTR models, aligner temp outputs, etc.
- `stat/{sample}.done` is `temp(touch(...))` and removed after downstream completion.

## `extract_exomes` behaviour (WGS vs WES)

- Rule `extract_exomes_gvcf` (in `gVCF.smk`) operates per `{sample, region}`.
  - For WGS (`"wgs" in SAMPLEINFO[sample]['sample_type']`), runs GATK `SelectVariants -L {region_to_file(... interval_list)} -ip 500` to clip the reblocked gVCF to exome capture bins, generate `gVCF/exome_extract/{region}/...`.
  - For WES samples, simply copies reblocked gVCF to exome_extract path, because it already contains only exome intervals.
- `Deepvariant.smk::extract_exomes_dv` mirrors behaviour on `DEEPVARIANT/gVCF/...` outputs using `DVWhatshapPhasingMerge` as source; handles sex-based skipping (`skipsex`).

## WGS region-level handling today

- Regions defined in `common.py` via `level0_range`, `level1_range`, etc.; `level1_regions` (~10 autosomal chunks + X/Y haploid splits) used for per-sample gVCF outputs.
- `Combine_gVCF.smk::combinegvcfs` merges per-region sample gVCFs into cohort-level `GVCF/MERGED/{samplefile}.{region}.wg.vcf.gz` using WGS interval (`region_to_file(..., wgs=True)`), dropping to level0 for WES samples.
- `DBImport.smk`/`GLnexus.smk` `generate_gvcf_input_*` functions map WGS vs WES to appropriate region granularity (e.g., `convert_to_level1` for WGS, `convert_to_level0` for WES) when building GenomicsDB or GLnexus inputs.

## Path to region-level tarballs (future work)

Target: after per-region gVCF (and exome subset) generation across cohort, create tar archives grouped by region, suitable for dCache upload.

### Current assets
- Region definitions (`level1_regions`, `level2_regions`, etc.) already used for per-sample gVCFs.
- Per-sample gVCFs live under `GVCF/reblock/{region}/{sample}.{region}.wg.vcf.gz` (and DeepVariant equivalents). These are the files to tar together, grouped by region at level2, with padded interval lists already baked into the upstream calling rules.

### Required steps

1. **Decide archive granularity**
   - Whole-genome gVCFs: tar per-sample level2 region files from `GVCF/reblock/{region}/` so each archive contains all samples for a region (gVCF + `.tbi`) using the padded bins already used upstream.
   - Exome gVCFs: use `GVCF/exome_extract/{region}/` for WGS-derived exomes plus `convert_to_level0` mapping for WES samples; DeepVariant equivalents sit under `Deepvariant/gVCF/...`.

2. **Snakemake rule additions**
   - New rule `tar_gvcf_region` iterating over `level1_regions` (or chosen level) to collect all `{sample}.{region}.wg.vcf.gz` plus indexes, output `GVCF/regions/{region}.gvcf.tar.gz`. Use Snakemake `expand()` to include WES copies (converted via `convert_to_level0` or `convert_to_level1` as needed).
   - Equivalent rule for exome subset e.g., `GVCF/exome_regions/{region}.gvcf.tar.gz` drawing from `gVCF/exome_extract/` (WGS) plus `gVCF/reblock/` for WES (converted region names).
   - Consider per-samplefile vs global dataset grouping; align with `Combine_gVCF.smk` naming for cohort-scope vs sample-scope deliverables.

3. **Dependency wiring**
   - Ensure tar rule inputs include combined stats or `gvcf_sample_done` to guarantee all sample-level gVCFs are ready. Example: `input: expand(pj(GVCF, "exome_extract", region, f"{sample}.{region}.wg.vcf.gz"), sample=sample_names)`.
   - For exome tar, ensure WES sample region mapping uses `convert_to_level0` to find correct file names.

4. **Use `temp()` to manage storage**
   - Keep original gVCFs persistent until remote transfer verified; mark derived tars as final deliverables. Optionally mark per-sample gVCFs as `temp()` once tar + upload succeed to reclaim space.

5. **Upload automation**
   - Mirror the verified-upload boundary in `Encrypt.smk::cram_encrypt_fused`: upload region tarballs and exome tarballs, verify checksums, and only then touch a `.copied` sentinel (e.g., `GVCF/{region}.tar.copied`).
   - Add configuration toggles for enabling gVCF/exome uploads.

6. **Update combined workflows**
   - `Snakefile` `END_RULE`/`CLEAN_RULE` should include new targets so `snakemake --cleanup-metadata` removes adhesives, and `END_POINT` gating ensures full DAG.

## Open questions / next actions
- Confirm desired grouping (per-cohort vs per-sample) for region tarballs to size outputs appropriately.
- Audit which upstream rules need `temp()` adjustments (e.g. `Stat.smk::stat_sample_done` inputs) before removing `--remove-files` in tar steps.
- Extend config schema with toggles for gVCF/exome tar creation and dCache uploads, mirroring existing CRAM encryption settings.

## TODO for implementation
- Define target file list for region tar (per-sample vs per-cohort). The current per-sample structure suggests tar-by-region across samples to support distributed retrieval.
- Implement Snakemake rules for tar creation and optional upload, ensuring they respect sex-specific region omissions (skip Y/H for female). Leverage existing `generate_gvcf_input_*` functions as reference for mapping sample types to region splits.
- Review chunking level: Level1 (10 autosomal + X/Y) may suit tar size; Level2 (100 autosomal partitions) yields smaller files but more tars.
- Evaluate concurrency impact on disk usage; ensure new tar rules declare resources to avoid oversubscription.
- Adjust `Stat.smk::tar_stats_per_sample` to drop `--remove-files` and mark upstream outputs `temp()` for consistent behaviour (optional but recommended for clarity).
