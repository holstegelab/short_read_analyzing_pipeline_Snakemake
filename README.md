# Main aim of the pipeline

The main aim of the pipeline is to provide a simple and easy-to-use tool for analyzing WGS and WES data.
The pipeline is designed to be used on the server with **SLURM** workload manager and **ZSLURM** add-on to it.
The pipeline is written in **Snakemake** and uses **Conda** for managing dependencies for portable and reproducible analysis.

The default production workflow that we use for cohort calling is:

**DRAGMAP** -> **DeepVariant** -> **GLnexus**

> DRAGMAP is run through `dragen-os`.
> DeepVariant is used for per-sample gVCF calling.
> GLnexus is used for joint cohort calling.


Analysis includes several steps:
* Removing adapters with **AdapterRemoval**
  > AdapterRemoval doesn't delete sequences from reads, so we always have a back-up copy
* Aligning reads with **DRAGMAP** (`dragen-os`)
* Merge different BAM files after DRAGMAP for a single sample *only in case if one sample has more than 1 pair of fastq files*
* Mark duplicates with **Samtools**

  > at this step we create additional **CRAM** files for storage. These files will be saved on the tape and used as back-up copy

* Check mapping quality and possible chimeric reads with custom scripts

  > **bam_stats.py** and **bam_stats_compare_hts.py** calculate alignment statistics
  > If the primary soft-clipped bp ratio is more than 0.5%, **bam_dechimer** is used during BAM merging

* Estimate sample contamination with **verifybamid2**
* Run **DeepVariant** for producing per-sample gVCF files
    > DeepVariant is the default caller used in the production workflow

    > **GATK HaplotypeCaller** is still available as an alternative caller with `--config caller=HaplotypeCaller`
* Combine gVCFs and perform joint calling with **GLnexus** for *cohort calling*
    > GLnexus is the default joint caller used after DeepVariant

    > Alternative version with **GenomicDBImport** or **CombineGVCF** for GATK-based pipeline followed by **GenotypeGVCFs**
* normalization of SNPs with patched version of **bcftools norm**
* get several statistics and produce combined statistic files in dir *stats/* and in tsv-files with extension *.bam_quality.tab*

## Default workflow options

For the default cohort workflow use:

```bash
--config END_POINT=Genotype caller=Deepvariant Combine_gVCF_method=GLnexus
```

Important defaults in the main `Snakefile`:
* `caller=Deepvariant`
* `Combine_gVCF_method=GLnexus`
* `glnexus_filtration=custom`

> `END_POINT=gVCF` stops after per-sample gVCF creation.
> Use `END_POINT=Genotype` or `END_POINT=Combine` when the GLnexus cohort-calling step should be included.

## Additional features

 * SV detection with **delly** (not updated)
 * CNV detection with **cnvkit** (not updated)
 * CNV detection with **GATK** (tested, not updated)
 * chrM analysis
 * Somatic calls for tumor analysis (in test)

## DeepVariant without Apptainer

The production `deepvariant` rule runs DeepVariant 1.9.0 from an extracted
official image, using the dynamic loader and libraries from that image. It
writes the existing `deepvariant/` outputs consumed by downstream phasing and
does not create a container or user namespace at run time.

The shared Snellius runtime is installed at
`/gpfs/work3/0/qtholstg/hg38_res_v2/software/deepvariant-1.9.0-native`. To
prepare it again at a different location (about 6 GB):

```bash
python scripts/prepare_deepvariant_native.py \
  --prefix /path/to/deepvariant-1.9.0-native
```

The default preparation route uses `skopeo`, not Apptainer. An existing Docker
archive can be supplied with `--docker-archive`; on Spider, an existing SIF can
be supplied with `--sif`. Do not move the prepared prefix because its launchers
record that absolute path. Override the configured runtime with
`deepvariant_native_prefix` or `DEEPVARIANT_NATIVE_PREFIX`.

The former Apptainer implementation remains available as an opt-in comparison
target:

```bash
snakemake --snakefile Snakefile DeepVariant_apptainer_all --use-singularity
```

Fallback outputs are isolated under `deepvariant_apptainer/` and are not
connected to the downstream production DAG. Change that comparison directory
with `deepvariant_apptainer_output`.

# HOW TO USE
1. clone this repo on server
2. *If you want use Zslurm*
    1. install Zslurm according to the manual page
    2. Open Zslurm
    3. run pipeline with **snakemake** command

      > snakemake --profile ~/.config/snakemake/zslurm/ --snakefile ~/short_read_analyzing_pipeline_Snakemake/Snakefile --use-conda --use-singularity --rerun-incomplete --retries 0 --config END_POINT=Genotype caller=Deepvariant Combine_gVCF_method=GLnexus

    > **NOTE ABOUT PROFILE**
    > copy zslurm.yaml to ~/.config/snakemake/zslurm/config.yaml and change conda prefix to your conda prefix

   The per-read-group alignment/merge/dechimer/sort fusion is opt-in:

   ```text
   --config fuse_alignment_phases=true alignment_lease_mode=required
   ```

   Keep `fuse_alignment_phases=false` (the default) until the active ZSlurm
   manager and chiefs expose dynamic lease support. `required` performs a
   lease preflight before DRAGMAP; `optional` retains the maximum reservation
   when leases are unavailable; `disabled` is intended only for local tests.

   The sample-level KMC/sex fusion is independently opt-in:

   ```text
   --config fuse_kmer_sex=true kmer_sex_lease_mode=required
   ```

   It keeps the temporary KMC database on assigned node SSD and shrinks from
   2 cores/36 GB to 0.5 core/3 GB before the sex-statistics phase.

   Regional DeepVariant plus Whatshap/merge is independently opt-in:

   ```text
   --config fuse_deepvariant_phasing=true deepvariant_lease_mode=required
   ```

   Raw regional VCF/gVCF files remain on assigned node SSD; the job shrinks
   from 8 cores/10 GB to 1 core/9 GB before phasing and publication.

   BAM/CRAM read-group extraction plus adapter removal is independently
   opt-in; native FASTQ samples remain on the legacy adapter rule:

   ```text
   --config fuse_external_adapter=true external_adapter_lease_mode=required
   ```

   The fused job publishes the legacy temporary raw FASTQ outputs because the
   BAM-tag merge still consumes them, but extraction and adapter processing
   happen only once and adapter processing reads the SSD-local copy.

   chrM/NUMT extraction plus the four realignments is independently opt-in:

   ```text
   --config fuse_chrm_extract_align=true
   ```

   The four intermediate paired FASTQs and alignment scratch files stay on
   the assigned node SSD. The eight legacy BAM/index paths remain unchanged.

   The downstream chrM/NUMT Mutect, merge/filter, and BP-resolution tail is a
   separate opt-in fusion:

   ```text
   --config fuse_chrm_mutect_tail=true
   ```

   It stages the four realigned BAMs once and publishes only both final
   annotated gVCFs and their indexes; all VCF intermediates stay on node SSD.

   The BAM-reading QC fan-out is independently opt-in:

   ```text
   --config fuse_bam_qc=true bam_qc_lease_mode=required
   ```

   VerifyBamID, HS metrics, artifact/OxoG metrics, samtools stats, sampled BAM
   stats, and mosdepth then share one staged markdup BAM and run in parallel on
   node SSD while retaining every legacy QC output path. Each of the six task
   groups atomically releases its relative CPU/memory share when it finishes,
   including on task failure. Use `disabled` only for local tests; `optional`
   retains unreleased capacity when the active chief lacks lease support.

4. If you want to run just several steps (for example only Alignment step) -
choose suitable **smk** file as **--snakefile** or use `END_POINT`

Available `END_POINT` values:
* `Align` - run only alignment/preprocessing workflow
* `gVCF` - run alignment and per-sample gVCF calling
* `Genotype` - run alignment, gVCF calling and GLnexus joint calling
* `Combine` - run GLnexus joint calling endpoint

We use `{cohort}.tsv` file as the start file for a pipeline and `{cohort}.source` file for additional paths to the files.
**These 2 files should be uploaded by the user to the server before starting the pipeline.**

### `{cohort}.tsv` - main file with filenames and paths (8 or 9 columns).

#### Column Descriptions

| Column Name             | Description                                                                                                                                                                                                                                                                                  |
|-------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **study**               | Study name, e.g., FR_Lille, DE_Bonn. It should contain a two-letter country code + a short study name (without underscores).                                                                                                                                                                 |
| **sample_id**           | Identifier of the sample in the form `studyname_sampleid` (to avoid duplicates between studies). Example: `FR_Lille_B00E9CH`. Subjects are linked to samples in the phenotype file. A subject can have multiple samples (e.g., WGS and WES), so `sample_id` is not the same as `subject_id`. |
| **file_type**           | Can be one of: `fastq_paired`, `bam`, `cram`, `recalibrated_bam`, `recalibrated_cram`, `sra`. Recalibrated files use GATK BQSR. The pipeline attempts to recover original quality data (OQ tags). Preferably, submit non-recalibrated data with BAM/CRAM files including unmapped reads.     |
| **sample_type**         | Can be one of: `illumina_exome`, `illumina_wgs`.                                                                                                                                                                                                                                             |
| **capture_kit**         | (For exomes only) Example: `Agilent_V5`. The corresponding capture kit file (`Agilent_V5.bed`) should be uploaded to the `resources/capture_kits` directory. For build 37 capture kits, upload as `{capture_kit}.b37.bed`. These files will be lifted over to build 38.                      |
| **sex**                 | Sex of the sample for validation and chrX/Y calling: `F` or `M`.                                                                                                                                                                                                                             |
| **filenames_read1**     | (For `fastq_paired` files) The file containing read1 sequences. If using BAM, CRAM, or interleaved FASTQ files, only this column is needed.                                                                                                                                                   |
| **filenames_read2**     | (For `fastq_paired` files) The file containing read2 sequences. For CRAM files, this column should contain the reference FASTA file needed for decoding.                                                                                                                                      |
| **additional_commands** | (Optional) Additional commands to run samples with.                                                                                                                                                                                                                                          |

Available additional commands:
`no_dedup=True` - skip deduplication (markdup) step

### `{cohort}.source` - additional file with paths to the files.
If your files are stored on dCache or archive, you can use this file to specify the paths to the files.
The file contains one line.

- NFS archive: `archive://archive/hulsmanm/source_files/UCL_NIH/`
- dCache: `dcache:<remote>:/<root>`, for example `dcache:mine_hiseq2000:/`

For dCache, put `<remote>.conf` next to the cohort listing. The remote name
must match the rclone section in that config. Absolute-looking paths in the
TSV are interpreted below the selected dCache remote root. The pipeline
stages stable batches directly from Snellius. Each sample is copied to active
storage with `dcache_cp`, Adler-32 verified, and then its stage pin is released.

### `{cohort}.target` - processed-data destination.

This optional one-line file controls where processed pipeline output is
written. A dCache target uses the same URI form:

`dcache:<remote>:/<processed-root>`

Put the corresponding `<remote>.conf` next to the listing. CRAMs, statistics,
Kraken results and region-level gVCF bundles retain their existing
subdirectory layout below this root. Uploads run directly from Snellius and are only
marked copied after the remote Adler-32 checksum matches.

Bulk staging uses the managed `ada` default pin lifetime of 7 days. Relevant
tuning keys are `dcache_stage_poll_seconds`, `dcache_stage_timeout`, and
`dcache_download_workers`.

Transfer jobs request `dcache_download_slots=1` or `dcache_upload_slots=1` from
zslurm. The manager-wide maxima are configured in zslurm, so downloads and
uploads are throttled independently. `dcache_download_lock_slots` (default 4)
is a separate, download-only advisory lock around `dcache_cp`; the old pipeline
key `dcache_transfer_slots` remains a fallback for that local lock.



#### File Handling Notes
- For `fastq_paired` files, paired read data is split across two files: `filenames_read1` for read1 and `filenames_read2` for read2.
- For BAM/CRAM or interleaved FASTQ files, only `filenames_read1` is used.
- If a CRAM file is used, the reference file must be specified in `filenames_read2` and delivered separately.
- For multiple read groups, files are comma-separated (`,`). Order should be maintained (e.g., lane 1 first, then lane 2, etc.).
- BAM/CRAM files containing multiple read groups do not need extra annotations. Read group information is parsed directly from the file metadata.

This documentation ensures that all file naming conventions and metadata structures are followed correctly for pipeline processing.
