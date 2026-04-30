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

# HOW TO USE
1. clone this repo on server
2. *If you want use Zslurm*
    1. install Zslurm according to the manual page
    2. Open Zslurm
    3. run pipeline with **snakemake** command

      > snakemake --profile ~/.config/snakemake/zslurm/ --snakefile ~/short_read_analyzing_pipeline_Snakemake/Snakefile --use-conda --use-singularity --rerun-incomplete --retries 0 --config END_POINT=Genotype caller=Deepvariant Combine_gVCF_method=GLnexus

    > **NOTE ABOUT PROFILE**
   > copy zslurm.yaml to ~/.config/snakemake/zslurm/config.yaml and change conda prefix to your conda prefix

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
File is 1-liner with the following structure:
`{dcache or archive}://{path to the files}` for example:  `archive://archive/hulsmanm/source_files/UCL_NIH/`
`/archive/hulsmanm/source_files/UCL_NIH/` will be added to the path from the `{cohort}.tsv` file.



#### File Handling Notes
- For `fastq_paired` files, paired read data is split across two files: `filenames_read1` for read1 and `filenames_read2` for read2.
- For BAM/CRAM or interleaved FASTQ files, only `filenames_read1` is used.
- If a CRAM file is used, the reference file must be specified in `filenames_read2` and delivered separately.
- For multiple read groups, files are comma-separated (`,`). Order should be maintained (e.g., lane 1 first, then lane 2, etc.).
- BAM/CRAM files containing multiple read groups do not need extra annotations. Read group information is parsed directly from the file metadata.

This documentation ensures that all file naming conventions and metadata structures are followed correctly for pipeline processing.
