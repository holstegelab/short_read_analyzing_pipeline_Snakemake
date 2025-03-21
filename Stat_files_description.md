# In this file, we will describe the Stat files



We use `{cohort}.tsv` file as the start file for a pipeline and `{cohort}.source` file for additional paths to the files.
**These 2 files should be uploaded by the user to the server before starting the pipeline.**

### `{cohort}.tsv` - main file with filenames and paths (7 or 8 columns).

#### Column Descriptions

| Column Name          | Description |
|----------------------|-------------|
| **study**           | Study name, e.g., FR_Lille, DE_Bonn. It should contain a two-letter country code + a short study name (without underscores). |
| **sample_id**       | Identifier of the sample in the form `studyname_sampleid` (to avoid duplicates between studies). Example: `FR_Lille_B00E9CH`. Subjects are linked to samples in the phenotype file. A subject can have multiple samples (e.g., WGS and WES), so `sample_id` is not the same as `subject_id`. |
| **file_type**       | Can be one of: `fastq_paired`, `bam`, `cram`, `recalibrated_bam`, `recalibrated_cram`, `sra`. Recalibrated files use GATK BQSR. The pipeline attempts to recover original quality data (OQ tags). Preferably, submit non-recalibrated data with BAM/CRAM files including unmapped reads. |
| **sample_type**     | Can be one of: `illumina_exome`, `illumina_wgs`. |
| **capture_kit**     | (For exomes only) Example: `Agilent_V5`. The corresponding capture kit file (`Agilent_V5.bed`) should be uploaded to the `resources/capture_kits` directory. For build 37 capture kits, upload as `{capture_kit}.b37.bed`. These files will be lifted over to build 38. |
| **sex**             | Sex of the sample for validation and chrX/Y calling: `F` or `M`. |
| **filenames_read1** | (For `fastq_paired` files) The file containing read1 sequences. If using BAM, CRAM, or interleaved FASTQ files, only this column is needed. |
| **filenames_read2** | (For `fastq_paired` files) The file containing read2 sequences. For CRAM files, this column should contain the reference FASTA file needed for decoding. |

#### File Handling Notes
- For `fastq_paired` files, paired read data is split across two files: `filenames_read1` for read1 and `filenames_read2` for read2.
- For BAM/CRAM or interleaved FASTQ files, only `filenames_read1` is used.
- If a CRAM file is used, the reference file must be specified in `filenames_read2` and delivered separately.
- For multiple read groups, files are comma-separated (`,`). Order should be maintained (e.g., lane 1 first, then lane 2, etc.).
- BAM/CRAM files containing multiple read groups do not need extra annotations. Read group information is parsed directly from the file metadata.

This documentation ensures that all file naming conventions and metadata structures are followed correctly for pipeline processing.

## Combined statistics files

### `{cohort}.bam_quality.tab` - BAM Quality Statistics
This file contains various quality metrics for BAM files used in the study. The columns are grouped as follows:

#### General Sample Information
- **sample**: Sample identifier.
- **total_sequences**: Total number of sequences in the sample.
- **avg_quality**: Average base quality score.
- **average_length, max_length**: Average and maximum read lengths.

#### Read Alignment Metrics
- **error_rate**: Estimated error rate in sequencing.
- **insert_size_avg, insert_size_std**: Mean and standard deviation of insert sizes.
- **reads_unmapped, reads_properly_paired, reads_mq0, reads_duplicated**: Read mapping and pairing statistics.
- **pairs_inward, pairs_outward, pairs_other, pairs_diff_chromosomes**: Pair orientation information.

#### Base and Coverage Statistics
- **bases_total, bases_mapped**: Total and mapped base counts.
- **qual_by_cycle, base_fractions**: Quality and base composition over cycles.
- **gc35_depth, gc50_depth**: GC content-based depth metrics.

#### Exome-Specific Metrics
- **exome_total_sequences, exome_reads_mq0**: Exome-specific sequence statistics.
- **exome_insert_size_avg, exome_insert_size_std**: Insert size metrics for exomes.
- **exome_pairs_inward, exome_pairs_outward, exome_pairs_other, exome_pairs_diff_chromosomes**: Exome-specific pairing details.
- **exome_insertions_frac, exome_deletions_frac**: Insertion and deletion fractions in exomes.
- **exome_coverage1_frac, exome_coverage1000plus_frac**: Coverage metrics at different depths.
- **exome_gc35_depth, exome_gc50_depth**: GC-based coverage statistics.

#### Quality Control Metrics
- **pca2_freemix**: Contamination estimate based on PCA.
- **unmapped_ratio_all, mqual20_ratio_all, duplicate_ratio_all**: Various mapping quality ratios.
- **soft_clipped_bp_ratio_all, aligned_bp_ratio_all**: Soft clipping and alignment percentages.
- **poly_a_all, poly_g_all, illumina_adapter_all, pcr_adapter_1_all, pcr_adapter_2_all, nextera_all**: Adapter content statistics.

#### Pre-Adapter and Bait Bias Metrics
- **pre_adapter_ac, pre_adapter_ag, ..., pre_adapter_tg**: Pre-adapter bias measurements.
- **bait_bias_ac, bait_bias_ag, ..., bait_bias_tg**: Bias introduced by bait hybridization.
- **bait_design_efficiency**: Efficiency of bait capture.

#### Targeted Sequencing Metrics
- **on_bait_bases, near_bait_bases, off_bait_bases**: Base distribution related to baits.
- **mean_bait_coverage, pct_usable_bases_on_bait, pct_usable_bases_on_target**: Bait-target efficiency metrics.
- **pct_target_bases_1x, pct_target_bases_10x, ..., pct_target_bases_500x**: Target base coverage at different depths.

#### Dropout Metrics
- **at_dropout, gc_dropout**: AT/GC dropout percentages.
- **het_snp_sensitivity**: Heterozygous SNP sensitivity measure.

### `{cohort}`.bam_rg_quality.tab - Read Group Quality Statistics
This file contains quality metrics for individual read groups in the study. The columns are grouped as follows:

#### General Information
- **sample**: The sample identifier.
- **readgroup**: The read group identifier.

#### Alignment Ratios and Read Retention
- **ar_read_pairs**: Number of read pairs analyzed.
- **ar_aligned_fraction**: Fraction of reads aligned to the reference genome.
- **ar_adapter_fraction**: Fraction of reads containing adapter sequences.
- **ar_retained_fraction**: Fraction of reads retained after filtering.
- **ar_average_read_length_retained**: Average read length after retention.

#### Adapter and Alignment Counts
- **ai_reads_aligned**: Number of reads successfully aligned.
- **ai_n_adapter**: Total number of detected adapters.
- **ai_adapter1(/2)**: Number of reads with adapter 1/2.

#### Merge Statistics
- **merge_fragments**: Number of merged fragments.
- **merge_alignments**: Number of merged alignments.
- **merge_supplementary_alignments_ratio**: Ratio of supplementary alignments.
- **merge_total_bp**: Total base pairs in merged alignments.
- **merge_primary_soft_clipped_bp_ratio**: Ratio of soft-clipped bases in primary alignments.
- **merge_supplementary_bp_ratio**: Ratio of base pairs in supplementary alignments.
- **merge_readded_fragments_ratio**: Ratio of re-added fragments.
- **merge_restored_read_ratio**: Ratio of restored reads.
- **merge_restored_bp_ratio**: Ratio of restored base pairs.

#### Mapping Statistics
- **dm_total_input_reads**: Total input reads.
- **dm_mapped_reads**: Number of reads mapped to the reference.
- **dm_properly_paired_reads**: Number of properly paired reads.
- **dm_total_bases**: Total number of bases in reads.
- **dm_mapped_bases**: Total number of bases mapped to the reference.
- **dm_unmapped_reads**: Number of unmapped reads.
- **dm_singleton_reads_itself_mapped_mate_unmapped**: Singleton reads where the mate is unmapped.
- **dm_paired_reads_itself_mate_mapped**: Number of paired reads where the mate is mapped.
- **dm_not_properly_paired_reads_discordant**: Number of discordant paired reads.
- **dm_paired_reads_mapped_to_different_chromosomes_mapq_10**: Number of paired reads mapped to different chromosomes with MAPQ ≥ 10.

#### Mapping Quality
- **dm_reads_with_mapq_40_inf**, **dm_reads_with_mapq_30_40**, **dm_reads_with_mapq_20_30**, **dm_reads_with_mapq_10_20**, **dm_reads_with_mapq_0_10**, **dm_reads_with_mapq_na_unmapped_reads**: Reads categorized by mapping quality ranges (≥40, 30-40, 20-30, 10-20, 0-10, and unmapped reads with no MAPQ).

#### Read Modification Statistics
- **dm_reads_with_indel**: Reads containing insertions or deletions.
- **dm_soft_clipped_bases**: Number of soft-clipped bases.
- **dm_mismatched_bases_excl_indels**: Number of mismatched bases, excluding indels.
- **dm_q30_bases**: Number of bases with quality score ≥ 30.
- **dm_soft_clipped_bases_diff**: Difference in soft-clipped bases before and after processing.
- **dm_reads_with_indel_diff**: Difference in number of reads with indels before and after processing.
- **dm_mismatched_bases_diff**: Difference in mismatched bases before and after processing.
- **dm_q30_bases_diff**: Difference in number of Q30 bases before and after processing.

#### Dechimerization Statistics
- **dechimer_fragment_modified_ratio**: Ratio of modified fragments after dechimerization.
- **dechimer_clip_ratio**: Ratio of clipped reads due to dechimerization.
- **dechimer_loose_end_clip_ratio**: Ratio of loose-end clipped reads.
- **dechimer_all_sup_discarded_ratio**: Ratio of all supplementary alignments discarded.
- **dechimer_min_alignment_length_unmap_ratio**: Ratio of alignments unmapped due to minimum alignment length filtering.

### Sex Chromosome Analysis (`sex_chrom.tab`)
This file contains statistics related to sex chromosome determination.

### General Information
- **sample**: The sample identifier.
- **reported**: Reported sex of the sample.
- **sex**: Predicted sex based on chromosomal data.

#### Chromosome Statistics
- **yyratio**: Ratio of Y chromosome reads to total reads.
- **auto**: Autosomal chromosome coverage.
- **chrx**: X chromosome coverage.
- **chry**: Y chromosome coverage.

#### Regression Analysis
- **intercept**: Intercept of the regression model for chromosome coverage.
- **slope**: Slope of the regression model.
- **rvalue**: Correlation coefficient indicating the fit of the model.

#### Distance Metrics
- **xratio**, **yratio**: Ratios of X and Y chromosome coverage to autosomal coverage.
- **dist_auto**: Distance metric for autosomal coverage.
- **dist_chrx**, **dist_chrx_ratio**: Distance metric and ratio for X chromosome coverage.
- **dist_chry**, **dist_chry_ratio**: Distance metric and ratio for Y chromosome coverage.

### Coverage metrics (`{cohort}.coverage.hdf5`)
This file contains coverage metrics for the study samples. The data is stored in HDF5 format for efficient storage and retrieval.

