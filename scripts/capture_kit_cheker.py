#!/usr/bin/env python3
"""
Capture Kit Finder - A script to identify which capture kit was used in targeted sequencing

This script analyzes a BAM file and compares its coverage against multiple BED files
representing different capture kits. It calculates various metrics to determine
which kit was most likely used in the experiment.

Designed to be integrated into a Snakemake pipeline, this script produces a simple TSV file
with coverage comparison metrics for all tested capture kits, ranked by likelihood.

This version supports precomputation of capture kit data for improved performance.

Usage:
    # Precompute mode:
    python capture_kit_finder.py --precompute_mode --kit_dir /path/to/bedfiles/ --genome /path/to/genome.fa --precompute_output precomputed_kits.json

    # Analysis mode:
    python capture_kit_finder.py --bam /path/to/sample.bam --precomputed precomputed_kits.json --output results/sample1 --expected_kit KitName

Output:
    - A TSV file with coverage metrics for all kits

Requirements:
    - pysam
    - polars
    - numpy
    - pybedtools
"""

import os
import argparse
import pysam
import polars as pl
import numpy as np
from pybedtools import BedTool
import time
import logging
import json
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('CaptureKitFinder')


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Identify which capture kit was used in a targeted sequencing experiment')

    # Common parameters
    parser.add_argument('--threads', type=int, default=multiprocessing.cpu_count(), help='Number of threads to use')
    parser.add_argument('--temp_dir', default='/tmp', help='Directory for temporary files')
    parser.add_argument('--min_mapping_quality', type=int, default=20, help='Minimum mapping quality for reads')

    # Mode selection
    parser.add_argument('--precompute_mode', action='store_true',
                        help='Run in precomputation mode to generate data for capture kits')

    # Precomputation mode parameters
    parser.add_argument('--kit_dir', help='Directory containing BED files for different capture kits')
    parser.add_argument('--genome', help='Reference genome file (genome.txt format with chromosome sizes)')
    parser.add_argument('--precompute_output', help='Output file for precomputed capture kit data')

    # Analysis mode parameters
    parser.add_argument('--bam', help='Input BAM file')
    parser.add_argument('--precomputed', help='Precomputed capture kit data file')
    parser.add_argument('--output', default='capture_kit_results',
                        help='Output prefix for results files (without extension)')
    parser.add_argument('--min_coverage', type=float, default=5.0,
                        help='Minimum mean coverage required to consider a kit (default: 5.0)')
    parser.add_argument('--min_percent_covered', type=float, default=30.0,
                        help='Minimum percentage of kit bases that must be covered to consider a kit (default: 30.0%)')
    parser.add_argument('--expected_kit', help='Direct input of the expected capture kit name')

    return parser.parse_args()


def check_bam_index(bam_file):
    """Check if BAM index exists, create if it doesn't"""
    index_file = bam_file + '.bai'
    if not os.path.exists(index_file):
        logger.info(f"BAM index not found, creating index for {bam_file}")
        pysam.index(bam_file)
    return True


def process_bed_file(args):
    """Process a single BED file to extract information for precomputation

    Args:
        args: Tuple containing (bed_file, genome_file, output_dir)

    Returns:
        Tuple of (kit_name, kit_data) or (None, None) if an error occurs
    """
    bed_file, genome_file, output_dir = args

    try:
        kit_name = os.path.basename(bed_file).replace('.bed', '')
        logger.info(f"Processing kit: {kit_name}")

        # Check if file exists and has content
        if not os.path.exists(bed_file):
            logger.error(f"BED file does not exist: {bed_file}")
            return None, None

        file_size = os.path.getsize(bed_file)
        if file_size == 0:
            logger.error(f"BED file is empty: {bed_file}")
            return None, None

        # Get absolute path to bed file
        abs_bed_file = os.path.abspath(bed_file)

        # Create BedTool object
        bed = BedTool(abs_bed_file)

        # Try to peek at the file contents for diagnostics
        try:
            with open(bed_file, 'r') as f:
                first_lines = [next(f) for _ in range(3)]
                logger.info(f"First few lines of {kit_name}: {first_lines}")
        except Exception as e:
            logger.warning(f"Could not read sample from file: {str(e)}")

        # First check if the bed file has any intervals
        if bed.count() == 0:
            logger.error(f"BED file has no valid intervals: {bed_file}")
            return None, None

        # Try to calculate basic statistics using a more robust approach
        try:
            regions = bed.to_dataframe()
            if regions.empty:
                logger.error(f"BED file converted to empty dataframe: {bed_file}")
                return None, None

            regions_count = len(regions)

            # Calculate region sizes
            regions['size'] = regions['end'] - regions['start']
            total_bases = regions['size'].sum()

            # Calculate interval statistics
            min_interval = regions['size'].min() if not regions.empty else 0
            max_interval = regions['size'].max() if not regions.empty else 0
            mean_interval = regions['size'].mean() if not regions.empty else 0

            # Count regions per chromosome
            chrom_counts = {}
            for chrom in regions['chrom'].unique():
                chrom_counts[chrom] = len(regions[regions['chrom'] == chrom])

        except Exception as e:
            logger.error(f"Error processing regions dataframe for {kit_name}: {str(e)}")

            # Fall back to direct BedTools operations without pandas
            logger.info(f"Trying alternative approach for {kit_name}")

            # Get region count directly
            regions_count = bed.count()

            # Calculate total bases directly
            total_bases = sum(end - start for _, start, end in bed)

            # Can't easily get min/max/mean directly, so use defaults
            min_interval = 0
            max_interval = 0
            mean_interval = 0

            # Count regions per chromosome - using BedTools directly for robustness
            chrom_counts = {}
            for interval in bed:
                chrom = interval.chrom
                if chrom in chrom_counts:
                    chrom_counts[chrom] += 1
                else:
                    chrom_counts[chrom] = 1

        # Generate complement BED file (off-target regions)
        complement_file = os.path.join(output_dir, f"{kit_name}_complement.bed")
        complement = bed.complement(g=genome_file)
        complement.saveas(complement_file)

        # Store results
        kit_data = {
            'bed_file': abs_bed_file,
            'complement_file': os.path.abspath(complement_file),
            'regions_count': regions_count,
            'total_bases': int(total_bases),
            'min_interval': int(min_interval),
            'max_interval': int(max_interval),
            'mean_interval': float(mean_interval),
            'chrom_counts': chrom_counts
        }

        return kit_name, kit_data

    except Exception as e:
        logger.error(f"Error processing {bed_file}: {str(e)}")
        # Print traceback for better debugging
        import traceback
        logger.error(traceback.format_exc())
        return None, None


def precompute_kit_data(kit_dir, genome_file, output_file, threads=1):
    """Precompute static information about capture kits to speed up future analyses"""
    logger.info(f"Starting precomputation of capture kit data from {kit_dir}")

    bed_files = [os.path.join(kit_dir, f) for f in os.listdir(kit_dir)
                 if f.endswith('.bed')]

    if not bed_files:
        logger.error(f"No BED files found in {kit_dir}")
        return False

    logger.info(f"Found {len(bed_files)} capture kit BED files")

    precomputed_data = {}

    # We'll use the provided genome file directly
    logger.info(f"Using provided genome file: {genome_file}")

    logger.info(f"Processing {len(bed_files)} BED files")

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(os.path.abspath(output_file))
    os.makedirs(output_dir, exist_ok=True)

    # Process bed files in parallel
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # Create arguments for each bed file: (bed_file, genome_file, output_dir)
        process_args = [(bed_file, genome_file, output_dir) for bed_file in bed_files]
        results = list(executor.map(process_bed_file, process_args))

    # Collect results
    for kit_name, kit_data in results:
        if kit_name is not None:
            precomputed_data[kit_name] = kit_data

    # Save precomputed data
    with open(output_file, 'w') as f:
        json.dump(precomputed_data, f, indent=2)

    logger.info(f"Precomputed data for {len(precomputed_data)} kits saved to {output_file}")
    return True


def calculate_coverage_metrics(bam_file, bed_file, min_mapping_quality=20, total_bases=None):
    """Calculate coverage metrics for a BAM file against a BED file with robust error handling"""
    kit_name = os.path.basename(bed_file).replace('.bed', '')
    logger.info(f"Calculating coverage metrics for kit: {kit_name}")

    # Convert to BedTool objects
    try:
        bam = BedTool(bam_file)
        bed = BedTool(bed_file)

        # First check if bed file has content
        if bed.count() == 0:
            logger.error(f"BED file appears to be empty: {bed_file}")
            return {'kit_name': kit_name, 'error': 'Empty BED file'}

        # Get coverage stats
        coverage = bed.coverage(bam, d=True)

        # Check if coverage file exists and has content
        if not os.path.exists(coverage.fn):
            logger.error(f"Coverage file was not created: {coverage.fn}")
            return {'kit_name': kit_name, 'error': 'Coverage file not created'}

        file_size = os.path.getsize(coverage.fn)
        if file_size == 0:
            logger.error(f"Coverage file is empty: {coverage.fn}")
            return {'kit_name': kit_name, 'error': 'Empty coverage file'}

        # Examine the raw content of the file
        try:
            with open(coverage.fn, 'r') as f:
                sample_lines = [next(f) for _ in range(min(5, sum(1 for _ in open(coverage.fn))))]
            logger.info(f"Coverage file sample for {kit_name}: {sample_lines}")
        except Exception as e:
            logger.warning(f"Could not read coverage file sample: {str(e)}")

        # Try using pandas first (which is what BedTools works with natively)
        try:
            logger.info(f"Attempting to read coverage with pandas for {kit_name}")
            import pandas as pd
            df_pd = pd.read_csv(coverage.fn, sep='\t', header=None,
                                names=['chrom', 'start', 'end', 'depth'])

            if df_pd.empty:
                logger.error(f"Pandas read resulted in empty dataframe for {kit_name}")
            else:
                logger.info(f"Successfully read coverage with pandas: {len(df_pd)} rows")

                # Calculate with pandas first
                df_pd['region_size'] = df_pd['end'] - df_pd['start']

                if total_bases is None:
                    total_bases = df_pd['region_size'].sum()

                covered_bases = df_pd[df_pd['depth'] > 0]['region_size'].sum()
                percent_covered = (covered_bases / total_bases) * 100 if total_bases > 0 else 0

                weighted_mean = (df_pd['depth'] * df_pd['region_size']).sum() / df_pd['region_size'].sum() if df_pd[
                                                                                                                  'region_size'].sum() > 0 else 0

                # Calculate uniformity
                if weighted_mean > 0:
                    weighted_variance = ((df_pd['depth'] - weighted_mean) ** 2 * df_pd['region_size']).sum() / df_pd[
                        'region_size'].sum()
                    cv = np.sqrt(weighted_variance) / weighted_mean
                    uniformity_score = 1 - cv
                else:
                    uniformity_score = 0

                # Calculate coverage thresholds
                coverage_thresholds = [1, 10, 20, 30, 50, 100]
                coverage_percentages = {}

                for threshold in coverage_thresholds:
                    bases_above = df_pd[df_pd['depth'] >= threshold]['region_size'].sum()
                    percentage = (bases_above / total_bases) * 100 if total_bases > 0 else 0
                    coverage_percentages[f'pct_above_{threshold}x'] = percentage

                # Return results
                results = {
                    'kit_name': kit_name,
                    'total_target_bases': total_bases,
                    'covered_bases': covered_bases,
                    'percent_covered': percent_covered,
                    'mean_coverage': weighted_mean,
                    'uniformity_score': uniformity_score
                }
                results.update(coverage_percentages)

                return results

        except Exception as e:
            logger.warning(f"Pandas approach failed for {kit_name}: {str(e)}")

        # If pandas approach failed, try with polars
        try:
            logger.info(f"Attempting to read coverage with polars for {kit_name}")
            # Read file as plain text first to inspect
            with open(coverage.fn, 'r') as f:
                lines = f.readlines()
                if not lines:
                    logger.error(f"Coverage file has no lines: {coverage.fn}")
                    return {'kit_name': kit_name, 'error': 'No content in coverage file'}

                # For diagnostic purposes, examine the first line
                first_line = lines[0].strip()
                logger.info(f"First line in coverage file: '{first_line}'")
                fields = first_line.split('\t')
                logger.info(f"Field count: {len(fields)}, Fields: {fields}")

            # Now try polars with identified schema
            df_columns = ['chrom', 'start', 'end', 'depth']
            df = pl.read_csv(
                coverage.fn,
                separator='\t',
                has_header=False,
                new_columns=df_columns,
                schema_overrides={
                    'chrom': pl.Utf8,
                    'start': pl.Int64,
                    'end': pl.Int64,
                    'depth': pl.Float64
                },
                null_values=["", "NA", "N/A", "null", "None"]
            )

            if df.height == 0:
                logger.error(f"Polars read resulted in empty dataframe for {kit_name}")
                return {'kit_name': kit_name, 'error': 'Empty dataframe from coverage file'}

            logger.info(f"Successfully read coverage with polars: {df.height} rows")

            # Continue with the original polars processing
            df = df.with_column((pl.col("end") - pl.col("start")).alias("region_size"))

            # Calculate basic metrics - use Polars' efficient aggregations
            if total_bases is None:
                total_bases = df.select(pl.sum("region_size")).item()

            # Calculate covered bases (depth > 0)
            covered_bases = df.filter(pl.col("depth") > 0).select(pl.sum("region_size")).item()
            percent_covered = (covered_bases / total_bases) * 100 if total_bases > 0 else 0

            # Calculate mean coverage (weighted by region size)
            mean_coverage = df.select(
                (pl.sum(pl.col("depth") * pl.col("region_size")) / pl.sum("region_size")).alias("mean_cov")
            ).item()

            # Calculate coverage uniformity (coefficient of variation)
            if mean_coverage > 0:
                # Calculate weighted variance
                weighted_variance = df.select(
                    (pl.sum((pl.col("depth") - mean_coverage) ** 2 * pl.col("region_size")) /
                     pl.sum("region_size")).alias("weighted_var")
                ).item()

                cv = np.sqrt(weighted_variance) / mean_coverage
                uniformity_score = 1 - cv  # Higher is better
            else:
                uniformity_score = 0

            # Calculate percentages of bases with coverage above thresholds
            coverage_thresholds = [1, 10, 20, 30, 50, 100]
            coverage_percentages = {}

            for threshold in coverage_thresholds:
                bases_above_threshold = df.filter(pl.col("depth") >= threshold).select(pl.sum("region_size")).item()
                percentage = (bases_above_threshold / total_bases) * 100 if total_bases > 0 else 0
                coverage_percentages[f'pct_above_{threshold}x'] = percentage

            # Calculate results
            results = {
                'kit_name': kit_name,
                'total_target_bases': total_bases,
                'covered_bases': covered_bases,
                'percent_covered': percent_covered,
                'mean_coverage': mean_coverage,
                'uniformity_score': uniformity_score
            }

            # Add coverage threshold percentages
            results.update(coverage_percentages)

            return results

        except Exception as e:
            logger.error(f"Polars approach also failed for {kit_name}: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())

            # Final fallback - direct calculation
            logger.info(f"Attempting direct calculation for {kit_name}")
            try:
                # Reset bed and get direct coverage without intermediate file
                bed = BedTool(bed_file)
                cov_summary = bed.coverage(bam)

                # Process using simple summing operations
                covered_regions = 0
                total_weighted_cov = 0
                total_size = 0

                for line in open(cov_summary.fn):
                    fields = line.strip().split('\t')
                    if len(fields) >= 7:  # coverage produces 7 fields
                        try:
                            start, end = int(fields[1]), int(fields[2])
                            size = end - start
                            cov = float(fields[3])
                            frac = float(fields[6])

                            total_size += size
                            total_weighted_cov += cov * size
                            if frac > 0:
                                covered_regions += size
                        except (ValueError, IndexError):
                            continue

                if total_bases is None:
                    total_bases = total_size

                if total_size > 0:
                    mean_cov = total_weighted_cov / total_size
                    percent_covered = (covered_regions / total_size) * 100
                else:
                    mean_cov = 0
                    percent_covered = 0

                # Create dummy values for coverage thresholds
                coverage_percentages = {}
                for threshold in [1, 10, 20, 30, 50, 100]:
                    coverage_percentages[f'pct_above_{threshold}x'] = 0

                # Return basic results
                basic_results = {
                    'kit_name': kit_name,
                    'total_target_bases': total_bases,
                    'covered_bases': covered_regions,
                    'percent_covered': percent_covered,
                    'mean_coverage': mean_cov,
                    'uniformity_score': 0.5  # Cannot calculate properly, use middle value
                }
                basic_results.update(coverage_percentages)
                return basic_results

            except Exception as e:
                logger.error(f"All approaches failed for {kit_name}: {str(e)}")
                import traceback
                logger.error(traceback.format_exc())

                # Return minimal results
                minimal_results = {
                    'kit_name': kit_name,
                    'error': str(e),
                    'total_target_bases': total_bases if total_bases else 0,
                    'covered_bases': 0,
                    'percent_covered': 0,
                    'mean_coverage': 0,
                    'uniformity_score': 0
                }
                for threshold in [1, 10, 20, 30, 50, 100]:
                    minimal_results[f'pct_above_{threshold}x'] = 0
                return minimal_results

    except Exception as outer_e:
        logger.error(f"Outer exception in coverage calculation for {kit_name}: {str(outer_e)}")
        import traceback
        logger.error(traceback.format_exc())

        # Return minimal results with error
        minimal_results = {
            'kit_name': kit_name,
            'error': str(outer_e),
            'total_target_bases': total_bases if total_bases else 0,
            'covered_bases': 0,
            'percent_covered': 0,
            'mean_coverage': 0,
            'uniformity_score': 0
        }
        for threshold in [1, 10, 20, 30, 50, 100]:
            minimal_results[f'pct_above_{threshold}x'] = 0
        return minimal_results


def calculate_off_target_metrics(bam_file, complement_file, min_mapping_quality=20):
    """Calculate off-target metrics using a precomputed complement BED file"""
    kit_name = os.path.basename(complement_file).replace('_complement.bed', '')
    logger.info(f"Calculating off-target metrics for kit: {kit_name}")

    try:
        # Calculate coverage on off-target regions
        bam = BedTool(bam_file)
        off_target = BedTool(complement_file)
        off_target_coverage = off_target.coverage(bam)

        # Convert to Polars DataFrame
        df_columns = ['chrom', 'start', 'end', 'coverage', 'bases_with_coverage', 'length', 'fraction_covered']

        # Use Polars to read the file
        df = pl.read_csv(
            off_target_coverage.fn,
            separator='\t',
            has_header=False,
            new_columns=df_columns
        )

        # Calculate basic metrics
        off_target_bases = df.select(pl.sum("length")).item()
        off_target_covered_bases = df.select(pl.sum(pl.col("fraction_covered") * pl.col("length"))).item()
        off_target_percent_covered = (off_target_covered_bases / off_target_bases) * 100 if off_target_bases > 0 else 0

        # Calculate mean off-target coverage
        off_target_mean_coverage = df.select(
            (pl.sum(pl.col("coverage") * pl.col("length")) / pl.sum("length")).alias("mean_cov")
        ).item() if off_target_bases > 0 else 0

        # On/off target ratio will be calculated when combining results

        # Calculate results
        results = {
            'kit_name': kit_name,
            'off_target_bases': off_target_bases,
            'off_target_covered_bases': off_target_covered_bases,
            'off_target_percent_covered': off_target_percent_covered,
            'off_target_mean_coverage': off_target_mean_coverage
        }

    except Exception as e:
        logger.error(f"Error calculating off-target metrics for {kit_name}: {str(e)}")
        results = {
            'kit_name': kit_name,
            'off_target_bases': 0,
            'off_target_covered_bases': 0,
            'off_target_percent_covered': 0,
            'off_target_mean_coverage': 0
        }

    return results


def process_kit(args):
    """Process a single kit (for parallel execution)"""
    bam_file, bed_file, genome_file, min_mapping_quality = args
    kit_name = os.path.basename(bed_file).replace('.bed', '')

    try:
        # Calculate on-target metrics
        on_target_metrics = calculate_coverage_metrics(bam_file, bed_file, min_mapping_quality)

        # Calculate off-target metrics
        off_target_metrics = calculate_off_target_metrics(bam_file, bed_file, genome_file, min_mapping_quality)

        # Combine results
        combined_results = {**on_target_metrics, **off_target_metrics}

        # Calculate on/off ratio
        if off_target_metrics['off_target_mean_coverage'] > 0:
            combined_results['on_off_ratio'] = on_target_metrics['mean_coverage'] / off_target_metrics[
                'off_target_mean_coverage']
        else:
            combined_results['on_off_ratio'] = float('inf')

        return combined_results

    except Exception as e:
        logger.error(f"Error processing kit {kit_name}: {str(e)}")
        return {
            'kit_name': kit_name,
            'error': str(e)
        }


def process_kit_precomputed(args):
    """Process a single kit using precomputed data"""
    bam_file, kit_name, kit_data, min_mapping_quality = args

    try:
        # Extract precomputed data
        bed_file = kit_data['bed_file']
        complement_file = kit_data['complement_file']
        total_bases = kit_data['total_bases']

        # Verify that complement file exists and has content
        if not os.path.exists(complement_file):
            logger.error(f"Complement file does not exist for {kit_name}: {complement_file}")
            return {
                'kit_name': kit_name,
                'error': 'Missing complement file',
                'total_target_bases': total_bases,
                'covered_bases': 0,
                'percent_covered': 0,
                'mean_coverage': 0,
                'uniformity_score': 0,
                'on_off_ratio': 0
            }

        # Check complement file size
        if os.path.getsize(complement_file) == 0:
            logger.error(f"Complement file is empty for {kit_name}: {complement_file}")
            return {
                'kit_name': kit_name,
                'error': 'Empty complement file',
                'total_target_bases': total_bases,
                'covered_bases': 0,
                'percent_covered': 0,
                'mean_coverage': 0,
                'uniformity_score': 0,
                'on_off_ratio': 0
            }

        # Try to examine chromosome naming in BAM and BED
        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
            bam_chroms = set(bam.references)
            bam.close()

            with open(bed_file, 'r') as f:
                bed_sample = [next(f).strip().split('\t')[0] for _ in range(min(5, sum(1 for _ in open(bed_file))))]
            bed_chroms = set(bed_sample)

            has_common_chroms = any(b in bam_chroms for b in bed_chroms)

            if not has_common_chroms:
                logger.warning(f"No common chromosome names between BAM and BED for {kit_name}")
                logger.warning(f"BAM chroms (sample): {list(bam_chroms)[:5]}")
                logger.warning(f"BED chroms (sample): {bed_chroms}")
        except Exception as e:
            logger.warning(f"Could not check chromosome naming: {str(e)}")

        # Calculate on-target metrics
        on_target_metrics = calculate_coverage_metrics(
            bam_file, bed_file, min_mapping_quality, total_bases
        )

        # Calculate off-target metrics
        off_target_metrics = calculate_off_target_metrics(
            bam_file, complement_file, min_mapping_quality
        )

        # Combine results
        combined_results = {**on_target_metrics, **off_target_metrics}

        # Calculate on/off ratio
        if off_target_metrics['off_target_mean_coverage'] > 0:
            combined_results['on_off_ratio'] = on_target_metrics['mean_coverage'] / off_target_metrics[
                'off_target_mean_coverage']
        else:
            combined_results['on_off_ratio'] = float('inf')

        return combined_results

    except Exception as e:
        logger.error(f"Error processing kit {kit_name}: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return {
            'kit_name': kit_name,
            'error': str(e)
        }

def generate_tsv_output(results, output_prefix, expected_kit=None):
    """Generate a simple TSV file with kit comparison results for pipeline integration"""
    logger.info("Generating TSV output file")

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not results or len(results) == 0:
        logger.warning("No results to output")
        # Create empty file
        with open(f"{output_prefix}.tsv", 'w') as f:
            f.write("No capture kits met the criteria\n")
        return "VALIDATION_FAILED"

    # Convert list of dicts to Polars DataFrame
    try:
        results_df = pl.from_dicts(results)
    except Exception as e:
        logger.error(f"Error creating DataFrame from results: {str(e)}")
        # Try to recover by converting to a more compatible format
        filtered_results = []
        for r in results:
            filtered_r = {k: (float('nan') if v == float('inf') else v) for k, v in r.items()}
            filtered_results.append(filtered_r)
        results_df = pl.from_dicts(filtered_results)

    # Sort results by percent covered and on/off ratio
    results_df = results_df.sort(
        by=[pl.col("percent_covered").rank(descending=True),
            pl.col("on_off_ratio").rank(descending=True)]
    )

    # Get the identified kit name (top ranked kit)
    identified_kit = results_df.select("kit_name").row(0)[0]

    # Select and rename columns for clarity in the output
    output_columns = [
        'kit_name',
        'percent_covered',
        'mean_coverage',
        'uniformity_score',
        'on_off_ratio'
    ]

    # Add all coverage threshold columns
    coverage_cols = [col for col in results_df.columns if col.startswith('pct_above_')]
    output_columns.extend(coverage_cols)

    # Select columns that exist in the dataframe
    valid_columns = [col for col in output_columns if col in results_df.columns]

    # Create output dataframe with selected columns
    output_df = results_df.select(valid_columns)

    # Add rank column at the beginning
    output_df = output_df.with_row_index(name="rank", offset=1)

    # Add a column indicating if this is the best match
    output_df = output_df.with_column(
        (pl.col("rank") == 1).cast(pl.Int8).alias("best_match")
    )

    # Handle expected kit validation if provided
    kit_validation_status = "UNKNOWN"
    if expected_kit:
        # Check if the expected kit exactly matches the identified kit
        is_expected_kit_exact_match = (expected_kit.lower() == identified_kit.lower())

        # Check if the expected kit is a partial match (substring)
        is_expected_kit_partial_match = (expected_kit.lower() in identified_kit.lower() or
                                         identified_kit.lower() in expected_kit.lower())

        # Add validation columns
        output_df = output_df.with_columns([
            pl.lit(expected_kit).alias("expected_kit"),
            (pl.col("kit_name").str.to_lowercase() == expected_kit.lower()).cast(pl.Int8).alias("exact_match"),
            ((pl.col("kit_name").str.to_lowercase().str.contains(expected_kit.lower())) |
             (pl.lit(expected_kit.lower()).str.contains(pl.col("kit_name").str.to_lowercase()))).cast(pl.Int8).alias(
                "partial_match")
        ])

        # Determine validation status
        if is_expected_kit_exact_match:
            kit_validation_status = "EXACT_MATCH"
        elif is_expected_kit_partial_match:
            kit_validation_status = "PARTIAL_MATCH"
        else:
            kit_validation_status = "NO_MATCH"

        # Add validation status as a column to all rows
        output_df = output_df.with_column(
            pl.lit(kit_validation_status).alias("validation_status")
        )

        logger.info(
            f"Kit validation status: {kit_validation_status} (Expected: {expected_kit}, Identified: {identified_kit})")

    # Save to TSV file
    output_file = f"{output_prefix}.tsv"
    output_df.write_csv(output_file, separator='\t')

    logger.info(f"Results saved to {output_file}")

    return kit_validation_status


def main():
    """Main function"""
    start_time = time.time()

    # Parse command line arguments
    args = parse_arguments()

    # Set the number of threads for Polars if possible
    # This is compatible with different Polars versions
    try:
        # Try newer Polars API first
        try:
            import polars as pl
            pl.Config.set_threading_enabled(True)
            logger.info(f"Set Polars threading to enabled")
        except AttributeError:
            # Try alternative API for older Polars versions
            try:
                pl.set_threading_enabled(True)
                logger.info(f"Set Polars threading to enabled (using alternative API)")
            except:
                pass

        # Try to set number of threads if that API exists
        try:
            pl.ThreadPool.set_num_threads(args.threads)
            logger.info(f"Set Polars to use {args.threads} threads")
        except (AttributeError, NameError):
            logger.info(f"Polars thread count configuration not available in this version")

    except Exception as e:
        logger.warning(f"Could not configure Polars threading: {str(e)}")

    # Check if we're running in precomputation mode
    if args.precompute_mode:
        if not args.kit_dir or not args.precompute_output or not args.genome:
            logger.error("Precomputation mode requires --kit_dir, --precompute_output, and --genome arguments")
            return 1

        # Run precomputation
        success = precompute_kit_data(
            args.kit_dir,
            args.genome,
            args.precompute_output,
            args.threads
        )

        if success:
            logger.info("Precomputation completed successfully")
            return 0
        else:
            logger.error("Precomputation failed")
            return 1

    # Analysis mode
    if not args.bam or not args.precomputed:
        logger.error("Analysis mode requires --bam and --precomputed arguments")
        return 1

    logger.info(f"Starting analysis of {args.bam} using precomputed kit data from {args.precomputed}")

    # Get expected kit information if provided directly
    expected_kit = args.expected_kit
    if expected_kit:
        logger.info(f"Expected capture kit: {expected_kit}")

    # Check if BAM index exists
    check_bam_index(args.bam)

    # Load precomputed kit data
    try:
        with open(args.precomputed, 'r') as f:
            precomputed_data = json.load(f)

        logger.info(f"Loaded data for {len(precomputed_data)} precomputed capture kits")

        if not precomputed_data:
            logger.error("No capture kit data found in precomputed file")
            with open(f"{args.output}.tsv", 'w') as f:
                f.write("No capture kit data found\n")
            return 1

    except Exception as e:
        logger.error(f"Error loading precomputed data: {str(e)}")
        with open(f"{args.output}.tsv", 'w') as f:
            f.write(f"Error loading precomputed kit data: {str(e)}\n")
        return 1

    # Process each kit in parallel
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        process_args = [(args.bam, kit_name, kit_data, args.min_mapping_quality)
                        for kit_name, kit_data in precomputed_data.items()]
        results = list(executor.map(process_kit_precomputed, process_args))

    # Filter out results with errors
    valid_results = [r for r in results if 'error' not in r]
    error_kits = [r['kit_name'] for r in results if 'error' in r]

    if error_kits:
        logger.warning(f"The following kits had errors: {', '.join(error_kits)}")

    # Filter kits based on minimum thresholds
    if valid_results:
        initial_count = len(valid_results)
        valid_results = [r for r in valid_results
                         if r['mean_coverage'] >= args.min_coverage
                         and r['percent_covered'] >= args.min_percent_covered]

        filtered_count = initial_count - len(valid_results)
        if filtered_count > 0:
            logger.info(f"Filtered out {filtered_count} kits that didn't meet minimum coverage thresholds")

    # Generate TSV output file with validation against expected kit
    validation_status = generate_tsv_output(valid_results, args.output, expected_kit)

    # Print most likely kit
    if valid_results:
        # Sort by percent covered and on/off ratio
        sorted_results = sorted(valid_results,
                                key=lambda x: (-x['percent_covered'], -x['on_off_ratio']))
        best_kit = sorted_results[0]

        logger.info(f"Most likely capture kit: {best_kit['kit_name']}")
        logger.info(f"Percent covered: {best_kit['percent_covered']:.2f}%")
        logger.info(f"Mean coverage: {best_kit['mean_coverage']:.2f}x")
        logger.info(f"On/off ratio: {best_kit['on_off_ratio']:.2f}")

        # Print validation status if we had an expected kit
        if expected_kit:
            logger.info(f"Validation status: {validation_status}")
            logger.info(f"Expected kit: {expected_kit}")
            logger.info(f"Identified kit: {best_kit['kit_name']}")
    else:
        logger.warning("No kits passed the filtering criteria or all had errors")
        # Write empty result to TSV file
        with open(f"{args.output}.tsv", 'w') as f:
            f.write("No kits passed filtering criteria\n")

    # Print execution time
    end_time = time.time()
    logger.info(f"Analysis completed in {end_time - start_time:.2f} seconds")

    return 0


if __name__ == "__main__":
    exit(main())