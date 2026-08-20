import utils
import yaml
import csv
import os
import math
import zlib
import time
import subprocess
import sys
import tempfile
import shutil
from shlex import quote

from constants import *
from read_samples import *
from pathlib import Path
import functools
from snakemake.shell import shell


def zslurm_lease_command(workflow_config):
    """Return the pipeline-local ZSlurm lease client by default.

    Fused rules activate different Conda environments, so relying on whichever
    ``zslurm_lease`` happens to be on the submit host's PATH is not portable.
    Keep the client with the workflow and pass that shared absolute path into
    every fused runner. An explicit config/environment override remains
    available for development and protocol compatibility tests.
    """
    pipeline_default = (
        Path(__file__).resolve().parent / 'scripts' / 'zslurm_lease_client.py'
    )
    configured = workflow_config.get(
        'zslurm_lease_command',
        os.environ.get('ZSLURM_LEASE_COMMAND', str(pipeline_default)),
    )
    configured = os.path.expanduser(str(configured))
    if os.sep in configured and not os.path.isabs(configured):
        configured = str(Path(__file__).resolve().parent / configured)
    resolved = (
        configured
        if os.sep in configured
        else shutil.which(configured)
    )
    if resolved and os.path.isfile(resolved) and os.access(resolved, os.X_OK):
        return os.path.realpath(resolved)
    return configured


# --- Active-storage reservation: shared helper + staged release --------------
# active_use_gb() reserves the input+bam PEAK at start_sample and is the shared basis
# for the staged RELEASE of the reservation, and lives here in common so every module
# (Aligner markdup, Encrypt copy_to_dcache, Snakefile finished_sample) can reach it.
#
# The peak collapses in stages, so we hand the reservation back in stages instead of
# holding the full peak until finished_sample:
#   * markdup done  -> fastqs + intermediate bams gone (only markdup.bam remains)
#   * cram uploaded -> mapped_hg38.cram gone
#   * finished      -> markdup.bam gone (deepvariant/stats/whatshap/chrM done)
# Fractions are deliberately CONSERVATIVE: releasing too much risks a real active
# disk overflow; releasing too little only over-reserves (safe). Tune here. They MUST
# sum to <= 1.0; a sample must traverse markdup + copy_to_dcache for the per-sample
# add/remove to balance exactly (partial paths only ever under-release = safe).
ACTIVE_RELEASE_FRAC_MARKDUP = 0.30
ACTIVE_RELEASE_FRAC_UPLOAD = 0.15

def active_use_gb(wildcards):
    """Reserve source-input lifecycle plus the estimated processing peak."""
    sample = SAMPLEINFO[wildcards['sample']]
    filesize = sample['filesize']
    capture_kit = sample['capture_kit']
    active_filesize = 2.0 * filesize if 'cram' in sample['file_type'] else filesize
    # Every route owns the input bytes for its whole sample lifecycle.  For an
    # external route these bytes are materialized here; for an active route
    # they already exist but must still be included so staged releases balance
    # the reservation and the scheduler sees the real active-storage pressure.
    res = filesize
    if capture_kit == 'WGS38_to_exome' or capture_kit == 'WGS37_to_exome':
        res += 2.0 * active_filesize * 0.15
    else:
        res += 2.0 * active_filesize
    return res

def active_release_markdup(wildcards):
    """Release the fastq + intermediate-bam share of the reservation at markdup."""
    return ACTIVE_RELEASE_FRAC_MARKDUP * active_use_gb(wildcards)

def active_release_upload(wildcards):
    """Release the cram share of the reservation once the cram has been uploaded."""
    return ACTIVE_RELEASE_FRAC_UPLOAD * active_use_gb(wildcards)

def active_release_finished(wildcards):
    """Release the remainder (markdup.bam share) at finished_sample."""
    return max(0.0, 1.0 - ACTIVE_RELEASE_FRAC_MARKDUP - ACTIVE_RELEASE_FRAC_UPLOAD) * active_use_gb(wildcards)


chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
main_chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
main_chrs_ploidy_male = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrXH','chrYH']
main_chrs_ploidy_female = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']


# Genome split levels (see Tools.smk). 

# 4 different components: A: autosomes + MT, X: X chromosome, Y: Y chromosome, F: full genome
# Autosomes have 5 different levels: 0: no split, 1: 10 splits, 2: 100 splits, 3: 1000 splits, 4: 10000 splits
# X and Y chromosome are split separately in 4 levels (levels: 0, 5, 50, 500 for X, 0, 2, 20, 200 for Y)
# full genome is not split (only level 0)

# 1. Split elements are 4-tuples of (component, level (nr. of digits), splitnr, ploidy)
# 2. A corresponding string region describer is generated from this tuple, e.g. ('A', 2, 3, 2) -> 'A23'
#    or ('X', 1, 3, 1) -> 'X3H' (see get_regions function below)
# 3. String region describers can be used to get the corresponding bed /interval_list files describing 
#    the region, e.g. A23 -> wes_bins_v2/merged.autosplit2.03.bed or wgs_bins/genome.autosplit2.03.bed (see region_to_file function below)
# 4. '.padded.bed' files have 1000bp of padding at the start and end of the bed file added to them, to enable accurate calling 
#    at the junctions between the bed files. Padded files are avilable for level 1-3 for A, X and Y. 




# Note: The split files are not disjoint. The split files of a higher level
# are a subset of the split files of a lower level. E.g. genome.autosplit1.3.bed combines
# genome.autosplit2.30.bed to genome.autosplit2.39.bed, and genome.autosplit2.30.bed
# combines genome.autosplit3.300.bed to genome.autosplit3.399.bed.
# Function convert_to_level0 (see below) can be used to convert a region describer to a level 0 region describer.

# Note 2: exome and wgs splits are synced. That is
# merged.<component>split<level>.<splitnr>.bed is always the exome region within
# genome.<component>split<level>.<splitnr>.bed.
# Note also that due to this at higher levels the exome split files can sporadically be empty.
# also, exome levels only go to level 3 (for auto) and level 2 (for X and Y).

# Note 3: autosomes of level n are usually combined with sex chromosomes of level n-1. This is described
#         by the level0_range, level1_range, level2_range, level3_range etc. lists below.


# Note 4: To deal with some peculiarities of GenomicDBImport, there is also a separate set of
# bam files in which A only covers the classical chromosomes and no other contigs.  These A files have
# the extension .classic.bed.
# All other contigs are stored in O files, and cover the whole contig. To use these, you need to make use
# of  levelX_regions_so lists, and enable the 'classic=True' on the region_to_file function.


#redefine the srcdir function to use the workflow.basedir
#(was removed from snakemake)
import inspect
def srcdir(path):
    
    frame = inspect.currentframe()
    try:
        # Walk back to the caller's frame (1 level up)
        caller_frame = frame.f_back
        # Access the 'workflow' variable from the caller's local variables
        workflow = caller_frame.f_globals.get("workflow")
        if workflow is None:
            raise ValueError("Could not find 'workflow' in the calling frame.")
        return Path(workflow.basedir) / Path(path)    
        # Add logic that uses the workflow variable
    finally:
        # Clean up the frame reference to avoid reference cycles
        del frame
        del caller_frame
    


level0_range = [('F', 0,0,2), ('X',0,0, 1), ('Y', 0,0, 1)]

level1_range = [('A', 1,x,2) for x in range(0,10)] + \
                [('X',0,0, 2), ('X', 0,0, 1),
                 ('Y',0,0, 2), ('Y', 0,0, 1)]

level2_range = [('A', 2,x,2) for x in range(0,100)] + \
               [('X', 1,x,2) for x in range(0,5)] + \
               [('X', 1, x, 1) for x in range(0, 5)] + \
               [('Y', 1, x,2) for x in range(0,2)] + \
               [('Y', 1, x, 1) for x in range(0, 2)] 

level3_range = [('A', 3, x,2) for x in range(0,1000)] + \
               [('X', 2, x,2) for x in range(0,50)] + \
               [('X', 2, x, 1) for x in range(0, 50)] + \
               [('Y', 2, x,2) for x in range(0,20)] + \
               [('Y', 2, x, 1) for x in range(0, 20)]               

level4_range = [('A', 4, x,2) for x in range(0,10000)] + \
                [('X', 3, x,2) for x in range(0,500)] + \
                [('X', 3, x, 1) for x in range(0, 500)] + \
                [('Y', 3, x,2) for x in range(0,200)] + \
                [('Y', 3, x, 1) for x in range(0, 200)]

level0_range_diploid_only = [('F', 0,0,2)]

level1_range_diploid_only = [('A', 1,x,2) for x in range(0,10)]

level2_range_diploid_only = [('A', 2,x,2) for x in range(0,100)] + \
               [('X', 1,x,2) for x in range(0,5)] + \
               [('Y', 1, x,2) for x in range(0,2)]

level3_range_diploid_only = [('A', 3, x,2) for x in range(0,1000)] + \
               [('X', 2, x,2) for x in range(0,50)] + \
               [('Y', 2, x,2) for x in range(0,20)]

level4_range_diploid_only = [('A', 4, x,2) for x in range(0,10000)] + \
                [('X', 3, x,2) for x in range(0,500)] + \
                [('Y', 3, x,2) for x in range(0,200)]
              


#levels were other contigs are separated, and bed files cover full contig. Intended for use with GenomicDBImport.
level2_range_so = [('A', 2,x,2) for x in range(0,99)] + \
               [('X', 1,x,2) for x in range(0,5)] + \
               [('X', 1, x, 1) for x in range(0, 5)] + \
               [('Y', 1, x,2) for x in range(0,2)] + \
               [('Y', 1, x, 1) for x in range(0, 2)] + \
               [('O', 1, x, 2) for x in range(0,2)]

level3_range_so = [('A', 3, x,2) for x in range(0,989)] + \
               [('X', 2, x,2) for x in range(0,50)] + \
               [('X', 2, x, 1) for x in range(0, 50)] + \
               [('Y', 2, x,2) for x in range(0,20)] + \
               [('Y', 2, x, 1) for x in range(0, 20)] + \
               [('O', 2, x, 2) for x in range(0,20)] 

level4_range_so= [('A', 4, x,2) for x in range(0,9895)] + \
                [('X', 3, x,2) for x in range(0,500)] + \
                [('X', 3, x, 1) for x in range(0, 500)] + \
                [('Y', 3, x,2) for x in range(0,200)] + \
                [('Y', 3, x, 1) for x in range(0, 200)] + \
                [('O', 3, x, 2) for x in range(0,200)]

def get_regions(lrange):
    """Converts a region describer (tuple format) to a list of regions (string format).
    E.g. [('A', 1, 3, 1), ('A', 2, 4, 2)] -> ['A3H', 'A04']
    """
    res = []
    
    for component, level,splitnr, ploidy in lrange:
        if ploidy == 1:
            ploidy = 'H'
        else:
            ploidy = ''
        if level == 0:
            region = f'{component}{ploidy}'
        else:
            region = f'{component}{splitnr:0{level}d}{ploidy}'
        res.append(region)
        
    return res

level0_regions = get_regions(level0_range)
level1_regions = get_regions(level1_range)
level2_regions = get_regions(level2_range)
level3_regions = get_regions(level3_range)
level4_regions = get_regions(level4_range)

level1_regions_diploid = get_regions(level1_range_diploid_only)
level2_regions_diploid = get_regions(level2_range_diploid_only)
level3_regions_diploid = get_regions(level3_range_diploid_only)
level4_regions_diploid = get_regions(level4_range_diploid_only)

level2_regions_so = get_regions(level2_range_so)
level3_regions_so = get_regions(level3_range_so)
level4_regions_so = get_regions(level4_range_so)

def convert_to_level0(region):
    """Converts a region describer of level >=0 to the corresponding level 0 region.
    E.g. A33 -> F, X22H -> XH, X1 -> F
    """
    if region in level0_regions:
        return region
    if region.startswith('A') or region.startswith('F') or region.startswith('O'):
        return 'F'
    elif region.startswith('X'):
        return 'F' if not region.endswith('H') else 'XH'
    elif region.startswith('Y'):
        return 'F' if not region.endswith('H') else 'YH'
    else:
        raise ValueError(f'Unknown region {region}')

def convert_to_level1(region):
    """Converts a region describer of level >=1 to the corresponding level 1 region.
    E.g. A33 -> A3, X22H -> X2H, X1 -> F
    """
    # Check if region ends with 'H'
    ends_with_H = region.endswith('H')
    if region in level1_regions:
        return region
    else:
        if ends_with_H:
            region = region[:-1]

        component = region[0]
        level = len(region[1:])
        split = region[1:]

        # If level is greater than 1, truncate the split number to 1 digit
        if component == 'X':
            level1_region = 'X'
        elif component == 'Y':
            level1_region = 'Y'
        elif component == 'O':
            level1_region = 'A9'
        else:
            if level > 1:
                split = split[:1]

            # Construct the level 1 region name
            level1_region = f'{component}{split}'

            # If original region ended with 'H', add it back
        if ends_with_H:
            level1_region += 'H'

    return level1_region




def region_to_file(region, wgs=False, classic=False, padding=False, extension='bed'):
    """ Converts a region describer to the filename of the file describing the region.
    
        E.g. A33H, wgs=False, extension=bed -> <interval_folder>/wes_bins_v2/merged.autosplit2.33.bed
        
        :param region: region describer
        :param wgs: whether to use the wgs or wes intervals
        :param classic: only include classic (chr1-chr22) autosomes. No effect for for O,X,Y components.
        :param padding: add padding (1000bp) to end of interval files.
        :param extension: extension of the file (bed or interval_list)
    """
    
    component = region[0]
    if region.endswith('H'):
        region = region[:-1]

    split = region[1:]
    
    if padding and component in ['A','X','Y'] and split:
        preextension = ['padded']
    else:
        preextension = []

    if component == 'A':
        component = 'auto'
        if classic:
            preextension.append('classic')
    elif component == 'F':
        component = 'full'
        if classic:
            preextension.append('classic')

    extension = '.'.join(preextension + [extension])

    level = len(split)
    if level > 0:
        split = '.' + split
    else:
        split = ''
    if wgs: 
        f = pj(INTERVALS_DIR, f'wgs_bins_v3/genome.{component}split{level}{split}.{extension}')
    else:
        f = pj(INTERVALS_DIR, f'wes_bins_v3/merged.{component}split{level}{split}.{extension}')
    return f
    


# OLD REGIONS

chr_p = [str('0') + str(e) for e in range(0, 10)] + [str(i) for i in range(10, 90)] + [str(a) for a in range(9000, 9762)]
main_chrs_db = []
main_chrs_db.extend(['chr1']*84)
main_chrs_db.extend(['chr2']*66)
main_chrs_db.extend(['chr3']*51)
main_chrs_db.extend(['chr4']*36)
main_chrs_db.extend(['chr5']*40)
main_chrs_db.extend(['chr6']*42)
main_chrs_db.extend(['chr7']*43)
main_chrs_db.extend(['chr8']*30)
main_chrs_db.extend(['chr9']*35)
main_chrs_db.extend(['chr10']*37)
main_chrs_db.extend(['chr11']*45)
main_chrs_db.extend(['chr12']*47)
main_chrs_db.extend(['chr13']*16)
main_chrs_db.extend(['chr14']*27)
main_chrs_db.extend(['chr15']*32)
main_chrs_db.extend(['chr16']*36)
main_chrs_db.extend(['chr17']*45)
main_chrs_db.extend(['chr18']*14)
main_chrs_db.extend(['chr19']*44)
main_chrs_db.extend(['chr20']*21)
main_chrs_db.extend(['chr21']*10)
main_chrs_db.extend(['chr22']*19)
main_chrs_db.extend(['chrX']*30)
main_chrs_db.extend(['chrY']*3)

valid_chr_p = {'chr1': chr_p[:84],
               'chr2': chr_p[84:150],
               'chr3': chr_p[150:201],
               'chr4': chr_p[201:237],
               'chr5': chr_p[237:277],
               'chr6': chr_p[277:319],
               'chr7': chr_p[319:362],
               'chr8': chr_p[362:392],
               'chr9': chr_p[392:427],
               'chr10': chr_p[427:464],
               'chr11': chr_p[464:509],
               'chr12': chr_p[509:556],
               'chr13': chr_p[556:572],
               'chr14': chr_p[572:599],
               'chr15': chr_p[599:631],
               'chr16': chr_p[631:667],
               'chr17': chr_p[667:712],
               'chr18': chr_p[712:726],
               'chr19': chr_p[726:770],
               'chr20': chr_p[770:791],
               'chr21': chr_p[791:801],
               'chr22': chr_p[801:820],
               'chrX': chr_p[820:850],
               'chrY': chr_p[850:]}


@functools.cache
def _read_sex_file(filename):
    sex = "UNK"
    with open(filename) as f:
        for line in f:
            if line.startswith('sex: '):
                parts = line.split(':')
                if len(parts) == 2:
                    sex = parts[1].strip()
                break

    assert sex == 'M' or sex == 'F', 'Unknown sex in sex detection result file.'
    return 'male' if sex == 'M' else 'female'

def get_validated_sex_file(input):
    #this file should exist after running 'get_validated_sex' job.
    #it should also certainly exist after the bam file is created,
    #as it relies on this file.
    filename = input['validated_sex']
    return _read_sex_file(filename)

def get_ref_by_validated_sex(wildcards, input):
    sex = get_validated_sex_file(input)
    return REF_FEMALE if sex == 'female' else REF_MALE

def get_refdir_by_validated_sex(wildcards, input):
    sex = get_validated_sex_file(input)
    return REF_FEMALE_DIR if sex == 'female' else REF_MALE_DIR

def get_strref_by_validated_sex(wildcards, input):
    sex = get_validated_sex_file(input)
    return REF_FEMALE_STR if sex == 'female' else REF_MALE_STR

def input_size_mb(paths):
    """Return the aggregate size of one or more Snakemake inputs in MiB."""
    size_mb = getattr(paths, 'size_mb', None)
    if size_mb is not None:
        if callable(size_mb):
            size_mb = size_mb()
        return float(size_mb)

    if isinstance(paths, (str, os.PathLike)):
        path = os.fspath(paths)
        if os.path.isdir(path):
            total = 0
            for root, _, filenames in os.walk(path):
                total += sum(os.path.getsize(os.path.join(root, name))
                             for name in filenames)
            return total / (1024.0 * 1024.0)
        return os.path.getsize(path) / (1024.0 * 1024.0)

    return sum(input_size_mb(path) for path in paths)


def ssd_gb_for_inputs(paths, factor=1.0, overhead_gb=1.0, minimum_gb=1):
    """Estimate node-local scratch, rounded up to whole GiB.

    ``factor`` describes the temporary-data/input-size ratio. ``overhead_gb``
    covers metadata, indexes and tools which briefly keep an extra small file.
    The ratios are calibrated against live /scratch-node use; keeping this in
    one helper makes future recalibration explicit.
    """
    estimated = (
        input_size_mb(paths) / 1024.0 * float(factor)
        + float(overhead_gb)
    )
    return max(int(minimum_gb), int(math.ceil(estimated)))


def node_ssd_base(tmpdir_fallback=None):
    user = os.environ.get('USER','')
    slurm_tmp = os.environ.get('SLURM_TMPDIR')
    job_id = os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_JOBID')
    if job_id:
        p = f"/scratch-node/{user}.{job_id}"
        if os.path.isdir(p) and os.access(p, os.W_OK):
            return p
    base = "/scratch-node"
    if os.path.isdir(base):
        try:
            entries = [os.path.join(base, d) for d in os.listdir(base) if d.startswith(user + '.')]
        except OSError:
            entries = []
        entries = [e for e in entries if os.path.isdir(e) and os.access(e, os.W_OK)]
        if entries:
            entries.sort(key=lambda q: os.stat(q).st_mtime, reverse=True)
            return entries[0]
    if tmpdir_fallback:
        return tmpdir_fallback
    if slurm_tmp and os.path.isdir(slurm_tmp) and os.access(slurm_tmp, os.W_OK):
        return slurm_tmp
    if os.path.isdir(TMPDIR_ALT) and os.access(TMPDIR_ALT, os.W_OK):
        return os.path.join(TMPDIR_ALT, user)
    return tmpdir

def node_tmp_path(*segments):
    base = node_ssd_base()
    parts = [str(s) for s in segments]
    return os.path.join(base, *parts)

def get_samplefile_folder(samplefile):
    return os.path.dirname(os.path.realpath(samplefile + '.tsv'))


def read_sexchrom(filename):
    result = {}
    with open(filename,'r') as f:
        r = csv.reader(f,delimiter='\t')
        lines = [row for row in r]
        for row in lines[1:]:
            
            result[row[0]] = row[2]
    return result            

cache = {}
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',cache)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()


def external_data_dir(sample, sinfo=None):
    """Return the protocol-specific active-storage directory for a sample."""
    info = SAMPLEINFO[sample] if sinfo is None else sinfo
    suffix = {
        "archive": ".data",
        "dcache": ".dcache_data",
        "s3": ".s3_data",
    }.get(info.get("from_external"), ".data")
    return pj(SOURCEDIR, str(sample) + suffix)


def get_time(rulename):
    """`resources: time = get_time('align_reads')` -- expected wall-clock seconds.

    Reads the measured estimate from constants.RUNTIME and picks the exome or the
    WGS column from this sample's sample_type. Rules without a `sample` wildcard
    (per-samplefile gathers) get the WGS column, which is the longer of the two.
    A retry gets 50% more time per attempt, mirroring the get_mem_mb_* helpers,
    since a job that hit the wall is exactly the one that needs more.
    """
    def _time(wildcards, attempt=1):
        wgs_seconds, exome_seconds = RUNTIME[rulename]
        seconds = wgs_seconds
        sample = getattr(wildcards, 'sample', None)
        if sample is not None and sample in SAMPLEINFO:
            if 'wgs' not in SAMPLEINFO[sample]['sample_type'].lower():
                seconds = exome_seconds
        return int(seconds * (1.0 + 0.5 * (int(attempt) - 1)))
    return _time

def remote_base_for_sample(sample):
    sinfo = SAMPLEINFO[sample]
    target = sinfo['target']
    samplefile = os.path.basename(sinfo['samplefile'])
    if target is None:
        target = os.path.join(sinfo['study'], samplefile)
    if target.endswith('/'):
        target = target[:-1]
    return target


def remote_base_for_samplefile(samplefile):
    samples = list(SAMPLEFILE_TO_SAMPLES[samplefile].keys())
    if not samples:
        raise ValueError(f"Samplefile {samplefile} has no samples defined")
    return remote_base_for_sample(samples[0])


def _dcache_endpoint(value):
    endpoint = parse_dcache_uri(value)
    if endpoint is None:
        return None
    remote, path = endpoint
    config_path = DCACHE_CONFIGS.get(remote)
    if not config_path:
        raise FileNotFoundError(
            f"No macaroon config registered for dCache remote {remote!r}; "
            "put <remote>.conf next to the sample listing"
        )
    return remote, path, config_path


def copy_from_dcache_uri(remote_uri, local_path, *, no_stage=False):
    """Download one dCache URI directly with checksum verification."""
    endpoint = _dcache_endpoint(remote_uri)
    if endpoint is None:
        raise ValueError(f"Not a dCache URI: {remote_uri!r}")
    remote, remote_path, config_path = endpoint
    local_path = os.path.abspath(str(local_path))
    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    fd, file_list = tempfile.mkstemp(prefix=".dcache-download-", suffix=".tsv", dir=os.path.dirname(local_path))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(f"{remote_path}\t{local_path}\n")
        cmd = [
            sys.executable,
            str(Path(__file__).resolve().parent / "scripts" / "dcache_transfer.py"),
            "download",
            "--config",
            str(config_path),
            "--remote",
            str(remote),
            "--file-list",
            file_list,
            "--workers",
            "1",
        ]
        if no_stage:
            cmd.append("--no-stage")
        subprocess.run(cmd, check=True)
    finally:
        try:
            os.unlink(file_list)
        except FileNotFoundError:
            pass


def copy_with_checksum(local_path, remote_dir, remote_name, checksum_path, config_path, ada_script, remote_profile='agh_processed'):
    endpoint = _dcache_endpoint(remote_dir)
    if endpoint is not None:
        remote, bare_remote_dir, endpoint_config = endpoint
        remote_file = os.path.join(bare_remote_dir, remote_name)
        cmd = [
            sys.executable,
            str(Path(__file__).resolve().parent / "scripts" / "dcache_transfer.py"),
            "upload",
            "--config",
            str(endpoint_config),
            "--remote",
            str(remote),
            "--source",
            str(local_path),
            "--destination",
            remote_file,
            "--checksum-output",
            str(checksum_path),
        ]
        subprocess.run(cmd, check=True)
        return

    adler_local = 1
    with open(local_path, 'rb') as fhandle:
        for chunk in iter(lambda: fhandle.read(16 * 1024 * 1024), b''):
            adler_local = zlib.adler32(chunk, adler_local)
    adler_local &= 0xffffffff

    adler_remote = ''
    retries = 0
    remote_dir_full = f"{remote_profile}:{remote_dir}"
    remote_file = f"{remote_dir}/{remote_name}"
    remote_file_full = f"{remote_dir_full}/{remote_name}"
    config_q = quote(config_path)
    local_q = quote(local_path)
    checksum_q = quote(checksum_path)
    remote_dir_full_q = quote(remote_dir_full)
    remote_file_full_q = quote(remote_file_full)
    remote_file_q = quote(remote_file)
    ada_q = quote(str(ada_script))

    shell(f"rclone --config {config_q} mkdir -v {remote_dir_full_q}")

    while f'{adler_local:08x}' != adler_remote and retries <= 3:
        shell(f"rclone --config {config_q} -v copyto {local_q} {remote_file_full_q}")
        shell(f"{ada_q} --tokenfile {config_q} --api https://dcacheview.grid.surfsara.nl:22880/api/v1 --checksum {remote_file_q} | awk '{{{{print $2}}}}' | awk -F '=' '{{{{print $2}}}}' > {checksum_q}")
        with open(checksum_path, 'r') as sum_file:
            adler_remote = sum_file.readline().rstrip('\n')

        retries += 1
        if f'{adler_local:08x}' != adler_remote:
            shell(f"rclone --config {config_q} -v deletefile {remote_file_full_q}")
            time.sleep(60)

    if f'{adler_local:08x}' != adler_remote:
        raise ValueError(f"Checksums do not match for {local_path} after 3 retries. Local: {adler_local:08x}, Remote: {adler_remote}")
