import utils
import yaml
import csv

from constants import *
from read_samples import *
from pathlib import Path
import functools

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
