import pickle
import os
import read_samples
import csv


#number of splits optimized for exome coverage
NSPLITS = dict([('1',81),('10',37),('11',42),('12',43),('13',17),('14',26),('15',31),('16',33),('17',41),('18',14),('19',38),('2',63),('20',20),('21',10),('22',19),('3',49),('4',35),('5',38),('6',40),('7',43),('8',30),('9',35),('MT',3),('X',29),('Y',4)])

NSPLITS_V2={'1': 799,
 '10': 331,
 '11': 430,
 '12': 424,
 '13': 144,
 '14': 251,
 '15': 293,
 '16': 318,
 '17': 416,
 '18': 127,
 '19': 403,
 '2': 591,
 '20': 192,
 '21': 87,
 '22': 177,
 '3': 467,
 '4': 323,
 '5': 371,
 '6': 391,
 '7': 396,
 '8': 282,
 '9': 326,
 'MT': 1,
 'X': 290,
 'Y': 30}

LINES_PER_SPLIT=250

PCR_MODEL_DEFAULT='CONSERVATIVE'   
PCR_MODEL_PCRFREE='NONE'

SAMPLE_CACHE = {}
SAMPLEINFO = {}


#HOW TO USE:
#1) cd to directory with sample file
#2) snakemake --snakefile 'path to this file'   <target>  
#   where <target> is equal to:
#   - empty: all sample files in this directory
#   - 'samplefile.done_v1': samples in samplefile.tsv
#   - 'samplefile.bqsr_done'
#   - bam_v1/<sample>.bam, gvcf_v1/<sample>.gvcf.gz, etc. 


def get_contamination_from_quality_file(filename):
    res = dict()
    with open(filename, 'r') as f:
        c = csv.reader(f,delimiter='\t')
        for row in c:
            res[row[0]] = (row[-2], row[-1])
    return res

#function to read in and cache a samplefile
def samplefile(filename):
   basename = os.path.basename(filename)
   if not basename in SAMPLE_CACHE: 
        datfilename = os.path.realpath(filename)[:-4] + '.dat'
        if not os.path.isfile(datfilename):
            if not os.path.isfile(filename):
                raise RuntimeError('Sample file ' + filename + ' does not exist')
            info = read_samples.read_samplefile(filename)
            sampleinfodict = dict([(a['sample'], a) for a in info])
            with open(datfilename,'wb') as f:
                pickle.dump(sampleinfodict, f)
        else:
            with open(datfilename, 'rb') as f:
                sampleinfodict = pickle.load(f)
            
        SAMPLE_CACHE[basename] = sampleinfodict
   return SAMPLE_CACHE[basename]


SAMPLE_FILES = []
SAMPLES_BY_FILE = {}

#read in all tsv files in current workdir as samplefiles
for f in os.listdir('.'):
    if f.endswith('.tsv') and not f.startswith('sample_'): #possible sample file
        f = os.path.join(os.getcwd(), f)
        with open(f, 'r') as fopen:
            ncount = len(fopen.readline().split('\t'))

        if ncount == 8 or ncount == 9: #sample file
            SAMPLE_FILES.append(f)
            w = samplefile(f)
            nw = {}
            for key,value in w.items():
                if len(value['readgroups']) > 0:
                    nw[key] = value
                else:
                    print('WARNING: %s has no readgroups' % key)

                    nw[key] = value
            SAMPLEINFO.update(nw)
            SAMPLES_BY_FILE[os.path.basename(f)] = nw


SAMPLE_FILES.sort()