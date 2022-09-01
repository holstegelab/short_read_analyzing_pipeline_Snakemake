import os
import csv
from collections import OrderedDict

import read_samples
import utils

chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
main_chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


SAMPLEFILE_TO_SAMPLES = {}

def get_contamination_from_quality_file(filename):
    res = dict()
    with open(filename, 'r') as f:
        c = csv.reader(f,delimiter='\t')
        for row in c:
            res[row[0]] = (row[-2], row[-1])
    return res


#function to read in and cache a samplefile
def samplefile(sfilename, config):
    basename = os.path.splitext(os.path.basename(sfilename))[0]
    if not basename in SAMPLEFILE_TO_SAMPLES: 
        if not os.path.isfile(sfilename):
            raise RuntimeError('Sample file ' + sfilename + ' does not exist')
       
        datfilename = os.path.realpath(sfilename)[:-4] + '.adat'
        if not os.path.isfile(datfilename):
            if not os.path.isfile(sfilename):
                raise RuntimeError('Sample file ' + sfilename + ' does not exist')
            samplelist = read_samples.read_samplefile_simple(sfilename, config)
            #using ordereddict to main stable order
            sampleinfodict = OrderedDict([(a['sample'], a) for a in samplelist])
            utils.save(sampleinfodict,datfilename)
        else:
            sampleinfodict = utils.load(datfilename)

        SAMPLEFILE_TO_SAMPLES[basename] = sampleinfodict
    return SAMPLEFILE_TO_SAMPLES[basename]


def load_samplefiles(filedir, config):
    if not 'SAMPLE_FILES' in config:
        SAMPLE_FILES = []
        SAMPLES_BY_FILE = {}
        SAMPLEINFO = {}

        #read in all tsv files in current workdir as samplefiles
        for f in os.listdir(filedir):
            if f.endswith('.tsv') and not f.startswith('sample_'): #possible sample file
                f = os.path.join(os.getcwd(), f)
                with open(f, 'r') as fopen:
                    ncount = len(fopen.readline().split('\t'))

                if ncount == 8 or ncount == 9: #sample file
                    SAMPLE_FILES.append(f)
                    w = samplefile(f, config)
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
        config['SAMPLE_FILES'] = SAMPLE_FILES
        config['SAMPLES_BY_FILE'] = SAMPLES_BY_FILE
        config['SAMPLEINFO'] = SAMPLEINFO

    return (config['SAMPLE_FILES'], config['SAMPLES_BY_FILE'], config['SAMPLEINFO'])
