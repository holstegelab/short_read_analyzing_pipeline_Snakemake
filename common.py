import pickle
import os
import read_samples
import csv


chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
main_chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


SAMPLE_CACHE = {}


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

def load_samplefiles(filedir, config):
    if not 'SAMPLE_FILES' in config:
        SAMPLE_FILES = []
        SAMPLES_BY_FILE = {}
        SAMPLEINFO = {}

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
        config['SAMPLE_FILES'] = SAMPLE_FILES
        config['SAMPLES_BY_FILE'] = SAMPLES_BY_FILE
        config['SAMPLEINFO'] = SAMPLEINFO

    return (config['SAMPLE_FILES'], config['SAMPLES_BY_FILE'], config['SAMPLEINFO'])
