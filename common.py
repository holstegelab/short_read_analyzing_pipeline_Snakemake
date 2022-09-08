import os
import csv
from collections import OrderedDict

import read_samples
import utils

chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
main_chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
chr_p = [str('0') + str(e) for e in range(0,10)] + [str(i) for i in range(10,90)] + [str(a) for a in range(9000,9009)]
main_chrs_db = []
main_chrs_db.extend(['chr1']*8)
main_chrs_db.extend(['chr2']*7)
main_chrs_db.extend(['chr3']*6)
main_chrs_db.extend(['chr4']*6)
main_chrs_db.extend(['chr5']*6)
main_chrs_db.extend(['chr6']*5)
main_chrs_db.extend(['chr7']*5)
main_chrs_db.extend(['chr8']*5)
main_chrs_db.extend(['chr9']*4)
main_chrs_db.extend(['chr10']*4)
main_chrs_db.extend(['chr11']*4)
main_chrs_db.extend(['chr12']*4)
main_chrs_db.extend(['chr13']*4)
main_chrs_db.extend(['chr14']*4)
main_chrs_db.extend(['chr15']*3)
main_chrs_db.extend(['chr16']*3)
main_chrs_db.extend(['chr17']*3)
main_chrs_db.extend(['chr18']*3)
main_chrs_db.extend(['chr19']*2)
main_chrs_db.extend(['chr20']*2)
main_chrs_db.extend(['chr21']*2)
main_chrs_db.extend(['chr22']*2)
main_chrs_db.extend(['chrX']*5)
main_chrs_db.extend(['chrY']*2)

SAMPLEFILE_TO_SAMPLES = {}

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
        SAMPLEINFO = {}

        #read in all tsv files in current workdir as samplefiles
        for f in os.listdir(filedir):
            if f.endswith('.tsv') and not f.startswith('sample_'): #possible sample file
                f = os.path.join(filedir, f)
                with open(f, 'r', encoding='utf-8') as fopen:
                    ncount = len(fopen.readline().split('\t'))

                if ncount == 8 or ncount == 9: #sample file
                    basename = os.path.splitext(os.path.basename(f))[0]
                    SAMPLE_FILES.append(basename)
                    w = samplefile(f, config)
                    
                    #generate some indices
                    for sample,info in w.items():
                        if sample in SAMPLEINFO:
                            print('WARNING!: Sample ' + sample + ' is defined in more than one sample files.')
                        SAMPLEINFO[sample] = info
                   
                    for key,value in w.items():
                        if len(value['readgroups']) == 0:
                            print('WARNING: %s has no readgroups' % key)


        SAMPLE_FILES.sort()
        config['SAMPLE_FILES'] = SAMPLE_FILES
        config['SAMPLEFILE_TO_SAMPLES'] = SAMPLEFILE_TO_SAMPLES
        config['SAMPLEINFO'] = SAMPLEINFO

    return (config['SAMPLE_FILES'], config['SAMPLEFILE_TO_SAMPLES'], config['SAMPLEINFO'])
