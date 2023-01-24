import os
import csv
from collections import OrderedDict

import read_samples
import utils

chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
main_chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
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
