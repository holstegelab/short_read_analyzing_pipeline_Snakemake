import os
import csv
from collections import OrderedDict

import read_samples
from read_samples import PROTOCOLS
import utils
import yaml

chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
main_chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
main_chrs_ploidy_male = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrXH','chrYH']
main_chrs_ploidy_female = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
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
MAX_BATCH_SIZE = 2 * 1024 #2TB


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
            samplelist = read_samples.read_samplefile(sfilename, config)
            #using ordereddict to main stable order
            sampleinfodict = OrderedDict([(a['sample'], a) for a in samplelist])
            utils.save(sampleinfodict,datfilename)
        else:
            sampleinfodict = utils.load(datfilename)

        SAMPLEFILE_TO_SAMPLES[basename] = sampleinfodict
    return SAMPLEFILE_TO_SAMPLES[basename]


def load_samplefiles(filedir, config):
    
    if not 'SAMPLE_FILES' in config:
        #already loaded

        SAMPLE_FILES = []
        SAMPLEINFO = {}
        SAMPLE_TO_BATCH = {}
        SAMPLEFILE_TO_BATCHES = {}

        #read in all tsv files in current workdir as samplefiles
        for f in os.listdir(filedir):
            if f.endswith('.tsv') and not f.startswith('sample_'): #possible sample file
                f = os.path.join(filedir, f)
                with open(f, 'r', encoding='utf-8') as fopen:
                    ncount = len(fopen.readline().split('\t'))
                if ncount == 7 or ncount == 8 or ncount == 9: #sample file
                    basename = os.path.splitext(os.path.basename(f))[0]
                    SAMPLE_FILES.append(basename)
                    w = samplefile(f, config)
                    
                    #generate some indices
                    no_readgroup = []
                    for sample,info in w.items():
                        if sample in SAMPLEINFO:
                            print('WARNING!: Sample ' + sample + ' is defined in more than one sample files.')
                        SAMPLEINFO[sample] = info
                        SAMPLE_TO_BATCH[sample] = None #default
                   
                        if len(info.get('readgroups',[])) == 0:
                            no_readgroup.append(sample)
                    if no_readgroup:
                        print('WARNING: %d/%d samples have no readgroups' % (len(no_readgroup), len(w)))
                    

                    #assign to batches for staging from archive or dcache
                    cursize_full = 0 #size if all files need to be staged
                    cursize_actual = 0 #size excluding files that are already retrieved
                    
                    #PROTOCOLS = ['dcache', 'archive']
                    cursize = {p:{'full':0, 'actual':0} for p in PROTOCOLS}

                    
                    new_batch = {p:[] for p in PROTOCOLS}
                    SAMPLEFILE_TO_BATCHES[basename] = {p:[] for p in PROTOCOLS}
                    
                    for sample,info in list(w.items()):
                        #check if sample is on active storage
                        if not info['from_external']: 
                            continue

                        #check if sample is already retrieved
                        info['need_retrieval'] = True
                        filesize = info['filesize']
                        
                        if os.path.exists(os.path.join(os.getcwd(), config['SOURCEDIR'], sample + '.retrieved')) or \
                            os.path.exists(os.path.join(os.getcwd(), config['SOURCEDIR'], sample + '.finished')):
                            info['need_retrieval'] = False
                            filesize=0
                        
                        
                        protocol = info['from_external']
                        
                        #check if batch is full
                        if (cursize[protocol]['full'] + info['filesize']) > MAX_BATCH_SIZE: #stable batch allocation
                            w[sample] = info
                            SAMPLEFILE_TO_BATCHES[basename][protocol].append({'samples':new_batch[protocol], 'size':cursize[protocol]['actual']})
                            new_batch[protocol] = []
                            cursize[protocol]['actual'] = 0
                            cursize[protocol]['full'] = 0
                        #add to batch
                        new_batch[protocol].append(sample)
                        #update batch size
                        cursize[protocol]['full'] += info['filesize']
                        cursize[protocol]['actual'] += filesize
                        

                    
                    
                    for p in PROTOCOLS:
                        #add last batch
                        if new_batch[p]:
                            SAMPLEFILE_TO_BATCHES[basename][p].append({'samples':new_batch[p], 'size':cursize[p]['actual']})

                        #assign batches to samples
                        for pos, batch in enumerate(SAMPLEFILE_TO_BATCHES[basename][p]):
                            for sample in batch['samples']:
                                SAMPLE_TO_BATCH[sample] = f'{p}_{pos}'

        SAMPLE_FILES.sort()


    return (SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES)

def get_validated_sex_file(input):
    #this file should exist after running 'get_validated_sex' job.
    #it should also certainly exist after the bam file is created,
    #as it relies on this file.
    filename = input['validated_sex']
    with open(filename) as f:
        xsample = yaml.load(f,Loader=yaml.FullLoader)
    return 'male' if  xsample['sex'] == 'M' else 'female'



