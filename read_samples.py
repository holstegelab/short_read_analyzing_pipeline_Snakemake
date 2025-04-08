#!/usr/bin/env python
import os
import os.path
import sys
import gzip
import datetime
import csv
import subprocess
import itertools
import struct
import stat
import bz2
import time
import argparse
import utils

from constants import *
from collections import OrderedDict


PROTOCOLS = ['archive','dcache']

def gzipFileSize(filename):
    """return UNCOMPRESSED filesize of a gzipped file.
    
    :param filename: name of the gzipped file
    :return: uncompressed filesize in bytes
    """
    rsize = os.path.getsize(filename)
    fo = open(filename, 'rb')
    fo.seek(-4, 2)
    r = fo.read()
    fo.close()
    res = struct.unpack('<I', r)[0]

    if res < rsize:
        return 1  # unknown length, gz format only supports max 2**32 (4 gb)
    else:
        return rsize


def file_size(fname):
    """Return the size of a file in bytes.
    
    :param fname: name of the file
    :return: size in bytes.
    """
    if fname.endswith('gz'):
        return gzipFileSize(fname)
    elif fname.endswith('bz2'):
        return 1
    else:
        return os.path.getsize(fname)


def get_bam_readgroups(fname, filetype):
    """Return the readgroups of a bam file.
    :param fname: name of the bam file
    :param filetype: type of the file (bam, cram, extracted_bam, recalibrated_bam, recalibrated_cram, extracted_cram)
    :return: list of readgroups."""
    p = subprocess.Popen('samtools view -H %s | grep ^@RG' % fname, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=True)
       
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(str(fname) + ': ' + str(err))
    if not isinstance(result, str):
        result = result.decode('utf-8')
    rows = result.rstrip('\n').split('\n')

    def rejoin(x):
        if len(x) == 2:
            return tuple(x)
        else:
            return (x[0], ':'.join(x[1:]))

    return [dict([rejoin(e.split(':')) for e in row.strip().split('\t')[1:]]) for row in rows]


def get_sra_readgroups(fname):
    """Return the readgroups of a sra file.
    :param fname: name of the sra file
    :return: list of readgroups."""
    attempts = 3
    success = False
    while not success and attempts > 0:
        p = subprocess.Popen(
            '(sam-dump %s || if [[ $? -eq 141 ]]; then true; else exit $?; fi) | samtools view -H | grep ^@RG' % fname,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, err = p.communicate()
        if p.returncode != 0:
            print('Retry...')
            time.sleep(15)
            attempts -= 1
        else:
            success = True

    if not success:
        print(fname.__class__)
        print(err.__class__)
        print(fname)
        print(err)

        raise IOError(fname + ': ' + err.decode('utf-8'))

    if not isinstance(result, str):
        result = result.decode('utf-8')
    rows = result.rstrip('\n').split('\n')
    res = [[tuple(e.split(':')) for e in row.strip().split('\t')[1:]] for row in rows]
    result = []
    for row in res:
        q = dict()
        for elem in row:
            if len(elem) > 2:
                elem = (elem[0], ':'.join(elem[1:]))
            q[elem[0]] = elem[1]
        result.append(q)
    return result


def warning(test, message):
    if not test:
        print("")
        print('* WARNING: ' + message)
    return not test


def error(test, message):
    if not test:
        print("")
        print('* ERROR: ' + message)
        sys.exit(0)
    return not test


filetypes = set(
    ['fastq_paired', 'bam', 'extracted_bam', 'recalibrated_bam', 'cram', 'recalibrated_cram', 'extracted_cram',
     'sra_paired', 'sra_single', 'gvcf'])
sampletypes = set(['illumina_exome', 'illumina_wgs', 'illumina_wgs_pcr_free', 'illumina_wgs_to_exome'])
sexes = set(['F', 'M'])




def read_samplefile(filename, prefixpath=None):
    """Read a sample file and return a list of samples.

    Sample file format:
    study, sample_id, file_type, sample_type, capture_kit, sex, filenames1[, filenames2[, sample_config]]
    
    if filenames2 is not present, it is assumed to be empty.
    if sample_config is not present, it is assumed to be empty.

    Prefixpath can be given by a .source file in the same directory as the sample file.
    Targetpath can be given by a .target file in the same directory as the sample file.

    Paths can be prepended by a protocol (e.g. archive:/path/to/data) to indicate that the data is not local.
    Accepted protocol values: archive, dcache


    :param filename: name of the sample file
    :param prefixpath: path to prepend to all filenames
    :return: list of sample dictionaries."""
    
    orig_filename = filename
    filename = os.path.realpath(filename)
    #TODO: add check for uniq sample names

    basename = os.path.splitext(filename)[0] 
    if os.path.exists(basename + '.source'):
        with open(basename + '.source','r') as fsource:
            prefixpath = fsource.readline().strip()
        print(f'SOURCE PATH OVERRIDE: by {basename}.source file to {prefixpath}')            

    if not prefixpath:
        prefixpath = os.path.dirname(filename)

    if os.path.exists(basename + '.target'):
        with open(basename + '.target','r') as ftarget:
            targetpath = ftarget.readline().strip()
        print(f'TARGET PATH set to {targetpath}') 
    else:
        targetpath = None
        print(f'TARGET PATH not set (no .target file)') 


    samples = []
    capture_kits = set()
    print('Reading sample file: ' + filename)
    with open(filename, 'r', encoding='utf-8') as f:
        c = csv.reader(f, delimiter='\t')
        c = list(c)
        s = os.stat(os.path.dirname(filename))
        warning(((s.st_mode & stat.S_IWGRP) > 0) & ((s.st_mode & stat.S_IRGRP) > 0),
                'Directory does not have group write/read permission: ' + os.path.dirname(filename))

        for rowpos, row in enumerate(c):
            alternative_names = set()

            error(len(row) == 7 or len(row) == 8 or len(row) == 9, 'Row encountered with != 8 or 9 fields: ' + str(row))
            if len(row) == 8:
                study, sample_id, file_type, sample_type, capture_kit, sex, filenames1, filenames2 = row
                sample_config = {}
            elif len(row) == 7:
                study, sample_id, file_type, sample_type, capture_kit, sex, filenames1 = row
                filenames2 = ""
                sample_config = {}
            else:
                study, sample_id, file_type, sample_type, capture_kit, sex, filenames1, filenames2, sample_config = row
                sample_config = dict([tuple(e.split('=')) for e in sample_config.split(',')])

            warning(sample_id.startswith(study),
                    'Sample id needs to start with study name to prevent sample name conflicts for ' + sample_id)
            sys.stdout.write('- Checking sample ' + sample_id + ' (%d/%d)\r' % (rowpos + 1, len(c)))
            warning(sample_type in sampletypes, 'Unknown sample_type: ' + sample_type)
            warning(sex in sexes, 'Unknown sex: ' + sex)
            error(file_type in filetypes, 'Unknown file_type: ' + file_type)

            if 'cram' in file_type:
                base_filesize_factor = 0.5
            else:
                base_filesize_factor = 1.0

            filesize = sample_config.get('filesize',None)
            if filesize is None:
                if sample_type == 'illumina_exome':
                    warning(capture_kit != '', 'Need to fill in capture kit for exome')
                    capture_kits.add(capture_kit)
                    filesize = 8.0 * base_filesize_factor
                else:
                    warning(capture_kit == '' or capture_kit.startswith('WGS'),
                            'Capture kit is not empty or WGS for a' + sample_type + ' sample')
                    filesize = 75.0 * base_filesize_factor

            filenames1 = [a.strip() for a in filenames1.split(',') if a.strip() != '']
            filenames2 = [a.strip() for a in filenames2.split(',') if a.strip() != '']
            cpref = os.path.splitext(os.path.commonprefix([os.path.basename(e) for e in (filenames1 + filenames2)]))[0].strip('_')
            if len(cpref) > 8:
                alternative_names.add(cpref)

            cram_refs = []
            if file_type == 'fastq_paired':
                warning(len(filenames1) > 0, 'No filename given for sample ' + sample_id)
                if len(filenames2) == 0:
                    file_type = 'fastq_interleaved'
                    error(False,
                          'Handling interleaved fastq files is not yet implemented, let me know if you need this')
                else:
                    warning(len(filenames1) == len(filenames2),
                            'Number of fastq files is not equal for filenames_read1 and filenames_read2 for sample ' + sample_id)
            elif file_type == 'sra_paired' or file_type == 'sra_single':
                warning(len(filenames2) == 0, 'No second filename can be given for SRA files: ' + sample_id)
            elif file_type == 'gvcf':
                warning(len(filenames2) == 0, 'No second filename can be given for GVCF files: ' + sample_id)
                warning(len(filenames1) == 1, 'Only a single GVCF file can be given: ' + sample_id)
            elif file_type == 'cram' or file_type == 'recalibrated_cram' or file_type == 'extracted_cram':
                error(len(filenames2) == 1 and (filenames2[0].endswith('fa') or filenames2[0].endswith('fasta')),
                      'Second filename for CRAM filetype should be fasta reference file')
                cram_refs = filenames2
                filenames2 = []
            else:
                warning(len(filenames2) == 0, 'No second filename can be given for bam files: ' + sample_id)

            res = {'samplefile': orig_filename[:-4], 'file1': filenames1, 'file2': filenames2, 'prefix': prefixpath,
                    'target':targetpath,
                   'sample': sample_id, 'filesize': filesize, 'alt_name': alternative_names, 'study': study,
                   'file_type': file_type, 'sample_type': sample_type, 'capture_kit': capture_kit, 'sex': sex, 'no_dedup':sample_config.get('no_dedup',False), 'cram_refs':cram_refs}
                
            all_files = [append_prefix(prefixpath,f) for f in itertools.chain(filenames1,filenames2)] 
            protocols = set([f.split(':')[0] for f in all_files if ':' in f])

            if protocols:
               assert len(protocols) == 1, f'Cannot combine multiple external data sources for sample {sample_id}: {protocols}'
               assert all([e in PROTOCOLS for e in protocols]), f'Unknown protocol for sample {sample_id}: {protocols}'

               res['from_external'] = protocols.pop()
            else:
                res, warnings = get_readgroups(res, prefixpath)
                for w in warnings:
                    warning(False, w)
                res['from_external'] = False

            samples.append(res)

    for capture_kit in capture_kits:
        f = pj(INTERVALS_DIR, capture_kit + '.bed')
        if not os.path.isfile(f):
            warning(False, 'capture kit file does not exist: ' + f)
    print("")
    print('Check complete')
    return samples


def check(warnings, condition, message):
    """Checks a condition and appends a warning if it is not met"""
    if not condition:
        warnings.append(message)


def append_prefix(prefix, filename):
    """Appends a prefix to a filename if it is not an absolute path"""
    if not os.path.isabs(filename):
        return os.path.join(prefix, filename)
    else:
        return filename


def get_readgroups(sample, sourcedir):
    """Obtains readgroups from fastq, SRA, or bam/cram files. Stores in the sample dictionary.
    :param sample: sample dictionary
    :param sourcedir: directory where the files are located
    
    :return: sample dictionary, list of warnings
    """
    # Make a copy of the sample dictionary to avoid modifying the original
    sample = sample.copy()
    
    # Extract relevant information from the sample dictionary
    sample_id = sample['sample']
    sample_type = sample['sample_type']
    study = sample['study']
    filenames1 = [append_prefix(sourcedir, f) for f in sample['file1']]
    filenames2 = [append_prefix(sourcedir, f) for f in sample['file2']]
    alternative_names = sample['alt_name']
    file_type = sample['file_type']
    warnings = []
    
    # Initialize an empty list to store the readgroups
    readgroups = []
    
    # Check the file type and obtain readgroups accordingly
    if file_type == 'fastq_paired':
        # Obtain readgroups from paired-end fastq files
        for pos, (f1, f2) in enumerate(zip(filenames1, filenames2)):
            info_f1 = read_fastqfile(f1)
            info_f2 = read_fastqfile(f2)
            
            # Check that the metadata of the fastq files match
            check(warnings,
                  info_f1['instrument'] == info_f2['instrument'] and info_f1['run'] == info_f2['run'] and info_f1[
                      'flowcell'] == info_f2['flowcell'] and info_f1['lane'] == info_f2['lane'] and \
                  info_f1['index'] == info_f2['index'],
                  'Fastq files for sample ' + sample_id + ' have nonmatching metadata')
            
             # Check that the fastq files are in the correct order
            check(warnings, info_f1['pair'] == '1' and info_f2['pair'] == '2',
                  'Fastq files are incorrect order for sample ' + sample_id)
            
             # Check that the paired fastq files have equal size
            fs1 = file_size(f1)
            fs2 = file_size(f2)
            check(warnings, fs1 == fs2,
                  'Paired fastq files are unequal in size (%d, %d) for sample ' % (fs1, fs2) + sample_id)

            # Define the readgroup information
            readgroup_info = {'ID': sample_id + '_rg%d' % pos, \
                              'PL': sample_type, \
                              'PU': info_f1['instrument'] + '.' + info_f1['flowcell'] + '.' + info_f1['lane'] + '.' +
                                    info_f1['index'], \
                              'LB': info_f1['instrument'] + '.' + info_f1['run'] + '.' + sample_id, \
                              'DT': datetime.datetime.fromtimestamp(int(info_f1['date'])).strftime(
                                  '%Y-%m-%dT%H:%M:%S+01:00'), \
                              'CN': study, \
                              'SM': sample_id}
            # Define the readgroup dictionary
            readgroup = {'info': readgroup_info, 'file_type': file_type, 'file1': f1, 'file2': f2, 'prefix': sourcedir}
            
            # Append the readgroup to the list of readgroups
            readgroups.append(readgroup)
    elif file_type == 'sra_paired' or file_type == 'sra_single':
        # Obtain readgroups from SRA file
        for sraid in filenames1:
            try:
                readgroups_info = get_sra_readgroups(sraid)
            except IOError as e:
                check(warnings, False, str(e))
                raise
            for readgroup_info in readgroups_info:
                if not 'DT' in readgroup_info:
                    readgroup_info['DT'] = ""
                if not 'CN' in readgroup_info:
                    readgroup_info['CN'] = study
                if not 'PL' in readgroup_info:
                    readgroup_info['PL'] = sample_type
                if 'SM' in readgroup_info:
                    alternative_names.add(readgroup_info['SM'])
                readgroup_info['SM'] = sample_id
                
                # Define the readgroup dictionary
                readgroup = {'info': readgroup_info, 'file_type': file_type, 'file': sraid,
                             'nreadgroups': len(readgroups_info), 'prefix': sourcedir}
                # Append the readgroup to the list of readgroups
                readgroups.append(readgroup)
    elif file_type == 'gvcf':
        for filename in filenames1:
            fstat = os.stat(filename)
            readgroup_info = {}
            filedate = fstat.st_mtime
            readgroup_info['ID'] = sample_id + '_1'
            readgroup_info['DT'] = datetime.datetime.fromtimestamp(filedate).strftime('%Y-%m-%dT%H:%M:%S+01:00')
            readgroup_info['CN'] = study
            readgroup_info['PL'] = sample_type
            readgroup_info['SM'] = sample_id
            
            # Define the readgroup dictionary
            readgroup = {'info': readgroup_info, 'file_type': file_type, 'file': filename, 'nreadgroups': 1}
            
            # Append the readgroup to the list of readgroups
            readgroups.append(readgroup)
    else:
        # Obtain readgroups from bam/cram files
        used_readgroups = set()
        for filename in [e for e in filenames1 if not e.endswith('crai') or e.endswith('bai')]:
            fstat = os.stat(filename)
            if 'cram' in file_type: #need fasta reference file
                if len(sample['cram_refs']) > 0:
                    assert len(sample['cram_refs']) == 1, 'No support (yet) for multiple cram reference files in "file2" column'
                    reference_file = sample['cram_refs'][0]
                else:                    
                    reference_file = "GRCh38_full_analysis_set_plus_decoy_hla.fa"  #default reference file

                if not reference_file.startswith('/'): #relative path, look in cram reference folder
                    reference_file = pj(CRAMREFS, sample['cram_refs'][0])
            else:
                reference_file = None

            readgroups_info = get_bam_readgroups(filename, file_type)

            for readgroup_info in readgroups_info:
                if not 'DT' in readgroup_info:
                    filedate = fstat.st_mtime
                    readgroup_info['DT'] = datetime.datetime.fromtimestamp(filedate).strftime('%Y-%m-%dT%H:%M:%S+01:00')
                if not 'CN' in readgroup_info:
                    readgroup_info['CN'] = study
                if not 'PL' in readgroup_info:
                    readgroup_info['PL'] = sample_type
                if 'SM' in readgroup_info:
                    alternative_names.add(readgroup_info['SM'])
                readgroup_info['SM'] = sample_id
                
                # Define the readgroup dictionary
                readgroup = {'info': readgroup_info, 'file_type': file_type, 'file': filename,
                             'nreadgroups': len(readgroups_info), 'prefix': sourcedir, 'reference_file': reference_file}
                
                # Check if the readgroup ID has already been used
                if readgroup_info['ID'] in used_readgroups:
                    readgroup_info['oldname'] = readgroup_info['ID']
                    counter = 2
                    while (readgroup_info['ID'] + '.' + str(counter)) in used_readgroups:
                        counter += 1
                    readgroup_info['ID'] = readgroup_info['ID'] + '.' + str(counter)

                used_readgroups.add(readgroup_info['ID'])
                
                # Append the readgroup to the list of readgroups
                readgroups.append(readgroup)
    
    # Check that at least one readgroup was defined
    check(warnings, len(readgroups) > 0, 'No readgroups defined for sample: ' + sample_id)

    # Update the sample dictionary with the alternative names and readgroups
    sample['alternative_names'] = alternative_names
    sample['readgroups'] = readgroups
    # Return the updated sample dictionary and the list of warnings
    return (sample, warnings)


def read_fastqfile(filename):
    """
    Reads a fastq file and returns a dictionary containing the information from the header line.
    :param filename: the fastq file
    :return: a dictionary containing the information from the header line

    """
    if filename.endswith('gz'):
        o = gzip.open
    elif filename.endswith('bz2'):
        o = bz2.BZ2File
    else:
        o = open
        warning(False, 'Fastq file ' + filename + ' uncompressed.')

    fstat = os.stat(filename)
    with o(filename, 'r') as f:
        header = f.readline()
        if not isinstance(header, str):
            header = header.decode('utf-8')
        

        header = header.strip()
        if '#' in header:
            machine_info, read_info = header.split('#')
            machine_fields = machine_info.split(':')
            read_fields = read_info.split('/')
            error(len(machine_fields) == 5, 'Unexpected fastq header format: ' + header)
            error(len(read_fields) == 2, 'Unexpected fastq header format: ' + header)
            instrument_name, flowcell_lane, tile_number, x, y = machine_fields
            index_seq, pairid = read_fields
            run_id = '0'
            flowcell_id = '0'

        elif ' ' in header:  # casava 1.8
            header = header.strip()
            error(' ' in header, 'Unexpected fastq header format: ' + header)
            machine_info, read_info = header.split(' ')
            machine_fields = machine_info.split(':')
            
            if len(machine_fields) == 8:
                instrument_name, run_id, flowcell_id, flowcell_lane, tile_number, x, y, dummy = machine_fields

            elif len(machine_fields) == 7:
                instrument_name, run_id, flowcell_id, flowcell_lane, tile_number, x, y = machine_fields
            elif len(machine_fields) == 6:
                instrument_name, flowcell_id, flowcell_lane, tile_number, x, y = machine_fields
            elif len(machine_fields) == 1:
                instrument_name = 'UNKNOWN'
                run_id = machine_fields[0]
                flowcell_id = '0'
                flowcell_lane = '0'
                tile_number = '0'
                x = '0'
                y = '0'                   
            else:     
                error(False, 'Unexpected fastq header format: ' + header)

            if '/' in read_info:   
                read_info, pairid = read_info.split('/')
                          
            read_fields = read_info.split(':')             
            
            if len(read_fields) == 1:
                index_seq = read_fields[0]

            elif len(read_fields) == 4:
                pairid, filtered, control_bits, index_seq = read_fields                
            elif len(read_fields) == 5:
                flowcell_id, flowcell_lane, tile_number, x,y = read_fields
                index_seq = ''
            elif len(read_fields) == 7:
                instrument_name, run_id, flowcell_id, flowcell_lane, tile_number, x, y = read_fields
                index_seq = ''
            else:
                error(False, 'Unexpected fastq header format: ' + header)
            
            
        elif '/' in header:  # french format, no multiplex
            machine_info, read_fields = header.split('/')
            machine_fields = machine_info.split(':')
            if len(machine_fields) == 6:
                instrument_name, flowcell_id, flowcell_lane, tile_number, x, y = machine_fields
                run_id = '0'
            elif len(machine_fields) == 7:
                instrument_name, run_id, flowcell_id, flowcell_lane, tile_number, x, y = machine_fields
            elif len(machine_fields) == 1:
                instrument_name = 'Artificial_reads'
                run_id = '0'
                flowcell_id = '0'
                flowcell_lane = '0'
                tile_number = '0'
                x = '0'
                y = '0'
            else:
                error(False, 'Unexpected fastq header format: ' + header)

            index_seq = '0'
            pairid = read_fields.strip('\n')
        else:
            pairid = '0'
            error(False, 'Unexpected fastq header formaat: ' + header)
        # BONN has pair id of 3 which means 2 .....???
        if pairid.strip() == '3':
            pairid = '2'

        res = {'file': str(filename), 'instrument': str(instrument_name[1:]), 'run': str(run_id),
               'flowcell': str(flowcell_id), 'lane': str(flowcell_lane), 'index': str(index_seq), 'pair': str(pairid),
               'date': fstat.st_mtime}

    return res



############# CODE used in Snakemake pipeline ######################
SAMPLEFILE_TO_SAMPLES = {}
MAX_BATCH_SIZE = 2 * 1024 #2TB


#function to read in and cache a samplefile 
def samplefile(sfilename):
    basename = os.path.splitext(os.path.basename(sfilename))[0]
    if not basename in SAMPLEFILE_TO_SAMPLES: 
        if not os.path.isfile(sfilename):
            raise RuntimeError('Sample file ' + sfilename + ' does not exist')
       
        datfilename = os.path.realpath(sfilename)[:-4] + '.adat'
        if not os.path.isfile(datfilename):
            if not os.path.isfile(sfilename):
                raise RuntimeError('Sample file ' + sfilename + ' does not exist')
            samplelist = read_samplefile(sfilename)
            #using ordereddict to main stable order
            sampleinfodict = OrderedDict([(a['sample'], a) for a in samplelist])
            utils.save(sampleinfodict,datfilename)
        else:
            sampleinfodict = utils.load(datfilename)

        SAMPLEFILE_TO_SAMPLES[basename] = sampleinfodict
    return SAMPLEFILE_TO_SAMPLES[basename]


def load_samplefiles(filedir, cache):
    
    if not 'SAMPLE_FILES' in cache:
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
                    w = samplefile(f)
                    
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

                        archive_retrieved = os.path.exists(os.path.join(os.getcwd(), SOURCEDIR, sample + '.archive_retrieved'))
                        dcache_retrieved = os.path.exists(os.path.join(os.getcwd(), SOURCEDIR, sample + '.dcache_retrieved'))

                        if (archive_retrieved or dcache_retrieved) or \
                            os.path.exists(os.path.join(os.getcwd(), SOURCEDIR, sample + '.finished')):
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

###############################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process sample description files for the preprocessing pipeline.')	
    parser.add_argument('samplefiles', metavar='SAMPLEFILE', nargs='+', help='Path to sample file.')
    
    parser.add_argument('--data_root', default='', help='Root folder in which to look for source files (referred to from sample file) (default=folder of sample file).')
    args = parser.parse_args()
    data_root = args.data_root


    samplefiles = [os.path.expanduser(sample_file) for sample_file in args.samplefiles]
    for samplefile in samplefiles:
        read_samplefile(samplefile, args.data_root)
