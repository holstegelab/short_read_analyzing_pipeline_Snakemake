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

def gzipFileSize(filename):
    """return UNCOMPRESSED filesize of a gzipped file.
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
    if fname.endswith('gz'):
        return gzipFileSize(fname)
    elif fname.endswith('bz2'):
        return 1
    else:
        return os.path.getsize(fname)


def get_bam_readgroups(fname, filetype):
    if 'bam' in filetype:
        p = subprocess.Popen('samtools view -H %s | grep ^@RG' % fname, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             shell=True)
    else:
        p = subprocess.Popen('samtools view -H %s -T /home/gozhegov/data/hg38_res/Ref/GRCh38_full_analysis_set_plus_decoy_hla.fa   | grep ^@RG' % fname,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

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




def read_samplefile_simple(filename, config, prefixpath=None):
    orig_filename = filename
    filename = os.path.realpath(filename)

    basename = os.path.splitext(filename)[0] 
    if os.path.exists(basename + '.source'):
        with open(basename + '.source','r') as fsource:
            prefixpath = fsource.readline().strip()
        print(f'Data path overridden from {basename}.source file to {prefixpath}')            

    if not prefixpath:
        prefixpath = os.path.dirname(filename)

    if os.path.exists(basename + '.target'):
        with open(basename + '.target','r') as ftarget:
            targetpath = ftarget.readline().strip()
        print(f'Data target path set to {targetpath}') 
    else:
        targetpath = None


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
                filesize = None
            elif len(row) == 7:
                study, sample_id, file_type, sample_type, capture_kit, sex, filenames1 = row
                filenames2 = ""
                filesize = None
            else:
                study, sample_id, file_type, sample_type, capture_kit, sex, filenames1, filenames2, filesize = row

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

            if filesize is None:
                if sample_type == 'illumina_exome':
                    warning(capture_kit != '', 'Need to fill in capture kit for exome')
                    capture_kits.add(capture_kit)
                    filesize = 8.0 * base_filesize_factor
                else:
                    warning(capture_kit == '' or capture_kit.startswith('WGS'),
                            'Capture kit filled in for ' + sample_type + ' sample')
                    filesize = 75.0 * base_filesize_factor

            filenames1 = [a.strip() for a in filenames1.split(',') if a.strip() != '']
            filenames2 = [a.strip() for a in filenames2.split(',') if a.strip() != '']
            alternative_names.add(
                os.path.splitext(os.path.commonprefix([os.path.basename(e) for e in (filenames1 + filenames2)]))[
                    0].strip('_'))

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
                error(os.path.exists(filenames2[0]),
                      'Fasta reference file for CRAM filetype does not exist: %s' % filenames2[0])
            else:
                warning(len(filenames2) == 0, 'No second filename can be given for bam files: ' + sample_id)

            res = {'samplefile': orig_filename[:-4], 'file1': filenames1, 'file2': filenames2, 'prefix': prefixpath,
                    'target':targetpath,
                   'sample': sample_id, 'filesize': filesize, 'alt_name': alternative_names, 'study': study,
                   'file_type': file_type, 'sample_type': sample_type, 'capture_kit': capture_kit, 'sex': sex}
            
            if any([not append_prefix(prefixpath,f).startswith("/archive") for f in itertools.chain(filenames1,filenames2)]):
                res, warnings = get_readgroups(res, prefixpath)
                for w in warnings:
                    warning(False, w)
            samples.append(res)

    for capture_kit in capture_kits:
        f = os.path.join(config['RES'], config['kit_folder'], capture_kit + '_hg38.bed')
        if not os.path.isfile(f):
            warning(False, 'capture kit file does not exist: ' + f)
    print('')
    print('Check complete')
    return samples


def check(warnings, condition, message):
    if not condition:
        warnings.append(message)


def append_prefix(prefix, filename):
    if not os.path.isabs(filename):
        return os.path.join(prefix, filename)
    else:
        return filename


def get_readgroups(sample, sourcedir):
    sample = sample.copy()
    sample_id = sample['sample']
    sample_type = sample['sample_type']
    study = sample['study']
    filenames1 = [append_prefix(sourcedir, f) for f in sample['file1']]
    filenames2 = [append_prefix(sourcedir, f) for f in sample['file2']]
    alternative_names = sample['alt_name']
    file_type = sample['file_type']
    warnings = []

    readgroups = []
    if file_type == 'fastq_paired':
        for pos, (f1, f2) in enumerate(zip(filenames1, filenames2)):
            info_f1 = read_fastqfile(f1)
            info_f2 = read_fastqfile(f2)

            check(warnings,
                  info_f1['instrument'] == info_f2['instrument'] and info_f1['run'] == info_f2['run'] and info_f1[
                      'flowcell'] == info_f2['flowcell'] and info_f1['lane'] == info_f2['lane'] and \
                  info_f1['index'] == info_f2['index'],
                  'Fastq files for sample ' + sample_id + ' have nonmatching metadata')
            check(warnings, info_f1['pair'] == '1' and info_f2['pair'] == '2',
                  'Fastq files are incorrect order for sample ' + sample_id)
            fs1 = file_size(f1)
            fs2 = file_size(f2)
            check(warnings, fs1 == fs2,
                  'Paired fastq files are unequal in size (%d, %d) for sample ' % (fs1, fs2) + sample_id)

            readgroup_info = {'ID': sample_id + '_rg%d' % pos, \
                              'PL': sample_type, \
                              'PU': info_f1['instrument'] + '.' + info_f1['flowcell'] + '.' + info_f1['lane'] + '.' +
                                    info_f1['index'], \
                              'LB': info_f1['instrument'] + '.' + info_f1['run'] + '.' + sample_id, \
                              'DT': datetime.datetime.fromtimestamp(int(info_f1['date'])).strftime(
                                  '%Y-%m-%dT%H:%M:%S+01:00'), \
                              'CN': study, \
                              'SM': sample_id}
            readgroup = {'info': readgroup_info, 'file_type': file_type, 'file1': f1, 'file2': f2, 'prefix': sourcedir}
            readgroups.append(readgroup)
    elif file_type == 'sra_paired' or file_type == 'sra_single':
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

                readgroup = {'info': readgroup_info, 'file_type': file_type, 'file': sraid,
                             'nreadgroups': len(readgroups_info), 'prefix': sourcedir}
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
            readgroup = {'info': readgroup_info, 'file_type': file_type, 'file': filename, 'nreadgroups': 1}
            readgroups.append(readgroup)
    else:
        used_readgroups = set()
        for filename in [e for e in filenames1 if not e.endswith('crai') or e.endswith('bai')]:
            fstat = os.stat(filename)
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

                readgroup = {'info': readgroup_info, 'file_type': file_type, 'file': filename,
                             'nreadgroups': len(readgroups_info), 'prefix': sourcedir}
                if readgroup_info['ID'] in used_readgroups:
                    readgroup_info['oldname'] = readgroup_info['ID']
                    counter = 2
                    while (readgroup_info['ID'] + '.' + str(counter)) in used_readgroups:
                        counter += 1
                    readgroup_info['ID'] = readgroup_info['ID'] + '.' + str(counter)

                if 'cram' in file_type:
                    assert len(filenames2) == 1
                    readgroup['reference_file'] = filenames2[0]

                used_readgroups.add(readgroup_info['ID'])

                readgroups.append(readgroup)
    check(warnings, len(readgroups) > 0, 'No readgroups defined for sample: ' + sample_id)

    sample['alternative_names'] = alternative_names
    sample['readgroups'] = readgroups
    return (sample, warnings)


def read_fastqfile(filename):
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
            read_fields = read_info.split(':')

            error(len(machine_fields) == 7, 'Unexpected fastq header format: ' + header)
            error(len(read_fields) == 4, 'Unexpected fastq header format: ' + header)

            instrument_name, run_id, flowcell_id, flowcell_lane, tile_number, x, y = machine_fields
            pairid, filtered, control_bits, index_seq = read_fields
        elif '/' in header:  # french format, no multiplex
            machine_info, read_fields = header.split('/')
            machine_fields = machine_info.split(':')
            if len(machine_fields) == 6:
                instrument_name, flowcell_id, flowcell_lane, tile_number, x, y = machine_fields
                run_id = '0'
            elif len(machine_fields) == 7:
                instrument_name, run_id, flowcell_id, flowcell_lane, tile_number, x, y = machine_fields
            else:
                error(len(machine_fields) == 6, 'Unexpected fastq header format: ' + header)

            index_seq = '0'
            pairid = read_fields.strip('\n')
        else:
            error(False, 'Unexpected fastq header formaat: ' + header)
        # BONN has pair id of 3 which means 2 .....???
        if pairid.strip() == '3':
            pairid = '2'

        res = {'file': str(filename), 'instrument': str(instrument_name[1:]), 'run': str(run_id),
               'flowcell': str(flowcell_id), 'lane': str(flowcell_lane), 'index': str(index_seq), 'pair': str(pairid),
               'date': fstat.st_mtime}

    return res

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process sample description files for the preprocessing pipeline.')	
    parser.add_argument('samplefiles', metavar='SAMPLEFILE', nargs='+', help='Path to sample file.')
    parser.add_argument('--config', default='', help='Location of the preprocessing pipeline path config file (default: Snakefile.paths.yaml)')
    parser.add_argument('--data_root', default='', help='Root folder in which to look for source files (referred to from sample file) (default=folder of sample file).')
    args = parser.parse_args()
    data_root = args.data_root

    #load config
    if not args.config:
        configfile = os.path.join(os.path.dirname(os.path.realpath(__file__)),'Snakefile.paths.yaml')
    else:
        configfile = os.path.expanduser(args.config)
    config = utils.read_yaml_config(configfile)

    samplefiles = [os.path.expanduser(sample_file) for sample_file in args.samplefiles]
    for samplefile in samplefiles:
        read_samplefile_simple(samplefile, config, args.data_root)

