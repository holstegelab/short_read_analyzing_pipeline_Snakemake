import pandas as pd
import read_stats
import os
import getpass

configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")

ref = os.path.join(config['RES'],config['ref'])
tmpdir = os.path.join(config['TMPDIR'],getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

wildcard_constraints:
    sample="[\w\d_\-@]+",
    extension='sam|bam|cram',
    filetype = 'fq|fastq',
    sex = 'male|female',
    batchnr='[\d]+'
    # readgroup="[\w\d_\-@]+"
module Reference_preparation:
    snakefile: "Reference_preparation.smk"
    config: config
from read_samples import *
import utils
from common import *

SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)
#ARCHIVE/DCACHE handling: it is not efficient to get files from tape file by file.
#For each file we would have to wait for the tape robot to get the tape and spin to the right position.
#Better is to first stage a batch of files together (preferably from the same tape). Once they are available,
# we can immediately copy them to active storage. 
#
#Therefore 'load_samplefiles' defines batches, which are stages/copied together. 
# this allows the tape robot to stage multiple files at once.
# and then to copy them to active storage

# extract all sample names from SAMPLEINFO dict to use it rule all


sample_names = SAMPLEINFO.keys()
sample_sex_names = get_sample_sex_combi(SAMPLEINFO)


def get_refdir_by_sex(wildcards):
    if wildcards['sex'] == 'female':
        ref_dir=os.path.join(config['RES'],config['ref_female_dir'])
    else:
        ref_dir=os.path.join(config['RES'],config['ref_male_dir'])

    return ref_dir

def sampleinfo(SAMPLEINFO, sample, checkpoint=False):
    """If samples are on tape, we do not have sample readgroup info.
    That is, the 'readgroups' field is empty.

    This function first checks if the readgroup info is available on disk,
    in the file config['SAMPLEINFODIR']/<sample>.dat. 

    Alternatively, the function injects a checkpoint rule to load this readgroup info.
    """

    sinfo = SAMPLEINFO[sample]
    if not 'readgroups' in sinfo:
        rgpath = os.path.join(config['SAMPLEINFODIR'], sample + ".dat")
        if os.path.exists(rgpath):
            xsample = utils.load(rgpath)
        elif checkpoint: 
            #no readgroup info yet
            filename = checkpoints.get_readgroups.get(sample=sample).output[0]
            xsample = utils.load(filename)
        sinfo = sinfo.copy()
        sinfo['readgroups'] = xsample['readgroups']
        sinfo['alternative_names'] = sinfo.get('alternative_names',set()).union(xsample['alternative_names'])
        SAMPLEINFO[sample] = sinfo
    return sinfo

def ensure_readgroup_info(wildcards):
    """Make sure the readgroup info for a sample is available."""
    sample = wildcards['sample']
    sinfo = SAMPLEINFO[sample]
    if not 'readgroups' in sinfo:       
        filename = checkpoints.get_readgroups.get(sample=sample).output[0]
        return [filename]
    return []


rule Aligner_all:
    input:
        expand("{cram}/{sample_sex}.mapped_hg38.cram",sample_sex=sample_sex_names, cram = config['CRAM']),
        rules.Reference_preparation_all.input
    default_target: True


def get_source_files(wildcards):
    """Make sure the source files for a sample are available.

    When they need to be obtained from archive or dcache (as indicated by the prefix 'archive:' or 'dcache:'), 
    request an indicator file which is written by the archive_to_active/dcache_to_active rules when all sample files are retrieved.
    """

    sinfo = SAMPLEINFO[wildcards['sample']]
    prefixpath = sinfo['prefix']
    files = []
    archive_retrieve_add = False
    dcache_retrieve_add = False
    for f in itertools.chain(sinfo['file1'],sinfo['file2']):
        if not f:
            continue
        f = append_prefix(prefixpath, f)

        if f.startswith('archive:'):
            if not archive_retrieve_add:
                files.append(ancient(os.path.join(config['SOURCEDIR'], wildcards['sample'] + '.archive_retrieved')))
                archive_retrieve_add = True
        elif f.startswith('dcache:'):
            if not archive_retrieve_add:
                files.append(ancient(os.path.join(config['SOURCEDIR'], wildcards['sample'] + '.dcache_retrieved')))
                archive_retrieve_add = True
        else:
            files.append(f)

    return files

checkpoint get_readgroups:
    """Get the readgroup info for a sample.
    
    Once the checkpoint rule is executed, the readgroup info is available in the file config['SAMPLEINFODIR']/<sample>.dat.

    Once readgroup info is available, snakemake will recalculate the DAG.
    """
    input:
        get_source_files
    output:
        os.path.join(config['SAMPLEINFODIR'], "{sample}.dat")
    run:
        sample = SAMPLEINFO[wildcards['sample']]
        if sample['from_external']:
            prefixpath = config['SOURCEDIR'] # location where the data is downloaded from the external data repository
        else:
            prefixpath = sample['prefix']               
        sample,warnings = read_samples.get_readgroups(sample, prefixpath, config)
        if warnings:
            warningfile = os.path.join(config['SAMPLEINFODIR'], wildcards['sample'] + '.warnings')
            if warnings:
                with open(warningfile,'w') as f:
                    for w in warnings:
                        print ("WARNING: " + w)
                        f.write(w + '\n')
        utils.save(sample, str(output))

rule archive_get:
    """Stage a batch of files from archive.

    The batch is defined in SAMPLEFILE_TO_BATCHES.
    Once the batch is retrieved, an indicator file is written to indicate that the batch is available.
    """

    output:
        os.path.join(config['FETCHDIR'], '{samplefile}.archive_{batchnr}.retrieved')
    resources:
        arch_use_add=lambda wildcards: SAMPLEFILE_TO_BATCHES[wildcards['samplefile']]['archive'][int(wildcards['batchnr'])]['size'],
        partition="archive",
        n=1,
        mem_mb=100
    run:
        dname = os.path.dirname(str(output))
        batch = SAMPLEFILE_TO_BATCHES[wildcards['samplefile']]['archive'][int(wildcards['batchnr'])]
        files = []
        for sample in batch['samples']:
            sinfo = SAMPLEINFO[sample]
            if not sinfo['need_retrieval']:
                continue

            files1 = sinfo['file1']
            files2 = sinfo['file2']
            prefixpath = sinfo['prefix']
            
            files.extend([append_prefix(prefixpath, e) for e in itertools.chain(files1,files2) if e and not os.path.isabs(e)])
            
        archive_files = " ".join([f.replace('archive:/','') for f in files if f.startswith('archive:/')])
        shell("dmget -a " + archive_files)
        shell("touch {output}")

def retrieve_batch(wildcards):
    """For a sample, determines which batch needs to be retrieved from tape."""

    sample = SAMPLEINFO[wildcards['sample']]
    batch = SAMPLE_TO_BATCH[wildcards['sample']]

    if sample['from_external']:
        #batch has format <protocol>_<batchnr>  (e.g. 'archive_0')
        if batch is None:
            # sample is not assigned to a batch, this should generate an error
            return ancient(os.path.join(config['FETCHDIR'], wildcards['sample'] + ".finished_samples_not_assigned_to_retrieval_batch"))
        else:
            return ancient(os.path.join(config['FETCHDIR'], os.path.basename(sample['samplefile'])+ f".{batch}.retrieved"))

    else:
        return []

def calculate_active_use(wildcards):
    """Calculate the amount of active storage that will be used by this sample."""

    sample = SAMPLEINFO[wildcards['sample']]
    filesize = sample['filesize']
    capture_kit = sample['capture_kit']
    active_filesize = 2.0 * filesize if 'cram' in sample['file_type'] else filesize
    res = 0
    if sample['from_external']:
        res += filesize
    if capture_kit == 'WGS38_to_exome' or capture_kit == 'WGS37_to_exome':
        res += 2.0 * active_filesize * 0.15
    else:
        res += 2.0 * active_filesize
    return res


rule start_sample:
    """Start processing a sample. 
    
    This rule can only start if the batch to which the sample belongs has been retrieved from tape.
    This rule will reserve the space on active storage.
    
    """
    input:
        retrieve_batch
    output:
        os.path.join(config['SOURCEDIR'], "{sample}.started")
    resources:
        active_use_add=calculate_active_use,
        mem_mb=50,
        n=1
    shell: """
        touch {output}
        """

rule archive_to_active:
    """Move the files of a sample from archive to active storage.

    Can only start after 'start_sample' has been executed. 
    Start_sample will reserve the space on active storage.
    
    """
    input:
        retrieve_batch,
        ancient(os.path.join(config['SOURCEDIR'], "{sample}.started"))
    resources:
        arch_use_remove=lambda wildcards: SAMPLEINFO[wildcards['sample']]['filesize'],
        partition="archive",
        n=1,
        mem_mb=50
    output:
        temp(os.path.join(config['SOURCEDIR'], "{sample}.archive_retrieved"))
    run:
        sample = SAMPLEINFO[wildcards['sample']]
        prefixpath = sample['prefix']
        destinationpath = os.path.dirname(str(output))
        files1 = sample['file1']
        files2 = sample['file2']
        for e in itertools.chain(files1,files2):
            if not e or os.path.isabs(e):
                continue
            source = os.path.join(prefixpath, e)
            if ':/' in source: #remove protocol
                source = source.split(':/')[1]

            destination = os.path.join(destinationpath, e)
            shell("""
            mkdir -p `dirname {destination}`
            rsync --size-only {source} {destination}
            dmput -r -w {source}
            """)
                
        shell("""
            touch {output}
        """)
        


def get_cram_ref(wildcards):
    """utility function to get cram reference file option for samtools
    the read_samples utility checks that this filename is set if the file type is cram

    :return: string with cram reference file option for samtools, e.g. "--reference /path/to/ref.fa" or "" if file type is not cram
    """
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)

    wildcards = dict(wildcards)
    if 'readgroup' in wildcards:
        readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    else:       
        readgroup = [readgroup for readgroup in sinfo['readgroups'] if wildcards['filename'] in readgroup['file']][0]
        
    if readgroup['file_type'] == 'cram':
        cram_options = '--reference ' + readgroup['reference_file']
    else:
        cram_options = ''
    return cram_options


def get_source_aligned_file(wildcards): 
    """utility function to get the path to the source file for a sample.

    Takes the wildcard filename and looks up the path to that filename in the sampleinfo dictionary.
    If the sample is from external, the path is relative to the SOURCEDIR,
    otherwise it is equal to the prefixpath, which is the path of the readgroup file or the path set in the .source file.

    :return: string with path to source file
    """
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    #there might be multiple bam/cram files as source for a sample (e.g. if sequenced multiple times)
    #look for a read group for which the source file matches wildcard 'filename'
    #return the full file path

    readgroup = [readgroup for readgroup in sinfo['readgroups'] if wildcards['filename'] in readgroup['file']][0]

    return readgroup['file']


rule split_alignments_by_readgroup:
    """Split a sample bam/cram file into multiple readgroups.
    Readgroup bam/cram files are stored in the READGROUPS/<sample>_<sourcefilename> folder.

    The filenames in this folder are equal <sample>.<readgroup_id>.<extension>
    """
    input:
        get_source_aligned_file,
        ancient(os.path.join(config['SOURCEDIR'], "{sample}.started"))
    output:
        #there can be multiple read groups in 'filename'. Store them in this folder.        
        readgroups=temp(directory(os.path.join(config['READGROUPS'], "{sample}_{filename}")))
    resources:
        n=8,
        mem_mb=1000
    params:
        cramref=get_cram_ref
    run:
        sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)        
        readgroups = [readgroup for readgroup in sinfo['readgroups'] if wildcards['filename'] in readgroup['file']]

        n = len(readgroups)
        if n == 1: #nothing to split, single readgroup, just link the source file
            extension = os.path.splitext(input[0])[1][1:].lower()
            readgroup_id = readgroups[0]['info']['ID']
            
            cmd = """
                mkdir -p {output}
                #switching to cp instead of hard link as hard links als update modification time of input[0]
                cp {input[0]} {output}/{wildcards.sample}.{readgroup_id}.{extension}
                """
            shell(cmd)                    
        else:
            #extract readgroups to output folder with {sample}.{readgroups ID}.{extension} format
            if readgroups[0]['file_type'] == 'cram': 
                #keep cram format, much more efficient, and even slighly faster
                output_fmt = 'cram,version=3.1'
                extension = 'cram'
            else:
                output_fmt = 'bam'
                extension = 'bam'

            cmd = """
                mkdir -p {output}
                samtools split -@ {resources.n} --output-fmt {output_fmt} {params.cramref} {input[0]} -f "{output}/{wildcards.sample}.%!.{extension}"

                """
            shell(cmd,conda_env=config['CONDA_MAIN'])



def get_aligned_readgroup_folder(wildcards):
    """Utility function to get the path to the folder containing the readgroup files for a sample.
    
    This is the READGROUPS/<sample>_<sourcefilename> folder.
    """
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    sfile = os.path.splitext(os.path.basename(readgroup['file']))[0]
    folder = os.path.join(config['READGROUPS'], wildcards['sample'] + '_' + sfile)
    return folder

def get_extension(wildcards):
    """Utility function to get the extension of the input file for a sample (bam/cram)."""
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    res =  os.path.splitext(readgroup['file'])[1][1:].lower()
    
    return res



rule external_alignments_to_fastq:
    """Convert a sample bam/cram file to fastq files.
    """
    input:
        get_aligned_readgroup_folder
    output:
        fq1 = temp(config['FQ'] + "/{sample}.{readgroup}_R1.fastq.gz"),
        fq2 = temp(config['FQ'] + "/{sample}.{readgroup}_R2.fastq.gz"),
        singletons = temp(config['FQ'] + "/{sample}.{readgroup}.extracted_singletons.fq.gz"),
    resources:
        n= 2,
        mem_mb = 5000,
        tmpdir=tmpdir
    log:
        fastq=os.path.join(config['LOG'],"{sample}.{readgroup}.align2fq.log"),
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.align2fq.txt")
    params:
        cramref=get_cram_ref,
        extension=get_extension
    priority: 10
    conda: config['CONDA_MAIN']
    #alternative: collate can run also in fast mode (e.g. -r 100000 -f), but this has potential impact on alignment (estimation of insert size becomes biased to genome location)
    #samtools fastq only uses ~2 threads, while collate uses much more. Therefore,pipe directly into adapterremoval to make use of the extra cores. Also save a store/load cycle.
    shell:
        """
            samtools collate -u -@ {resources.n} {params.cramref} -O {input}/{wildcards.sample}.{wildcards.readgroup}.{params.extension} {resources.tmpdir}/{wildcards.sample}.{wildcards.readgroup}.collate |
            samtools fastq -O -N -@ {resources.n} -0 /dev/null -1 {output.fq1} -2 {output.fq2} -s {output.singletons}  2> {log.fastq} 
        """


rule fastq_bz2togz:
    """Convert a bz2 compressed fastq file to a gz compressed fastq file.
    filetype can be 'fq' or 'fastq'
    """
    input:
        "{path}.{filetype}.bz2"
    output:
        temp("{path}.{filetype}.gz")
    shell: """
        bzcat {input} | bgzip > {output}
        """



def get_fastqpaired(wildcards):
    """Utility function to get the path to the fastq files for a sample."""
    
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    # check readgroups
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    
    if sinfo['file_type'] == 'fastq_paired' or sinfo['file_type'] == 'fastq':
        file1 = readgroup['file1']
        if file1.endswith('.bz2'):
            file1 = file1[:-4] + '.gz'
        file2 = readgroup['file2']
        if file2.endswith('.bz2'):
            file2 = file2[:-4] + '.gz'
    else:    
        file1 = config['FQ'] + f"/{wildcards['sample']}.{wildcards['readgroup']}_R1.fastq.gz"
        file2 = config['FQ'] + f"/{wildcards['sample']}.{wildcards['readgroup']}_R2.fastq.gz"
    return [file1, file2]


rule adapter_removal:
    """Remove adapters from fastq files."""
    input:
        get_fastqpaired,
        ancient(os.path.join(config['SOURCEDIR'], "{sample}.started"))
    output:
        for_f=temp(os.path.join(config['FQ'],"{sample}.{readgroup}.fastq.cut_1.fq.gz")),
        rev_f=temp(os.path.join(config['FQ'],"{sample}.{readgroup}.fastq.cut_2.fq.gz"))
    # log file in this case contain some stats about removed seqs
    log:
        adapter_removal=os.path.join(config['STAT'],"{sample}.{readgroup}.adapter_removal.log"),
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.adapter_removal.txt")
    priority: 10
    conda: "envs/preprocess.yaml"
    resources: 
        n=4,
        mem_mb=200
    ##FIXME: combine adapter removal and adapter identify, use paste <(pigz -cd  test_r1cut.f1.gz | paste - - - -) <(pigz -cd test_r2cut.fq.gz | paste - - - -) |  tr '\t' '\n' |
    shell:
        """
		    AdapterRemoval --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --file1 {input[0]} --file2 {input[1]} --gzip --gzip-level 1 --output1 {output.for_f} --output2 {output.rev_f} --settings {log.adapter_removal} --minlength 40 --singleton /dev/null --discarded /dev/null --threads {resources.n} 
		"""

rule adapter_removal_identify:
    """Identify adapters in fastq files."""
    input:
        get_fastqpaired,
        ancient(os.path.join(config['SOURCEDIR'], "{sample}.started"))
    output:
        stats=os.path.join(config['STAT'],"{sample}.{readgroup}.fastq.adapters"),
    priority: 5
    conda: "envs/preprocess.yaml"
    resources: 
        n=4,
        mem_mb=200
    benchmark: os.path.join(config['BENCH'],"{sample}.{readgroup}.adapter_removal_identify.txt")
    shell:
        """
		AdapterRemoval --identify-adapters --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --file1 {input[0]} --file2 {input[1]}  --threads {resources.n} > {output.stats}
		"""


def get_readgroup_params(wildcards):
    """Utility function to get the readgroup params for a sample.
       Fills in missing values with 'unknown' to avoid errors in downstream tools.
    """
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    res = [rg for rg in sinfo['readgroups'] if rg['info']['ID'] == wildcards['readgroup']][0][
        'info']

    return {'ID': res['ID'], 'LB': res.get('LB','unknown'), 'PL': res.get('PL','unknown'),
            'PU': res.get('PU','unknown'), \
            'CN': res.get('CN','unknown'), 'DT': res.get('DT','unknown')}


# rule to align reads from cutted fq on hg38 ref
# use dragmap aligner
# samtools fixmate for future step with samtools mark duplicates

def get_mem_mb_align_reads(wildcards, attempt):
    """Utility function to get the memory for the align_reads rule.

    """
    MEM_DEFAULT_USAGE = 44000
    return (attempt - 1) * 0.5 * int(MEM_DEFAULT_USAGE) + int(MEM_DEFAULT_USAGE)

def get_prepared_fastq(wildcards):
    """Utility function to get the path to the (adapter-removed) fastq files for a sample."""
    file1 = os.path.join(config['FQ'], wildcards['sample'] + '.' + wildcards['readgroup'] + '.fastq.cut_1.fq.gz')
    file2 = os.path.join(config['FQ'], wildcards['sample'] + '.' + wildcards['readgroup'] + '.fastq.cut_2.fq.gz')
    return [file1,file2]


rule align_reads:
    """Align reads to reference genome."""
    input:
        get_prepared_fastq,
        ensure_readgroup_info
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.{sex}.aligned.bam")
    params:
        ref_dir=get_refdir_by_sex,
        dragmap=os.path.join(config['RES'], config['SOFTWARE'],'dragen-os'),
        rg_params=get_readgroup_params
    conda: "envs/preprocess.yaml"
    log:
        dragmap_log=os.path.join(config['LOG'],"{sample}.{readgroup}.{sex}.dragmap.log"),
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.{sex}.dragmap.txt")
    priority: 15
    resources:
        n=24, #reducing thread count, as first part of dragmap is single threaded
        mem_mb=get_mem_mb_align_reads,
    shell:
        "({params.dragmap} -r {params.ref_dir} -1 {input[0]} -2 {input[1]} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --num-threads {resources.n}  | samtools view -@ {resources.n} -o {output.bam}) 2> {log.dragmap_log} "
# --enable-sampling true used for (unmapped) bam input. It prevents bugs when in output bam information about whicj read is 1st or 2nd in pair.
#--preserve-map-align-order 1 was tested, so that unaligned and aligned bam have sam read order (requires thread synchronization). But reduces performance by 1/3.  Better to let mergebam job deal with the issue.

def get_mem_mb_merge_bam_alignment(wildcards, attempt):
    MEM_DEFAULT_USAGE=12000
    return attempt*(MEM_DEFAULT_USAGE)

rule merge_bam_alignment:
    """Merge bam alignment with original fastq files."""
    input:
        get_fastqpaired,
        rules.align_reads.output.bam,
        rules.adapter_removal_identify.output.stats
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.{sex}.merged.bam"),
        stats=os.path.join(config['STAT'],"{sample}.{readgroup}.{sex}.merge_stats.tsv")
    conda: "envs/pypy.yaml"
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.{sex}.mergebam.txt")
    log: os.path.join(config['LOG'],"{sample}.{readgroup}.{sex}.mergebamaligment.log")
    params:
        bam_merge=srcdir(config['BAMMERGE'])
    resources: 
        n = 2,
        mem_mb = get_mem_mb_merge_bam_alignment
    priority: 15
    shell:
        """
         (samtools view -h --threads {resources.n} {input[2]} | \
         pypy {params.bam_merge} -a  {input[0]} -b {input[1]} -s {output.stats}  |\
         samtools fixmate -@ {resources.n} -u -O BAM -m - {output.bam}) 2> {log}
        """

rule dechimer:
    """Dechimer reads.

    Dechimer is only executed if the fraction of soft-clipped bases is > DECHIMER_THRESHOLD.
    Otherwise, makes an hard link to the merged bam file.
    """
    input:
        bam=rules.merge_bam_alignment.output.bam,
        stats=rules.merge_bam_alignment.output.stats
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.{sex}.dechimer.bam"),
        stats=os.path.join(config['STAT'],"{sample}.{readgroup}.{sex}.dechimer_stats.tsv")
    priority: 16
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.{sex}.dechimer.txt")
    params:
        dechimer=srcdir(config['DECHIMER'])
    resources:
        n=4,
        mem_mb=1000
    run:
        with open(input['stats'],'r') as f:
            stats = [l for l in f.readlines()]
        primary_aligned_bp = [e.split('\t')[1].strip() for e in stats if e.startswith('primary_aligned_bp')][0]
        primary_soft_clipped_bp = [e.split('\t')[1].strip() for e in stats if e.startswith('primary_soft_clipped_bp')][0]
        res = float(primary_soft_clipped_bp) / float(primary_aligned_bp + primary_soft_clipped_bp)

        if res > float(config['DECHIMER_THRESHOLD']):
            cmd = """(samtools view -h --threads {resources.n} {input.bam} | pypy {params.dechimer} --min_align_length 40 --loose-ends -i {input.bam} -s {output.stats} |
             samtools fixmate -@ {resources.n} -u -O BAM -m - {output.bam})"""
            shell(cmd,conda_env='envs/pypy.yaml')
        else:
            cmd = """
                    #switching to copy instead of hard link, as hard link also updates modification time input.bam
                    cp {input.bam} {output.bam}
                    touch {output.stats}
                   """
            shell(cmd)


rule sort_bam_alignment:
    """Sort bam alignment by chromosome and position."""
    input:
        in_bam=rules.dechimer.output.bam
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.{sex}.sorted.bam")
    conda: "envs/preprocess.yaml"
    log:
        samtools_sort=os.path.join(config['LOG'],"{sample}.{readgroup}.{sex}.samtools_sort.log"),
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.{sex}.sort.txt")
    priority: 17
    resources:
        tmpdir=tmpdir,
        n = 2,
        mem_mb = 13000
    params:
        temp_sort=os.path.join("sort_temporary_{sample}_{readgroup}_{sex}"),
        memory_per_core= lambda wildcards, resources: int(((resources['mem_mb'] - 2000) / float(resources['n'])))        
    shell:
        """
            (samtools sort -T {resources.tmpdir}/{params.temp_sort} -@ {resources.n} -l 1 -m {params.memory_per_core}M -o {output.bam} {input}) 2> {log.samtools_sort}            
        """

rule index_sort:
    """Index sorted bam alignment."""
    input: bam = rules.sort_bam_alignment.output.bam
    output: bai = os.path.join(config['BAM'],"{sample}.{readgroup}.{sex}.sorted.bam.bai")
    conda: "envs/preprocess.yaml"
    log: samtools_index=os.path.join(config['LOG'],"{sample}.{readgroup}.{sex}.samtools_index.log"),
    priority: 18
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.{sex}.index_sorted.txt")
    resources:
        n=2,
        mem_mb=100   
    shell: "samtools index -@ {resources.n} {input.bam} 2> {log.samtools_index}"

# # function to get information about readgroups
# # needed if sample contain more than 1 fastq files
def get_readgroups_bam(wildcards):
    """Get sorted bam files for all readgroups for a given sample."""
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []
    
    for readgroup in readgroups_b:
        files.append(os.path.join(config['BAM'],wildcards['sample'] + '.' + readgroup['info']['ID'] + '.' + wildcards['sex'] + '.sorted.bam'))
    return files

def get_readgroups_bai(wildcards):
    """Get sorted bam index files for all readgroups for a given sample."""
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []
    for readgroup in readgroups_b:
        files.append(os.path.join(config['BAM'],wildcards['sample'] + '.' + readgroup['info']['ID'] + '.' + wildcards['sex'] + '.sorted.bam.bai'))
    return files

# merge different readgroups bam files for same sample
rule merge_rgs:
    """Merge bam files for different readgroups of the same sample.
    If there is only one readgroup, just link the bam file."""
    input:
        bam = get_readgroups_bam,
        bai = get_readgroups_bai,
    output:
        mer_bam=os.path.join(config['BAM'],"{sample}.{sex}.merged.bam")
    log: os.path.join(config['LOG'],"{sample}.{sex}.mergereadgroups.log")
    benchmark: "benchmark/{sample}.{sex}.merge_rgs.txt"
    resources:
        n=1,
        mem_mb=2000
    priority: 19
    run:
        if len(input.bam) > 1:
            cmd = "samtools merge -@ {resources.n} {output} {input.bam} 2> {log}"
            shell(cmd,conda_env='envs/preprocess.yaml')
        else:
            #switching to copy as hard link updates also time of input.bam
            cmd = "cp {input.bam} {output}"
            shell(cmd)

rule markdup:
    """Mark duplicates using samtools markdup."""
    input:
        bam = rules.merge_rgs.output.mer_bam
    output:
        mdbams=os.path.join(config['BAM'],"{sample}.{sex}.markdup.bam"),
        mdbams_bai=os.path.join(config['BAM'],"{sample}.{sex}.markdup.bam.bai"),
        MD_stat=os.path.join(config['STAT'],"{sample}.{sex}.markdup.stat")
    benchmark: "benchmark/{sample}.{sex}.markdup.txt"
    priority: 20
    params:
        machine=2500  #change to function
    # 100 for HiSeq and 2500 for NovaSeq
    # if fastq header not in known machines?
    # if not Illumina?
    log:
        samtools_markdup=os.path.join(config['LOG'],"{sample}.{sex}.markdup.log"),
        samtools_index_md=os.path.join(config['LOG'],"{sample}.{sex}.markdup_index.log")
    resources:
        n=4,
        mem_mb=500
    conda: "envs/preprocess.yaml"
    shell:
        """
            samtools markdup -f {output.MD_stat} -S -d {params.machine} -@ {resources.n} {input.bam} {output.mdbams} 2> {log.samtools_markdup}
            samtools index -@ {resources.n} {output.mdbams} 2> {log.samtools_index_md}
        """
    

rule mCRAM:
    """Convert bam to mapped cram."""
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai
    output:
        CRAM=os.path.join(config['CRAM'],"{sample}.{sex}.mapped_hg38.cram")
    resources:
        n=2,
        mem_mb=1000
    benchmark: os.path.join(config['BENCH'],'{sample}.{sex}.mCRAM.txt')
    priority: 30
    conda: "envs/preprocess.yaml"
    log:
        os.path.join(config['LOG'],"{sample}.{sex}.mCRAM.log")
    shell:
        "samtools view --output-fmt cram,version=3.1,archive --reference {ref} -@ {resources.n} -o {output.CRAM} {input.bam} 2> {log}"
