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
    filetype = 'fq|fastq'
    # readgroup="[\w\d_\-@]+"

from read_samples import *
from common import *

SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()
current_dir = os.getcwd()

def sampleinfo(SAMPLEINFO, sample, checkpoint=False):
    """If samples are on tape, we do not have sample readgroup info.
    This function injects a checkpoint to load this readgroup info."""

    sinfo = SAMPLEINFO[sample]
    if not 'readgroups' in sinfo:
        rgpath = os.path.join(config['SAMPLEINFODIR'], sample + ".adat")
        if os.path.exists(rgpath):
            xsample = load(rgpath)
        elif checkpoint: 
            #no readgroup info yet
            filename = checkpoints.get_readgroups.get(sample=sample).output[0]
            xsample = load(filename)
        sinfo = sinfo.copy()
        sinfo['readgroups'] = xsample['readgroups']
        sinfo['alternative_names'] = sinfo['alternative_names'].union(xsample['alternative_names'])
        SAMPLEINFO[sample] = sinfo
    return sinfo



def get_source_files(wildcards):
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
            with open(warningfile,'w') as f:
                print ("WARNING: " + w)
                for w in warnings:
                    f.write(w + '\n')
        save(sample, str(output))

rule archive_get:
    output:
        temp(os.path.join(config['FETCHDIR'], '{samplefile}.archive_{batchnr}.retrieved'))
    resources:
        arch_use_add=lambda wildcards: SAMPLEFILE_TO_BATCHES[wildcards['samplefile']]['archive'][int(wildcards['batchnr'])]['size'],
        partition="archive",
        n=1,
        mem=10
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
        shell("mkdir -p `dirname {output}`")
        shell("touch {output}")

def retrieve_batch(wildcards):
    sample = SAMPLEINFO[wildcards['sample']]
    batch = SAMPLE_TO_BATCH[wildcards['sample']]

    if sample['from_external']:
        if batch is None:
            return ancient(os.path.join(config['FETCHDIR'], wildcards['sample'] + ".finished_samples_not_assigned_to_retrieval_batch"))
        else:
            return ancient(os.path.join(config['FETCHDIR'], os.path.basename(sample['samplefile'])+ "." + str(batch) + ".retrieved"))

    else:
        return []

def calculate_active_use(wildcards):
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
    print('ACTIVE: ' + str(res))
    return res


rule start_sample:
    input:
        retrieve_batch
    output:
        temp(os.path.join(config['SOURCEDIR'], "{sample}.started"))
    resources:
        active_use_add=calculate_active_use,
        mem=1,
        n=1,
        partition="archive"
    shell: """
        touch {output}
        """

rule archive_to_active:
    input:
        retrieve_batch,
        ancient(os.path.join(config['SOURCEDIR'], "{sample}.started"))
    resources:
        arch_use_remove=lambda wildcards:   SAMPLEINFO[wildcards['sample']]['filesize'],
        partition="archive"
    output:
        temp(os.path.join(config['SOURCEDIR'], "{sample}.retrieved"))
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
            destination = os.path.join(destinationpath, e)
            shell("""
            mkdir -p `dirname {destination}`
            rsync --size-only {source} {destination}
            dmput -r -w {source}
            """)
                
        shell("""
            mkdir -p `dirname {output}`
            touch {output}
        """)
        






rule Aligner_all:
    input:
        expand("{cram}/{sample}_mapped_hg38.cram",cram=config['CRAM'],sample=sample_names),
        #expand('{stat}/{sample}.bam_all.tsv',stat=config['STAT'],sample=sample_names)
        #expand("{bams}/{sample}.markdup.bam.bai", sample=sample_names, bams=config['BAM']),
        #expand("{bams}/{sample}.merged.bam",sample=sample_names,bams=config['BAM']),
    default_target: True



#utility function to get cram reference file option for samtools
# maybe change somehow? If we don't have acces to original reference
def get_cram_ref(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    if readgroup['file_type'] == 'cram':
        cram_options = '--reference ' + readgroup['reference_file']
    else:
        cram_options = ''
    return cram_options


def get_source_aligned_file(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    #there might be multiple bam/cram files as source for a sample (e.g. if sequenced multiple times)
    #look for a read group for which the source file matches wildcard 'filename'
    #return the full file path
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if wildcards['filename'] in readgroup['info']['file']][0]
    
    if sinfo['from_external']:
        prefixpath = config['SOURCEDIR']
    else:
        prefixpath = readgroup['prefix']
   
    return os.path.join(prefixpath, readgroup['file'])


rule split_alignments_by_readgroup:
    input:
        get_source_aligned_file,
        ancient(os.path.join(config['SOURCEDIR'], "{sample}.started"))
    output:
        #there can be multiple read groups in 'filename'. Store them in this folder.        
        readgroups=temp(directory(os.path.join(config['READGROUPS'], "{sample}_{filename}")))
    resources:
        n=8,
        mem=1000
    params:
        cramref=get_cram_ref
    run:
        sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
        readgroups = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['file'].endswith(wildcards['filename'])]

        n = len(readgroups)
        if n == 1: #nothing to split, single readgropu, just link the source file
            extension = os.path.splitext(wildcards['filename'])[1]
            cmd = """
                mkdir -p {output}
                ln {input[0]} {output}/{wildcards.sample}.{readgroups[0]['info']['ID']}.{extension}
                """
            shell(cmd)                    
        else:
            #extract readgroups to output folder with {sample}.{readgroups ID}.{extension} format
            if readgroups[0]['file_type'] == 'cram': 
                #keep cram format, much more efficient, and even slighly faster
                output_fmt = 'cram,version=3.1'
            else:
                output_fmt = 'bam'

            cmd = """
                mkdir -p {output}
                samtools split -@ {resources.n} --output-fmt {output_fmt} {params.cramref} {input[0]} -f '{output}/{wildcards.sample}.%!.%.'
                """
            shell(cmd,conda_env=config['CONDA_MAIN'])



def get_aligned_readgroup_folder(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    sfile = os.path.splitext(os.path.basename(readgroup['file']))[0]
    folder = os.path.join(config['READGROUPS'], wildcards['sample'] + '_' + sfile)
    return folder



rule external_alignments_to_fastq:
    input:
        get_aligned_readgroup_folder
    output:
        fq1 = temp(config['FQ'] + "/{sample}.{readgroup}.{extension}.cut_1.fq.gz"),
        fq2 = temp(config['FQ'] + "/{sample}.{readgroup}.{extension}.cut_2.fq.gz"),
        singletons = temp(config['FQ'] + "/{sample}.{readgroup}.{extension}.extracted_singletons.fq.gz"),
        stats=os.path.join(config['STAT'],"{sample}.{readgroup}.{extension}.adapters"),
    resources:
        n: 4,
        mem: 10000
    log:
        fastq=os.path.join(config['LOG'],"{sample}.{readgroup}.{extension}.cram2fq.log"),
        adapter_removal=os.path.join(config['STAT'],"{sample}.{readgroup}.{extension}.adapter_removal.log"),
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.{extension}.cram2fq.txt")
    params:
        cramref=get_cram_ref,
        threads_ar=lambda wildcards,resources: max(int(resources['n'] / 2.0),1)
    priority: 10
    conda: config['CONDA_MAIN']
    #alternative: collate can run also in fast mode (e.g. -r 100000 -f), but this has potential impact on alignment (estimation of insert size becomes biased to genome location)
    #samtools fastq only uses ~2 threads, while collate uses much more. Therefore,pipe directly into adapterremoval to make use of the extra cores. Also save a store/load cycle.
    shell:
        """
            mkdir -p {TMPDIR}
            samtools collate -u -@ {resources.n} {cramref} -O {input}/{wildcards.sample}.{wildcards.readgroup}.{extension} {TMPDIR}/{wilcards.sample}.{wildcards.readgroup}.collate |
            samtools fastq -O -N -@ {resources.n} -0 /dev/null -s {output.singletons}  2> {log.fastq} |
            tee >(AdapterRemoval --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --interleaved-input --file1 /dev/stdin --gzip --gzip-level 1 --output1 {output.fq1} --output2 {output.fq2} --settings {log.adapter_removal} --minlength 40  --threads {params.threads_ar + 1}) |
            AdapterRemoval --identify-adapters --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --interleaved-input --file1 /dev/stdin --minalignmentlength 30 --threads {params.threads_ar - 1} > {output.stats}
        """

rule fastq_bz2togz:
    input:
        "{path}.{filetype}.bz2"
    output:
        "{path}.{filetype}.gz"
    shell: """
        bzcat {input} | bgzip > {output}
        """



#just alignment and convert to bams

def get_fastqpaired(wildcards):
    # command to extract path to fastq files from samplefile
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    # check readgroups
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    
    print('get_fastqpaired', wildcards['sample'], readgroup) 
    if sinfo['file_type'] == 'fastq_paired' or sinfo['file_type'] == 'fastq':
        file1 = readgroup['file1']
        if file1.endswith('.bz2'):
            file1 = file1[:-4] + '.gz'
        file2 = readgroup['file2']
        if file2.endswith('.bz2'):
            file2 = file2[:-4] + '.gz'
    else:
        raise RuntimeError('Unexected call to get_fastqpaired')

    return [file1, file2]

# cut adapters from inout
rule adapter_removal:
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
    threads: config["adapter_removal"]["n"]
    ##FIXME: combine adapter removal and adapter identify, use paste <(pigz -cd  test_r1cut.f1.gz | paste - - - -) <(pigz -cd test_r2cut.fq.gz | paste - - - -) |  tr '\t' '\n' |
    shell:
        """
		    AdapterRemoval --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --file1 {input[0]} --file2 {input[1]} --gzip --gzip-level 1 --output1 {output.for_f} --output2 {output.rev_f} --settings {log.adapter_removal} --minlength 40 --discarded /dev/null --threads {threads} 
		"""

rule adapter_removal_identify:
    input:
        get_fastqpaired
    output:
        stats=os.path.join(config['STAT'],"{sample}.{readgroup}.fastq.adapters"),
    priority: 5
    conda: "envs/preprocess.yaml"
    threads: config["adapter_removal_identify"]["n"]
    benchmark: os.path.join(config['BENCH'],"{sample}.{readgroup}.adapter_removal_identify.txt")
    shell:
        """
		AdapterRemoval --identify-adapters --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --file1 {input[0]} --file2 {input[1]}  --threads {threads} > {output.stats}
		"""


def get_readgroup_params(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    res = [rg for rg in sinfo['readgroups'] if rg['info']['ID'] == wildcards['readgroup']][0][
        'info']

    return {'ID': res['ID'], 'LB': res.get('LB','unknown'), 'PL': res.get('PL','unknown'),
            'PU': res.get('PU','unknown'), \
            'CN': res.get('CN','unknown'), 'DT': res.get('DT','unknown')}


# rule to align reads from cutted fq on hg38 ref
# use dragmap aligner
# samtools fixmate for future step with samtools mark duplicates

def get_mem_mb_align_reads(wildcrads, attempt):
    return (attempt - 1) * 0.5 * int(config['align_reads']['mem']) + int(config['align_reads']['mem'])

def check_sfile(wildcards):
    sfile = SAMPLEINFO[wildcards['sample']]['samplefile']
    return sfile

def get_prepared_fastq(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)

    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    if sinfo['file_type'] == 'fastq_paired' or sinfo['file_type'] == 'fastq':
        extension = 'fastq'
    elif 'cram' in sinfo['file_type']:
        extension = 'cram'
    elif 'bam' in sinfo['file_type']:
        extension = 'bam'
    else:
        raise NotImplementedError('File type not yet supported')

    file1 = os.path.join(config['FQ'], wildcards['sample'] + '.' + wildcards['readgroup'] + '.' + extension + '.cut_1.fq.gz')
    file2 = os.path.join(config['FQ'], wildcards['sample'] + '.' + wildcards['readgroup'] + '.' + extension + '.cut_2.fq.gz')
    return [file1,file2]



rule align_reads:
    input:
        get_prepared_fastq
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.aligned.bam")
    params:
        ref_dir=os.path.join(config['RES'],config['ref_dir']),
        # mask bed for current reference genome
        mask_bed=os.path.join(config['RES'],config['mask_bed']),
        temp_sort=os.path.join("sort_temporary_{sample}_{readgroup}"),
    conda: "envs/preprocess.yaml"
    threads: config["align_reads"]["n"]
    log:
        dragmap_log=os.path.join(config['LOG'],"{sample}.{readgroup}.dragmap.log"),
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.dragmap.txt")
    priority: 15
    resources:
        mem_mb=get_mem_mb_align_reads,
    shell:
        "(dragen-os -r {params.ref_dir} -1 {input[0]} -2 {input[1]} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --ht-mask-bed {params.mask_bed} --num-threads {threads}  | samtools view -@ {threads} -o {output.bam}) 2> {log.dragmap_log} "
# --enable-sampling true used for (unmapped) bam input. It prevents bugs when in output bam information about whicj read is 1st or 2nd in pair.
#--preserve-map-align-order 1 was tested, so that unaligned and aligned bam have sam read order (requires thread synchronization). But reduces performance by 1/3.  Better to let mergebam job deal with the issue.

def get_mem_mb_merge_bam_alignment(wildcrads, attempt):
    return attempt*(config['merge_bam_alignment']['mem'])

rule merge_bam_alignment:
    input:
        get_fastqpaired,
        rules.align_reads.output.bam,
        rules.adapter_removal_identify.output.stats
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.merged.bam"),
        stats=os.path.join(config['STAT'],"{sample}.{readgroup}.merge_stats.tsv")
    conda: "envs/pypy.yaml"
    threads: config["merge_bam_alignment"]["n"]
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.mergebam.txt")
    log: os.path.join(config['LOG'],"{sample}.{readgroup}.mergebamaligment.log")
    params:
        bam_merge=srcdir(config['BAMMERGE'])
    resources: mem_mb = get_mem_mb_merge_bam_alignment
    priority: 15
    shell:
        """
         (samtools view -h --threads {threads} {input[2]} | \
         pypy {params.bam_merge} -a  {input[0]} -b {input[1]} -s {output.stats}  |\
         samtools fixmate -@ {threads} -u -O BAM -m - {output.bam}) 2> {log}
        """

rule dechimer:
    input:
        bam=rules.merge_bam_alignment.output.bam,
        stats=rules.merge_bam_alignment.output.stats
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.dechimer.bam"),
        stats=os.path.join(config['STAT'],"{sample}.{readgroup}.dechimer_stats.tsv")
    priority: 16
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.dechimer.txt")
    params:
        dechimer=srcdir(config['DECHIMER'])
    threads: config['dechimer']['n']
    run:
        with open(input['stats'],'r') as f:
            stats = [l for l in f.readlines()]
        print(stats)
        primary_aligned_bp = [e.split('\t')[1].strip() for e in stats if e.startswith('primary_aligned_bp')][0]
        primary_soft_clipped_bp = [e.split('\t')[1].strip() for e in stats if e.startswith('primary_soft_clipped_bp')][0]
        res = float(primary_soft_clipped_bp) / float(primary_aligned_bp + primary_soft_clipped_bp)

        if res > float(config['DECHIMER_THRESHOLD']):
            cmd = """(samtools view -h --threads {threads} {input.bam} | pypy {params.dechimer} --min_align_length 40 --loose-ends -i {input.bam} -s {output.stats} |
             samtools fixmate -@ {threads} -u -O BAM -m - {output.bam})"""
            shell(cmd,conda_env='envs/pypy.yaml')
        else:
            cmd = """
                    ln {input.bam} {output.bam}
                    touch {output.stats}
                   """
            shell(cmd)


rule sort_bam_alignment:
    input:
        in_bam=rules.dechimer.output.bam
    output:
        bam=os.path.join(config['BAM'],"{sample}.{readgroup}.sorted.bam")
    params:
        # mask bed for current reference genome
        temp_sort=os.path.join("sort_temporary_{sample}_{readgroup}"),
        memory_per_core= int(((config["sort_bam_alignment"]["mem"] * 1000) / (config['sort_bam_alignment']['n'] * 1024))-500)
    conda: "envs/preprocess.yaml"
    threads: config["sort_bam_alignment"]["n"]
    log:
        samtools_sort=os.path.join(config['LOG'],"{sample}.{readgroup}.samtools_sort.log"),
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.sort.txt")
    priority: 17
    resources:
        tmpdir=tmpdir,
        mem_mb = config['sort_bam_alignment']['mem']
    shell:
        """
            (samtools sort -T {resources.tmpdir}/{params.temp_sort} -@ {threads} -l 1 -m {params.memory_per_core}M -o {output.bam} {input}) 2> {log.samtools_sort}            
        """

# something here (after re-running snakemake siad that below steps have to run because input was updated by above jobs)
rule index_sort:
    input: bam = ancient(rules.sort_bam_alignment.output.bam)
    output: bai = os.path.join(config['BAM'],"{sample}.{readgroup}.sorted.bam.bai")
    conda: "envs/preprocess.yaml"
    log: samtools_index=os.path.join(config['LOG'],"{sample}.{readgroup}.samtools_index.log"),
    threads: config["index_sort"]["n"]
    priority: 18
    benchmark:
        os.path.join(config['BENCH'],"{sample}.{readgroup}.index_sorted.txt")
    shell: "samtools index -@ {threads} {input.bam} 2> {log.samtools_index}"

# # function to get information about readgroups
# # needed if sample contain more than 1 fastq files
def get_readgroups_bam(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []
    for readgroup in readgroups_b:
        files.append(os.path.join(config['BAM'],wildcards['sample'] + '.' + readgroup['info']['ID'] + '.sorted.bam'))
    return files

def get_readgroups_bai(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []
    for readgroup in readgroups_b:
        files.append(os.path.join(config['BAM'],wildcards['sample'] + '.' + readgroup['info']['ID'] + '.sorted.bam.bai'))
    return files

# merge different readgroups bam files for same sample
rule merge_rgs:
    input:
        bam = get_readgroups_bam,
        bai = get_readgroups_bai,
    output:
        mer_bam=os.path.join(config['BAM'],"{sample}.merged.bam")
    log: os.path.join(config['LOG'],"{sample}.mergereadgroups.log")
    benchmark: "benchmark/{sample}.merge_rgs.txt"
    threads: config['merge_rgs']['n']
    priority: 19
    run:
        if len(input.bam) > 1:
            cmd = "samtools merge -@ {threads} {output} {input.bam} 2> {log}"
            shell(cmd,conda_env='envs/preprocess.yaml')
        else:
            cmd = "ln {input.bam} {output}"
            shell(cmd)

rule markdup:
    input:
        bam = rules.merge_rgs.output.mer_bam
    output:
        mdbams=os.path.join(config['BAM'],"{sample}.markdup.bam"),
        MD_stat=os.path.join(config['STAT'],"{sample}.markdup.stat")
    benchmark: "benchmark/{sample}.markdup.txt"
    priority: 20
    params:
        machine=2500  #change to function
    # 100 for HiSeq and 2500 for NovaSeq
    # if fastq header not in known machines?
    # if not Illumina?
    log:
        samtools_markdup=os.path.join(config['LOG'],"{sample}.markdup.log"),
    threads: config['markdup']['n']
    conda: "envs/preprocess.yaml"
    shell:
        """
            samtools markdup -f {output.MD_stat} -S -d {params.machine} -@ {threads} {input.bam} {output.mdbams} 2> {log.samtools_markdup}
        """

rule markdup_index:
    input:
        bam = rules.markdup.output.mdbams
    output:
        mdbams_bai=os.path.join(config['BAM'],"{sample}.markdup.bam.bai"),
    benchmark: "benchmark/{sample}.index_markduped.txt"
    priority: 25
    log:
        samtools_index_md=os.path.join(config['LOG'],"{sample}.markdup_index.log")
    threads: config['markdup_index']['n']
    conda: "envs/preprocess.yaml"
    shell:
        """
            samtools index -@ {threads} {input.bam} 2> {log.samtools_index_md}
        """

# mapped cram
rule mCRAM:
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup_index.output.mdbams_bai
    output:
        CRAM=os.path.join(config['CRAM'],"{sample}_mapped_hg38.cram")
    resources:
        n=2,
        mem=1000
    benchmark: os.path.join(config['BENCH'],'{sample}_mCRAM.txt')
    priority: 30
    conda: "envs/preprocess.yaml"
    log:
        os.path.join(config['LOG'],"{sample}.mCRAM.log")
    shell:
        "samtools view --output-fmt cram,version=3.1,archive --reference {ref} -@ {resources.n} -o {output.CRAM} {input.bam} 2> {log}"
