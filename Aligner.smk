import os
import os.path
import itertools
import subprocess
import time
import re
import shlex
import sys
import traceback
import json
from collections.abc import Mapping

import read_samples
from common import *
import utils

FAILED_JOBS_LOG = pj(LOG, "failed_jobs.jsonl")


def _normalize(value):
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if isinstance(value, Mapping):
        return {str(key): _normalize(val) for key, val in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_normalize(item) for item in value]
    try:
        return _normalize(dict(value))
    except Exception:
        return str(value)


def _serialize_iterable(values):
    if values is None:
        return []
    if isinstance(values, (str, os.PathLike)):
        return [_normalize(values)]
    try:
        return [_normalize(v) for v in list(values)]
    except TypeError:
        return [_normalize(values)]


def _serialize_mapping(obj):
    if obj is None:
        return {}
    if isinstance(obj, Mapping):
        return {str(key): _normalize(val) for key, val in obj.items()}
    try:
        return {str(key): _normalize(val) for key, val in dict(obj).items()}
    except Exception:
        return {"value": _normalize(obj)}


def log_failure(
    wildcards,
    input,
    output,
    params,
    log,
    threads,
    resources,
    rule_name,
    exception=None,
    **extra,
):
    if exception is None:
        exception = extra.get("exception")
    record = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "rule": rule_name,
        "wildcards": {str(key): _normalize(val) for key, val in dict(wildcards).items()},
        "input": _serialize_iterable(input),
        "output": _serialize_iterable(output),
        "params": _serialize_mapping(params),
        "log": _serialize_iterable(log),
        "threads": threads,
        "resources": _serialize_mapping(resources),
    }
    attempt = extra.get("attempt")
    if attempt is not None:
        record["attempt"] = attempt
    if exception is not None:
        record["exception"] = "".join(
            traceback.format_exception_only(type(exception), exception)
        ).strip()
    try:
        target_dir = os.path.dirname(FAILED_JOBS_LOG)
        if target_dir:
            os.makedirs(target_dir, exist_ok=True)
        with open(FAILED_JOBS_LOG, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(record, ensure_ascii=False) + "\n")
    except Exception:
        print(f"[log_failure] Failed to record failure for rule {rule_name}", file=sys.stderr)
        traceback.print_exc()


def failure_logger(rule_name):
    def _handler(wildcards, input, output, params, log, threads, resources, **extra):
        log_failure(
            wildcards,
            input,
            output,
            params,
            log,
            threads,
            resources,
            rule_name=rule_name,
            **extra,
        )

    return _handler


onsuccess: shell("rm -fr logs/Aligner/*")

wildcard_constraints:
    sample=r"[\w\d_\-@]+",
    extension=r'sam|bam|cram',
    filetype=r'fq|fastq',
    batchnr=r'[\d]+',
    readid=r'R1|R2',

# readgroup="[\w\d_\-@]+"


#ARCHIVE/DCACHE handling: it is not efficient to get files from tape file by file.
#For each file we would have to wait for the tape robot to get the tape and spin to the right position.
#Better is to first stage a batch of files together (preferably from the same tape). Once they are available,
# we can immediately copy them to active storage.
#
#Therefore 'load_samplefiles' defines batches, which are stages/copied together.
# this allows the tape robot to stage multiple files at once.
# and then to copy them to active storage

# extract all sample names from SAMPLEINFO dict to use it rule all





##THIS FUNCION IS COPIED ALSO in Stat.smk and Kraken.smk
def sampleinfo(SAMPLEINFO, sample, checkpoint=False):  #{{{
    """If samples are on tape, we do not have sample readgroup info.
    That is, the 'readgroups' field is empty.

    This function first checks if the readgroup info is available on disk,
    in the file SAMPLEINFODIR/<sample>.dat.

    Alternatively, the function injects a checkpoint rule to load this readgroup info.
    """

    sinfo = SAMPLEINFO[sample]
    if not 'readgroups' in sinfo:
        rgpath = pj(SAMPLEINFODIR,sample + ".dat")
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


#}}}

def get_source_files(wildcards):  #{{{
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
        f = append_prefix(prefixpath,f)

        if f.startswith('archive:'):
            if not archive_retrieve_add:
                files.append(ancient(pj(SOURCEDIR,wildcards['sample'] + '.archive_retrieved')))
                files.append(ancient(pj(SOURCEDIR,wildcards['sample'] + '.data')))
                archive_retrieve_add = True
        elif f.startswith('dcache:'):
            if not archive_retrieve_add:
                files.append(ancient(pj(SOURCEDIR,wildcards['sample'] + '.dcache_retrieved')))
                files.append(ancient(pj(SOURCEDIR,wildcards['sample'] + '.data')))
                archive_retrieve_add = True
        else:
            files.append(f)

    return files


#}}}

rule Aligner_all:
    input:
        expand("{cram}/{sample}.mapped_hg38.cram",sample=sample_names,cram=CRAM)   #default target removed, as it keeps all cram files on disk till end of pipeline

checkpoint get_readgroups:
    """Get the readgroup info for a sample.

    Once the checkpoint rule is executed, the readgroup info is available in the file SAMPLEINFODIR/<sample>.dat.

    Once readgroup info is available, snakemake will recalculate the DAG.
    """
    input:
        get_source_files,
        ancient(pj(SOURCEDIR,"{sample}.started"))
    output:
        temp(pj(SAMPLEINFODIR,"{sample}.dat"))
    resources:
        n="1",
        mem_mb=150
    run:
        sample = SAMPLEINFO[wildcards['sample']]
        if sample['from_external']:
            prefixpath = pj(SOURCEDIR,wildcards['sample'] + '.data')  # location where the data is downloaded from the external data repository
        else:
            prefixpath = sample['prefix']
        sample, warnings = read_samples.get_readgroups(sample,prefixpath)
        if warnings:
            warningfile = pj(SAMPLEINFODIR,wildcards['sample'] + '.warnings')
            if warnings:
                with open(warningfile,'w') as f:
                    for w in warnings:
                        print("WARNING: " + w)
                        f.write(w + '\n')
        utils.save(sample,str(output))

rule archive_get:
    """Stage a batch of files from archive.

    The batch is defined in SAMPLEFILE_TO_BATCHES.
    Once the batch is retrieved, an indicator file is written to indicate that the batch is available.
    """
    output:
        temp(pj(FETCHDIR,'{samplefile}.archive_{batchnr}.retrieved'))
    resources:
        arch_use_add=lambda wildcards:
        SAMPLEFILE_TO_BATCHES[wildcards['samplefile']]['archive'][int(wildcards['batchnr'])]['size'],
        partition="archive",
        n="0.1",
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

            files.extend([append_prefix(prefixpath, e) for e in itertools.chain(files1, files2) if e and not os.path.isabs(e)])

        # Collect archive-managed paths (strip protocol for commands)
        check_paths = [f.replace("archive:/", "") for f in files if f.startswith("archive:/")]

        # Request staging for all files in this batch (if any)
        if check_paths:
            print(f"[archive_get] Batch {wildcards['samplefile']}#{wildcards['batchnr']}: staging {len(check_paths)} file(s)")
            total = len(check_paths)
            for i in range(0, total, 100):
                chunk = check_paths[i:i+100]
                print(f"[archive_get] daget chunk {i//100 + 1}/{(total + 99)//100}: {len(chunk)} file(s)", flush=True)
                try:
                    print("RUNNING DAGET")
                    res = subprocess.run(["/opt/dacommands/bin/daget", "-av", *chunk], capture_output=True, text=True)
                    if res.stdout:
                        print(res.stdout, end="", flush=True)
                    if res.stderr:
                        print(res.stderr, end="", file=sys.stderr, flush=True)
                    if res.returncode != 0:
                        raise RuntimeError(f"[archive_get] daget exited with code {res.returncode} for this chunk")
                except Exception:
                    traceback.print_exc(file=sys.stderr)
                    raise

            # Poll until all are online (DUL or REG)
            done = set()
            while True:
                pending = [p for p in check_paths if p not in done]
                if not pending:
                    print("[archive_get] All files online. Proceeding.", flush=True)
                    break
                for i in range(0, len(pending), 100):
                    chunk = pending[i:i+100]
                    try:
                        res = subprocess.run(["/opt/dacommands/bin/dals", "-l", *chunk], capture_output=True, text=True)
                        out = (res.stdout or "") + (res.stderr or "")
                    except Exception:
                        traceback.print_exc(file=sys.stderr)
                        out = ""
                    for line in out.splitlines():
                        m = re.search(r"\(([A-Z]{3})\)\s+(.*)$", line)
                        if not m:
                            continue
                        status = m.group(1)
                        filepath = m.group(2).strip()
                        if status in ("DUL", "REG"):
                            done.add(filepath)
                print(f"[archive_get] Poll: {len(done)}/{len(check_paths)} online (DUL/REG)", flush=True)
                time.sleep(30)

        # Mark batch as retrieved only when staged
        print(f"[archive_get] Writing retrieved flag: {str(output)}", flush=True)
        try:
            os.makedirs(os.path.dirname(str(output)), exist_ok=True)
        except Exception:
            pass
        with open(str(output), "w") as _f:
            _f.write("")


def retrieve_batch(wildcards):  #{{{
    """For a sample, determines which batch needs to be retrieved from tape."""

    sample = SAMPLEINFO[wildcards['sample']]
    batch = SAMPLE_TO_BATCH[wildcards['sample']]

    if sample['from_external']:
        #batch has format <protocol>_<batchnr>  (e.g. 'archive_0')
        if batch is None:
            # sample is not assigned to a batch, this non-existing file should generate an error
            return ancient(pj(FETCHDIR,wildcards['sample'] + ".finished_samples_not_assigned_to_retrieval_batch"))
        else:
            return ancient(pj(FETCHDIR,os.path.basename(sample['samplefile']) + f".{batch}.retrieved"))

    else:
        return []


#}}}

def calculate_active_use(wildcards):  #{{{
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


#}}}

rule start_sample:
    """Start processing a sample. 

    This rule can only start if the batch to which the sample belongs has been retrieved from tape.
    This rule will reserve the space on active storage.

    Rule 'finished_sample' in Snakefile will free the active storage space.
    """
    input:
        retrieve_batch
    output:
        pj(SOURCEDIR,"{sample}.started")
    resources:
        # storsge memory is reserved for the sample
        active_use_add=calculate_active_use,
        mem_mb=50,
        n="0.1"
    shell: """
        touch {output}
        """



rule archive_to_active:
    """Move the files of a sample from archive to active storage.

    Can only start after 'start_sample' has been executed. 
    Start_sample will reserve the space on active storage.

    """
    input:
        ancient(pj(SOURCEDIR,"{sample}.started"))
    resources:
        arch_use_remove=lambda wildcards: SAMPLEINFO[wildcards['sample']]['filesize'],
        partition="archive",
        n="0.6",
        mem_mb=50
    output:
        flag=temp(pj(SOURCEDIR,"{sample}.archive_retrieved")),
        fpath=temp(directory(pj(SOURCEDIR,"{sample}.data")))
    run:
        sample = SAMPLEINFO[wildcards['sample']]
        prefixpath = sample['prefix']
        destinationpath = str(output.fpath)
        files1 = sample['file1']
        files2 = sample['file2']
        for e in itertools.chain(files1,files2):
            if not e or os.path.isabs(e):
                continue
            source = pj(prefixpath, e)
            if ':/' in source:  #remove protocol
                source = source.split(':/')[1]

            destination = pj(destinationpath,e)
            shell("""
            set -euo pipefail
            echo "[archive_to_active] RSYNC copy"
            echo "[archive_to_active] src: {source}"
            echo "[archive_to_active] dst: {destination}"
            mkdir -p "$(dirname "{destination}")"
            rsync --size-only --progress "{source}" "{destination}"
            """)
            if destination.endswith(".bz2"):
                new_filename = destination[:-4] + '.gz'
                #convert to gzip
                shell("""
                set -euo pipefail
                echo "[archive_to_active] Converting bz2->gz"
                echo "[archive_to_active] src: {destination}"
                echo "[archive_to_active] dst: {new_filename}"
                ls -l "{destination}" || true
                if command -v pigz >/dev/null 2>&1; then
                  cc="pigz -p {threads} -c"
                else
                  cc="gzip -c"
                fi
                if command -v pbzip2 >/dev/null 2>&1; then
                  dc="pbzip2 -dc -p {threads}"
                elif command -v lbzip2 >/dev/null 2>&1; then
                  dc="lbzip2 -dc -n {threads}"
                else
                  dc="bzcat"
                fi
                echo "[archive_to_active] using decompressor: $dc"
                echo "[archive_to_active] using compressor:   $cc"
                tmp="{new_filename}.tmp.$$"
                set -x
                eval "$dc" "{destination}" | eval "$cc" > "$tmp"
                set +x
                mv -f "$tmp" "{new_filename}"
                ls -l "{new_filename}" || true
                """)
                
                
        
        for e in itertools.chain(files1,files2):
            if not e or os.path.isabs(e):
                continue
            source = pj(prefixpath, e)
            if ':/' in source:  #remove protocol
                source = source.split(':/')[1]
            shell("""
            /opt/dacommands/bin/darelease "{source}" || true
            """)
        shell("""
            touch {output.flag}
        """)


def get_cram_ref(wildcards):  #{{{
    """utility function to get cram reference file option for samtools
    the read_samples utility checks that this filename is set if the file type is cram

    :return: string with cram reference file option for samtools, e.g. "--reference /path/to/ref.fa" or "" if file type is not cram
    """
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)

    wildcards = dict(wildcards)
    if 'readgroup' in wildcards:
        readgroup = \
        [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    else:
        readgroup = [readgroup for readgroup in sinfo['readgroups'] if wildcards['filename'] in readgroup['file']][0]

    if readgroup['file_type'] == 'cram':
        cram_options = '--reference ' + readgroup['reference_file']
    else:
        cram_options = ''
    return cram_options


#}}}

def ensure_source_aligned_file(wildcards):  #{{{
    """utility function to get the path to the source file for a sample.

    Takes the wildcard filename and looks up the path to that filename in the sampleinfo dictionary.
    If the sample is from external, the path is relative to the SOURCEDIR,
    otherwise it is equal to the prefixpath, which is the path of the readgroup file or the path set in the .source file.

    :return: string with path to source file
    """
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    #there might be multiple bam/cram files as source for a sample (e.g. if sequenced multiple times)
    #look for a read group for which the source file matches wildcard 'filename'
    #return the full file path

    readgroup = [readgroup for readgroup in sinfo['readgroups'] if wildcards['filename'] in readgroup['file']][0]

    #raise error if file does not exist
    result = []
    if sinfo['from_external']:
        result.append(ancient(pj(SOURCEDIR,wildcards['sample'] + '.' + sinfo['from_external'] + '_retrieved')))
        result.append(ancient(pj(SOURCEDIR,wildcards['sample'] + '.data')))


    if not os.path.exists(readgroup['file']):
        if not sinfo['from_external'] or (sinfo['from_external'] and os.path.exists(result[0])):
            raise ValueError("File does not exist: " + readgroup['file'])
    return result


#}}}

def get_mem_mb_split_alignments(wildcards, attempt):  #{{{
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    readgroups_b = sinfo['readgroups']
    if len(readgroups_b) <= 1:
        return 150
    else:
        res = 3000
    return attempt * res


#}}}

rule split_alignments_by_readgroup:
    """Split a sample bam/cram file into multiple readgroups.
    Readgroup bam/cram files are stored in the READGROUPS/<sample>_<sourcefilename> folder.

    The filenames in this folder are equal <sample>.<readgroup_id>.<extension>
    """
    input:
        ancient(pj(SOURCEDIR,"{sample}.started")),
        ensure_source_aligned_file
    output:
        #there can be multiple read groups in 'filename'. Store them in this folder.
        readgroups=temp(directory(pj(READGROUPS,"{sample}.sourcefile.{filename}"))),
        done=temp(pj(READGROUPS,"{sample}.sourcefile.{filename}.checks_done"))
    resources:
        n="1",
        mem_mb=get_mem_mb_split_alignments
    conda: CONDA_MAIN
    priority: 99
    params:
        cramref=get_cram_ref,
        fixer=srcdir('scripts/fix_bam_rg_pairs'),    
    run:
        # All branching in Python; shell executes a single, fixed command string
        sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
        readgroups = [rg for rg in sinfo['readgroups'] if wildcards['filename'] in rg['file']]

        readfile = readgroups[0]['file']
        erf_correct = SAMPLEINFO[wildcards['sample']].get('erf_correct', False)
        if erf_correct:
            print(f"[split_alignments_by_readgroup] ERF correct enabled for sample {wildcards.sample}, source {wildcards.filename} using {params.fixer}", file=sys.stderr)

        # Determine formats and reference flags
        extension_in = os.path.splitext(readfile)[1][1:].lower()
        file_type = readgroups[0]['file_type']
        reference_file = readgroups[0].get('reference_file', None)
        if file_type == 'cram':
            rflag = f"-r {reference_file}" if reference_file else ""
            output_fmt = 'cram,version=3.1'
            extension = 'cram'
        else:
            rflag = ""
            output_fmt = 'bam'
            extension = 'bam'

        sanitized = f"{output.readgroups}/{wildcards.sample}.sanitized.{extension_in}"
        n = len(readgroups)

        if n == 1:
            # Single RG path: optionally sanitize, then link
            readgroup_id = readgroups[0]['info']['ID']
            if erf_correct:
                cmd = f"""
                    set -euo pipefail
                    mkdir -p {output.readgroups}
                    {params.fixer} -i {readfile} -o {output.readgroups}/{wildcards.sample}.{readgroup_id}.{extension_in} {rflag} --threads {resources.n}                    
                    touch {output.done}
                """
            else:
                cmd = f"""
                    set -euo pipefail
                    mkdir -p {output.readgroups}
                    ln {readfile} {output.readgroups}/{wildcards.sample}.{readgroup_id}.{extension_in}
                    touch {output.done}
                """
            shell(cmd)
        else:
            # Multi-RG path: optionally sanitize, then split
            if erf_correct:
                pre = f"{params.fixer} -i {readfile} -o {sanitized} {rflag} --threads {resources.n}\n                "
                inpath = sanitized
            else:
                pre = ""
                inpath = readfile
            cmd = f"""
                set -euo pipefail
                mkdir -p {output.readgroups}
                {pre}samtools split -@ {resources.n} --output-fmt {output_fmt} {params.cramref} {inpath} -f "{output.readgroups}/{wildcards.sample}.%!.{extension}"
                touch {output.done}
            """
            shell(cmd)
            if erf_correct:
                shell(f"rm {sanitized}")


def get_aligned_readgroup_folder(wildcards):  #{{{
    """Utility function to get the path to the folder containing the readgroup files for a sample.

    This is the READGROUPS/<sample>_<sourcefilename> folder.
    """
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    sfile = os.path.splitext(os.path.basename(readgroup['file']))[0]
    folder = pj(READGROUPS,wildcards['sample'] + '.sourcefile.' + sfile)
    checkfile = pj(READGROUPS,wildcards['sample'] + '.sourcefile.' + sfile + '.checks_done')

    return [folder, checkfile]
#}}}


def get_extension(wildcards):  #{{{
    """Utility function to get the extension of the input file for a sample (bam/cram)."""
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]
    res = os.path.splitext(readgroup['file'])[1][1:].lower()

    return res
#}}}




rule external_alignments_to_fastq:
    """Convert a sample bam/cram file to fastq files.
    """
    input: get_aligned_readgroup_folder
    output:
        fq1=temp(FQ + "/{sample}.{readgroup}_R1.fastq.gz"),
        fq2=temp(FQ + "/{sample}.{readgroup}_R2.fastq.gz"),
        singletons=temp(FQ + "/{sample}.{readgroup}.extracted_singletons.fq.gz"),
    resources:
        n="1.5",
        mem_mb=lambda wildcards, attempt: (attempt - 1) * 14250 * 0.5 + 14250,
        tmpdir=tmpdir
    params:
        cramref=get_cram_ref,
        extension=get_extension,
        temp_sort=pj("external_sort_temporary_{sample}_{readgroup}_"),
        memory_per_core=6000,
        dir=lambda wildcards: get_aligned_readgroup_folder(wildcards)[0],
    priority: 10
    conda: CONDA_MAIN
    #replaced samtools collate with samtools sort due to weird memory usage behaviour of collate.
    #samtools collate -u -@ {resources.n} {params.cramref} -O {input}/{wildcards.sample}.{wildcards.readgroup}.{params.extension} {resources.tmpdir}/{wildcards.sample}.{wildcards.readgroup}.collate |
    #alternative: collate can run also in fast mode (e.g. -r 100000 -f), but this has potential impact on alignment (estimation of insert size in aligner becomes biased to genome location)
    shell:
        """
            TMP_SSD="/scratch-node/${{USER}}.${{SLURM_JOB_ID}}"
            if [ ! -d "$TMP_SSD" ] || [ ! -w "$TMP_SSD" ]; then CAND=$(ls -1dt /scratch-node/${{USER}}.* 2>/dev/null | head -n1); if [ -n "$CAND" ] && [ -d "$CAND" ] && [ -w "$CAND" ]; then TMP_SSD="$CAND"; fi; fi
            if [ -d "$TMP_SSD" ] && [ -w "$TMP_SSD" ]; then TMPDIR_USE="$TMP_SSD"; elif [ -n "$SLURM_TMPDIR" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then TMPDIR_USE="$SLURM_TMPDIR"; elif [ -d "/tmp" ] && [ -w "/tmp" ]; then TMPDIR_USE="/tmp/${{USER}}"; else TMPDIR_USE="{resources.tmpdir}"; fi
            JOB_ID="${{SLURM_JOB_ID}}"
            if [ -z "$JOB_ID" ]; then JOB_ID="${{SLURM_JOBID}}"; fi
            if [ -z "$JOB_ID" ]; then JOB_ID="$$"; fi
            JOB_TMP="$TMPDIR_USE/aligner_fastq/$JOB_ID/{wildcards.sample}.{wildcards.readgroup}"
            echo "SSD base: $TMP_SSD" >&2
            echo "TMPDIR_USE: $TMPDIR_USE" >&2
            echo "JOB_ID: $JOB_ID" >&2
            echo "JOB_TMP: $JOB_TMP" >&2
            mkdir -p "$JOB_TMP"
            trap '/bin/rm -rf "$JOB_TMP" 2>/dev/null || true' EXIT INT TERM
            CRAM_REF="{params.cramref}"
            samtools view -@ 2 -u -h $CRAM_REF {params.dir}/{wildcards.sample}.{wildcards.readgroup}.{params.extension} |\
            samtools reset -@ 2 --output-fmt BAM,level=0 --no-PG --no-RG --keep-tag OQ  |\
            samtools sort -T "$JOB_TMP"/{params.temp_sort} -@ 2 -u -n  -m {params.memory_per_core}M | \
            samtools fastq -O -N -@ 2 -0 /dev/null -1 {output.fq1} -2 {output.fq2} -s {output.singletons} 
        """


# rule fastq_bz2togz:
#     """Convert a bz2 compressed fastq file to a gz compressed fastq file.
#     filetype can be 'fq' or 'fastq'
#     """
#     input:
#         lambda wildcards: f"{path}.{filetype}.bz2" if os.path.exists(f"{path}.{filetype}.bz2") else None
#     output:
#         temp("{path}.{filetype}.gz")
#     resources:
#         n="1",
#         mem_mb=150
#     shell: """
#         bzcat {input} | bgzip > {output}
#         """


def get_fastqpaired(wildcards):  #{{{
    """Utility function to get the path to the fastq files for a sample."""

    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    # check readgroups
    readgroup = [readgroup for readgroup in sinfo['readgroups'] if readgroup['info']['ID'] == wildcards['readgroup']][0]

    if sinfo['file_type'] == 'fastq_paired' or sinfo['file_type'] == 'fastq':
        file1 = readgroup['file1']
        if file1.endswith('.bz2'):
            file1 = file1[:-4] + '.gz'
        file2 = readgroup['file2']
        if file2.endswith('.bz2'):
            file2 = file2[:-4] + '.gz'
        files = [file1, file2]
        if sinfo['from_external']:  #ensure the data folder is available if this data is retrieved from tape.
            files = files + [ancient(pj(SOURCEDIR,wildcards['sample'] + '.' + sinfo['from_external'] + '_retrieved')),
                             ancient(pj(SOURCEDIR,wildcards['sample'] + '.data'))]

    else:  #source file is a bam /cram file. We will extract fastq files with the following names:
        file1 = FQ + f"/{wildcards['sample']}.{wildcards['readgroup']}_R1.fastq.gz"
        file2 = FQ + f"/{wildcards['sample']}.{wildcards['readgroup']}_R2.fastq.gz"
    return [file1, file2]


#}}}

rule adapter_removal:
    """Remove adapters from fastq files."""
    input:
        get_fastqpaired,
        ancient(pj(SOURCEDIR,"{sample}.started"))
    output:
        for_f=temp(pj(FQ,"{sample}.{readgroup}.fastq.cut_1.fq.gz")),
        rev_f=temp(pj(FQ,"{sample}.{readgroup}.fastq.cut_2.fq.gz")),
        adapter_removal=ensure(pj(STAT,"{sample}.{readgroup}.adapter_removal.log"), non_empty=True),
        fastq_stats=pj(STAT,"{sample}.{readgroup}.fastq.stats.tsv"),
        adapters=pj(STAT,"{sample}.{readgroup}.fastq.adapters"),
    # log file in this case contain some stats about removed seqs
    priority: 10
    conda: CONDA_MAIN
    params: adapters=ADAPTERS, fastq_stats=srcdir('scripts/fastq_stats.py'),
            rmdups=srcdir('scripts/remove_interleaved_duplicates.py'), 
            rescuer=srcdir('scripts/fastq_pair_rescue.py'),
            remove_duplicated_reads=lambda wildcards: 1 if SAMPLEINFO[wildcards['sample']].get('remove_duplicated_reads', False) else 0
    resources:
        n="5",
        mem_mb=200,
        attempt=lambda wildcards, attempt: attempt
    ##FIXME: slight efficiency gain (?) if we combine adapter removal and adapter identify, use paste <(pigz -cd  test_r1cut.f1.gz | paste - - - -) <(pigz -cd test_r2cut.fq.gz | paste - - - -) |  tr '\t' '\n' |
    run:
        f1 = input[0]
        f2 = input[1]
        import gzip, bz2
        f1_q = shlex.quote(str(f1))
        f2_q = shlex.quote(str(f2))
        adapters_q = shlex.quote(str(params.adapters))
        out_fastq_stats_q = shlex.quote(str(output.fastq_stats))
        out_adapters_q = shlex.quote(str(output.adapters))
        out_for_f_q = shlex.quote(str(output.for_f))
        out_rev_f_q = shlex.quote(str(output.rev_f))
        out_adapter_removal_q = shlex.quote(str(output.adapter_removal))
        def _detect_fastq_lines(fn):
            if fn.endswith('.gz'):
                fh = gzip.open(fn, 'rt', encoding='utf-8', errors='replace')
            elif fn.endswith('.bz2'):
                fh = bz2.open(fn, 'rt', encoding='utf-8', errors='replace')
            else:
                fh = open(fn, 'rt', encoding='utf-8', errors='replace')
            try:
                l1 = fh.readline()
                l2 = fh.readline()
                l3 = fh.readline()
            finally:
                fh.close()
            if l3.startswith('+'):
                return 4
            else:
                return 2
        def _detect_quality_range(fn, max_records=1000):
            if fn.endswith('.gz'):
                fh = gzip.open(fn, 'rt', encoding='utf-8', errors='replace')
            elif fn.endswith('.bz2'):
                fh = bz2.open(fn, 'rt', encoding='utf-8', errors='replace')
            else:
                fh = open(fn, 'rt', encoding='utf-8', errors='replace')
            qmin = 10**9
            qmax = -1
            seen = 0
            try:
                while seen < max_records:
                    h = fh.readline()
                    if not h:
                        break
                    s = fh.readline()
                    p = fh.readline()
                    q = fh.readline()
                    if not q:
                        break
                    for c in q.rstrip('\n\r'):
                        oc = ord(c)
                        if oc < qmin:
                            qmin = oc
                        if oc > qmax:
                            qmax = oc
                    seen += 1
            finally:
                fh.close()
            if qmax < 0 or qmin > qmax:
                return (33, 42)
            base = 64 if qmin >= 64 else 33
            return (base, max(0, qmax - base))
        flatten1 = "- - - -"
        flatten2 = "- - - -"
        b1, m1 = _detect_quality_range(f1)
        b2, m2 = _detect_quality_range(f2)
        qualitybase = 64 if (b1 == 64 and b2 == 64) else 33
        ascii_max1 = m1 + b1
        ascii_max2 = m2 + b2
        observed_max_phred = max(ascii_max1, ascii_max2) - qualitybase
        qualitymax = 62
        qualitybase_flag = f"--qualitybase {qualitybase}" if qualitybase == 64 else ""
        qualitybase_output_flag = "--qualitybase-output 33" if qualitybase == 64 else ""
        dedup_line = f"| python3 {params.rmdups} --input - --output - --quiet " if int(params.remove_duplicated_reads) == 1 else ""
        use_rescue = 1 if int(resources.attempt) >= 2 else 0
        if use_rescue == 1:
            print(f"[adapter_removal] Using fastq_pair_rescue (attempt={int(resources.attempt)}; threads={resources.n}; external decompress)", file=sys.stderr)
            import fcntl
            err_path = str(SAMPLEINFO[wildcards['sample']]['samplefile']) + '.errors'
            with open(err_path, 'a+', encoding='utf-8') as ef:
                fcntl.flock(ef, fcntl.LOCK_EX)
                ef.seek(0)
                existing = set(line.strip() for line in ef if line.strip())
                if str(wildcards['sample']) not in existing:
                    ef.write(str(wildcards['sample']) + "\n")
                    ef.flush()
                fcntl.flock(ef, fcntl.LOCK_UN)
            src_line = f"python3 {params.rescuer} --r1 {f1_q} --r2 {f2_q} --decompress external --threads {resources.n} --buffer-size 2048 --quiet \\" 
        else:
            src_line = (                
                f"paste <(pigz -cd {f1} | paste {flatten1}) <(pigz -cd {f2} | paste {flatten2}) \\\n            | tr '\\t' '\\n' \\"
            )
        if dedup_line:
            print(f"[adapter_removal] Enabling interleaved duplicate removal", file=sys.stderr)
        cmd = f"""
            set -o pipefail
            {src_line}
            {dedup_line}| tee >(python3 {params.fastq_stats} --interleaved --input - -s {out_fastq_stats_q}) \
            | tee >(AdapterRemoval --identify-adapters --adapter-list {adapters_q} --interleaved-input --file1 /dev/stdin --threads 4 {qualitybase_flag} {qualitybase_output_flag} > {out_adapters_q}) \
            | AdapterRemoval --adapter-list {adapters_q} --interleaved-input --file1 /dev/stdin --gzip --gzip-level 1 --output1 {out_for_f_q} --output2 {out_rev_f_q} --settings {out_adapter_removal_q} --minlength 40 --singleton /dev/null --discarded /dev/null --threads 4 {qualitybase_flag} {qualitybase_output_flag} --qualitymax {qualitymax}
        """
        cmd = cmd.replace('{', '{{').replace('}', '}}')
        shell(cmd)

def get_readgroup_params(wildcards):  #{{{
    """Utility function to get the readgroup params for a sample.
       Fills in missing values with 'unknown' to avoid errors in downstream tools.
    """
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    res = [rg for rg in sinfo['readgroups'] if rg['info']['ID'] == wildcards['readgroup']][0][
        'info']

    return {'ID': res['ID'], 'LB': res.get('LB','unknown'), 'PL': res.get('PL','unknown'),
            'PU': res.get('PU','unknown'), \
            'CN': res.get('CN','unknown'), 'DT': res.get('DT','unknown')}


#}}}

def get_all_prepared_fastq(wildcards):  #{{{
    """Utility function to get the path to all (adapter-removed) fastq files for a sample (for all readgroups)."""
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []
    for readgroup in readgroups_b:
        files.append(pj(FQ,wildcards['sample'] + '.' + readgroup['info']['ID'] + '.fastq.cut_1.fq.gz'))
        files.append(pj(FQ,wildcards['sample'] + '.' + readgroup['info']['ID'] + '.fastq.cut_2.fq.gz'))
    return files


#}}}

rule kmer_reads:
    input:
        fastq=get_all_prepared_fastq
    output:
        out1=temp(pj(KMER,"{sample}.kmc_pre")),
        out2=temp(pj(KMER,"{sample}.kmc_suf")),
        lst=temp(pj(KMER,"{sample}.lst"))
    params:
        tmpdir=pj(tmpdir,"kmer_{sample}"),
        kmerdir=KMER
    conda: CONDA_KMC
    log:
        kmer_log=pj(LOG,"Aligner","{sample}.kmer.log"),
    priority: 15
    resources:
        n="2",
        mem_mb=lambda wildcards, attempt: (attempt - 1) * 0.5 * int(36000) + int(36000)
    run:
        with open(output.lst,'w') as f:
            for file in input.fastq:
                f.write(file + '\n')
        job_id = os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_JOBID') or str(os.getpid())
        kmer_dir = str(params.kmerdir)
        tmpdir_fallback = str(params.tmpdir)
        sample = str(wildcards.sample)
        list_path = str(output.lst)
        ssd_base = node_ssd_base(tmpdir_fallback)
        job_tmp = os.path.join(ssd_base, 'kmc', job_id, sample)

        shell(f"""
        mkdir -p "{kmer_dir}"
        mkdir -p "{job_tmp}"

        echo "Using temporary directory: {job_tmp}"
        echo "Temporary directory fallback: {tmpdir_fallback}"
        echo "SSD: {ssd_base}"
        echo "Using kmer directory: {kmer_dir}"
        trap '/bin/rm -rf "{job_tmp}" 2>/dev/null || true' EXIT INT TERM
        kmc -fq -k32 -cs8192 -sf12 -sp12 -sr1 -m36 @{list_path} {kmer_dir}/{sample} "{job_tmp}"
        """)

rule get_validated_sex:
    input:
        out1=pj(KMER,"{sample}.kmc_pre"),
        out2=pj(KMER,"{sample}.kmc_suf")
    output:
        yaml=temp(pj(KMER,"{sample}.result.yaml")),
        chry=temp(pj(KMER,"{sample}.chry.tsv")),
        chrx=temp(pj(KMER,"{sample}.chrx.tsv")),
        chrm=temp(pj(KMER,"{sample}.chrm.tsv")),
        auto=temp(pj(KMER,"{sample}.auto.tsv"))
    resources:
        n="0.5",
        mem_mb=lambda wildcards, attempt: (attempt - 1) * 0.5 * 4500 + 2500 if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else (attempt - 1) * 0.5 * 4500 + 2500
    params:
        kmerdir=KMER,
        process_sex=srcdir('scripts/process_sex.py')
    conda: CONDA_KMC
    shell: """
           #intersecting fastq kmers with kmers that are unique to chrY
           kmc_tools -t1 simple {params.kmerdir}/{wildcards.sample} {KMER_CHRY} -cx1 intersect {output.chry}.tmp  -ocleft
           #make sure that all kmers in chrY are present (as chrY is not always present in samples)
           kmc_tools -t1 simple {output.chry}.tmp {KMER_CHRY} union {output.chry}  -ocsum

           #intersecting fastq kmers with kmers that are unique to chrX
           kmc_tools -t1 simple {params.kmerdir}/{wildcards.sample} {KMER_CHRX} -cx1 intersect {output.chrx} -ocleft
           #intersecting fastq kmers with kmers that are unique to chrM
           kmc_tools -t1 simple {params.kmerdir}/{wildcards.sample} {KMER_CHRM} -cx1 intersect {output.chrm}  -ocleft
           #intersecting fastq kmers with kmers that are unique to autosomes
           kmc_tools -t1 simple {params.kmerdir}/{wildcards.sample} {KMER_AUTO} -cx1 intersect {output.auto} -ocleft

           #dumping kmers to tsv files            
           kmc_tools -t1 transform {output.chry} dump {output.chry}
           kmc_tools -t1 transform {output.chrx} dump {output.chrx}
           kmc_tools -t1 transform {output.chrm} dump {output.chrm}
           kmc_tools -t1 transform {output.auto} dump {output.auto}

           #removing temporary files (kmer databases)
           rm {output.chry}.*
           rm {output.chrm}.*
           rm {output.chrx}.*
           rm {output.auto}.*

           #calculate summary statistics
           python {params.process_sex} {output.auto} {output.chry} {output.chrx} {output.chrm} {output.yaml}
        """


# rule to align reads from cutted fq on hg38 ref
# use dragmap aligner
# samtools fixmate for future step with samtools mark duplicates


def get_prepared_fastq(wildcards):  #{{{
    """Utility function to get the path to the (adapter-removed) fastq files for a sample."""
    file1 = pj(FQ,wildcards['sample'] + '.' + wildcards['readgroup'] + '.fastq.cut_1.fq.gz')
    file2 = pj(FQ,wildcards['sample'] + '.' + wildcards['readgroup'] + '.fastq.cut_2.fq.gz')
    return [file1, file2]


#}}}

rule align_reads:
    """Align reads to reference genome."""
    input:
        fastq=get_prepared_fastq,
        validated_sex=rules.get_validated_sex.output.yaml
    output:
        bam=temp(pj(BAM,"{sample}.{readgroup}.aligned.bam")),
        dragmap_log=pj(STAT,"{sample}.{readgroup}.dragmap.log")            
    params:
        ref_dir=get_refdir_by_validated_sex,
        rg_params=get_readgroup_params
    conda:  CONDA_DRAGMAP
    priority: 15
    resources:
        n="22.75",#reducing thread count, as first part of dragmap is single threaded
        use_threads=24,
        mem_mb=lambda wildcards, attempt: (attempt - 1) * 0.25 * int(38000) + int(38000),
    shell:
        "(dragen-os -r {params.ref_dir} -1 {input.fastq[0]} -2 {input.fastq[1]} --RGID {wildcards.readgroup} --RGSM {wildcards.sample}  --num-threads {resources.use_threads}  | samtools view -@ 2 -o {output.bam}) 2> {output.dragmap_log} "
# --enable-sampling true used for (unmapped) bam input. It prevents bugs when in output bam information about whicj read is 1st or 2nd in pair.
#--preserve-map-align-order 1 was tested, so that unaligned and aligned bam have sam read order (requires thread synchronization). But reduces performance by 1/3.  Better to let mergebam job deal with the issue.



rule merge_bam_alignment_dechimer:
    """Merge + optional Dechimer in one streaming step.

    Runs bam_merge to restore tags and compute merge_stats, then conditionally runs
    dechimer based on primary_soft_clipped_bp_ratio. Always writes the final BAM via fixmate.
    Also produces badmap FASTQs from the merge step.
    """
    input:
        fastq=get_fastqpaired,
        bam=rules.align_reads.output.bam,
        fastq_stats=rules.adapter_removal.output.fastq_stats
    output:
        bam=temp(pj(BAM,"{sample}.{readgroup}.dechimer.bam")),
        stats=pj(STAT,"{sample}.{readgroup}.dechimer_stats.tsv"),
        badmap_fastq1=pj(FQ_BADMAP,"{sample}.{readgroup}.badmap_R1.fastq.gz"),
        badmap_fastq2=pj(FQ_BADMAP,"{sample}.{readgroup}.badmap_R2.fastq.gz"),
        merge_stats=ensure(pj(STAT,"{sample}.{readgroup}.merge_stats.tsv"),non_empty=True),
        checked=temp(pj(BAM,"{sample}.{readgroup}.bam_checked")),
        check_stats=pj(STAT,"{sample}.{readgroup}.bam_check_stats.tsv")
    priority: 16
    params:
        bam_merge=srcdir(BAMMERGE),
        dechimer=srcdir(DECHIMER),
        bam_stats_compare_hts=srcdir('scripts/bam_stats_compare_hts.py')
    resources:
        n="1.5",
        mem_mb=lambda wildcards, attempt: attempt * 5500 if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else attempt * 4500
    conda: CONDA_PYPY
    run:
        import os, shlex, tempfile
        from snakemake.shell import shell
        job_id = os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_JOBID') or str(os.getpid())
        tmpdir_fallback = os.path.dirname(str(output.bam))
        sample = str(wildcards.sample)
        tmp_dir = node_ssd_base(tmpdir_fallback)

        # Stage 1: merge + fixmate to temporary BAM, teeing to checker
        fd, merged_tmp = tempfile.mkstemp(prefix=f"{wildcards.sample}.{wildcards.readgroup}.merged.", suffix=".bam", dir=tmp_dir)
        os.close(fd)

        bami_q = shlex.quote(str(input.bam))
        fq1 = shlex.quote(str(input.fastq[0]))
        fq2 = shlex.quote(str(input.fastq[1]))
        ua = shlex.quote(str(output.badmap_fastq1))
        ub = shlex.quote(str(output.badmap_fastq2))
        merge_stats = shlex.quote(str(output.merge_stats))
        check_stats = shlex.quote(str(output.check_stats))
        checked = shlex.quote(str(output.checked))
        out_stats = shlex.quote(str(output.stats))
        outbam_final = shlex.quote(str(output.bam))
        # Build bam_merge command: run with python3 when using the Python script
        bam_merge_path = str(params.bam_merge)
        bam_merge_q = shlex.quote(bam_merge_path)
        bam_merge = bam_merge_q
        if bam_merge_path.endswith('.py'):
            bam_merge = f"python3 {bam_merge_q}"
        dechimer = shlex.quote(str(params.dechimer))
        bam_stats = shlex.quote(str(params.bam_stats_compare_hts))
        fastq_stats = shlex.quote(str(input.fastq_stats))

        cmd_stage1 = (
            "set -o pipefail; "
            f"samtools view -h --threads 2 {bami_q} "
            f"| {bam_merge} -a {fq1} -b {fq2} -ua {ua} -ub {ub} -s {merge_stats} "
            f"| tee >(python3 {bam_stats} -i - --threads 2 --fastq-stats {fastq_stats} -s {check_stats} -c {checked} > /dev/null) "
            f"| samtools fixmate -@ 2 -u -O BAM -m - {shlex.quote(merged_tmp)}"
        )
        shell(cmd_stage1)

        # Decide from merge_stats whether to run dechimer
        ratio = 0.0
        with open(str(output.merge_stats), "rt") as f:
            for line in f:
                if line.startswith("primary_soft_clipped_bp_ratio"):
                    try:
                        ratio = float(line.split("\t")[1].strip())
                    except Exception:
                        ratio = 0.0
                    break
        need_dechimer = (ratio > float(DECHIMER_THRESHOLD))

        if need_dechimer:
            fd, dechimer_tmp = tempfile.mkstemp(prefix=f"{wildcards.sample}.{wildcards.readgroup}.dechimer.", suffix=".bam", dir=tmp_dir)
            os.close(fd)
            fix_threads = 4
            cmd_stage2 = (
                "set -o pipefail; "
                f"samtools view -h --threads 2 {shlex.quote(merged_tmp)} "
                f"| {dechimer} --min_align_length 40 --loose_ends -i - -s {out_stats} "
                f"| tee >(python3 {bam_stats} -i - --threads 2 --fastq-stats {fastq_stats} -s {check_stats} -c {checked} > /dev/null) "
                f"| samtools fixmate -@ {fix_threads} -u -O BAM -m - {shlex.quote(dechimer_tmp)}"
            )
            shell(cmd_stage2)
            shell(f"mv -f {shlex.quote(dechimer_tmp)} {outbam_final}")
            try:
                os.unlink(merged_tmp)
            except Exception:
                pass
        else:
            shell(f"touch {out_stats}")
            shell(f"mv -f {shlex.quote(merged_tmp)} {outbam_final}")


rule sort_bam_alignment:
    """Sort bam alignment by chromosome and position."""
    input:
        in_bam=rules.merge_bam_alignment_dechimer.output.bam
    output:
        bam=temp(pj(BAM,"{sample}.{readgroup}.sorted.bam")),
        bai=temp(pj(BAM,"{sample}.{readgroup}.sorted.bam.bai"))
    conda: CONDA_MAIN
    log:
        samtools_sort=pj(LOG,"Aligner","{sample}.{readgroup}.samtools_sort.log"),
    priority: 17
    resources:
        tmpdir=tmpdir,
        n="1.3",
        mem_mb=13000
    params:
        temp_sort=pj("sort_temporary_{sample}_{readgroup}"),
        memory_per_core=6000
    shell:
        """
            TMP_SSD="/scratch-node/${{USER}}.${{SLURM_JOB_ID}}"
            if [ ! -d "$TMP_SSD" ] || [ ! -w "$TMP_SSD" ]; then CAND=$(ls -1dt /scratch-node/${{USER}}.* 2>/dev/null | head -n1); if [ -n "$CAND" ] && [ -d "$CAND" ] && [ -w "$CAND" ]; then TMP_SSD="$CAND"; fi; fi
            if [ -d "$TMP_SSD" ] && [ -w "$TMP_SSD" ]; then TMPDIR_USE="$TMP_SSD"; elif [ -n "$SLURM_TMPDIR" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then TMPDIR_USE="$SLURM_TMPDIR"; elif [ -d "/tmp" ] && [ -w "/tmp" ]; then TMPDIR_USE="/tmp/${{USER}}"; else TMPDIR_USE="{resources.tmpdir}"; fi
            JOB_ID="${{SLURM_JOB_ID}}"
            if [ -z "$JOB_ID" ]; then JOB_ID="${{SLURM_JOBID}}"; fi
            if [ -z "$JOB_ID" ]; then JOB_ID="$$"; fi
            JOB_TMP="$TMPDIR_USE/aligner_sort/$JOB_ID/{wildcards.sample}.{wildcards.readgroup}"
            echo "SSD base: $TMP_SSD" >&2
            echo "TMPDIR_USE: $TMPDIR_USE" >&2
            echo "JOB_ID: $JOB_ID" >&2
            echo "JOB_TMP: $JOB_TMP" >&2
            mkdir -p "$JOB_TMP"
            trap '/bin/rm -rf "$JOB_TMP" 2>/dev/null || true' EXIT INT TERM
            (samtools sort -T "$JOB_TMP"/{params.temp_sort} -@ 2 -l 1 -m {params.memory_per_core}M --write-index -o {output.bam}##idx##{output.bai} {input}) 2> {log.samtools_sort}
        """


# # function to get information about readgroups
# # needed if sample contain more than 1 fastq files
def get_readgroups_bam(wildcards):  #{{{
    """Get sorted bam files for all readgroups for a given sample."""
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []

    for readgroup in readgroups_b:
        files.append(pj(BAM,wildcards['sample'] + '.' + readgroup['info']['ID'] + '.sorted.bam'))
    return files


#}}}

def get_readgroups_bai(wildcards):  #{{{
    """Get sorted bam index files for all readgroups for a given sample."""
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []
    for readgroup in readgroups_b:
        files.append(pj(BAM,wildcards['sample'] + '.' + readgroup['info']['ID'] + '.sorted.bam.bai'))
    return files


#}}}


 


def get_readgroup_checks(wildcards):
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []
    for readgroup in readgroups_b:
        files.append(pj(BAM,wildcards['sample'] + '.' + readgroup['info']['ID'] + '.bam_checked'))
    return files



# merge different readgroups bam files for same sample
rule merge_rgs:
    """Merge bam files for different readgroups of the same sample.
    If there is only one readgroup, just link the bam file."""
    input:
        bam=get_readgroups_bam,
        bai=get_readgroups_bai,
        checks=get_readgroup_checks
    output:
        mer_bam=temp(pj(BAM,"{sample}.merged.bam"))
    log: pj(LOG,"Aligner","{sample}.mergereadgroups.log")
    resources:
        n="1",
        mem_mb=150
    priority: 19
    conda: CONDA_MAIN
    run:
        if len(input.bam) > 1:
            cmd = "samtools merge -@ {resources.n} {output} {input.bam} 2> {log}"
            shell(cmd)
        else:
            #switching to copy as hard link updates also time of input.bam
            cmd = "cp {input.bam} {output}"
            shell(cmd)


def get_badmap_fastq(wildcards):  #{{{
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    readgroups_b = sinfo['readgroups']
    files = []

    for readgroup in readgroups_b:
        files.append(pj(FQ_BADMAP,
            wildcards['sample'] + '.' + readgroup['info']['ID'] + '.badmap_' + wildcards['readid'] + '.fastq.gz'))
    return files


#}}}


rule merge_rgs_badmap:
    """Combines fastq files across readgroups from unmapped/badly mapped read for contamination check."""
    input:
        fastq=get_badmap_fastq
    output:
        fastq=temp(pj(FQ_BADMAP,"{sample}.badmap.{readid}.fastq.gz"))
    conda: CONDA_MAIN
    resources:
        n="1",
        mem_mb=150
    shell:
        """
        zcat {input.fastq} | bgzip > {output.fastq} 
        """


def get_mem_mb_markdup(wildcards, attempt):  #{{{
    res = 1500 if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else 150
    #large range of memory usage for markdup
    return (attempt - 1) * res * 3 + res


#}}}

def get_markdup_input_bam(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    rgs = sinfo['readgroups']
    if len(rgs) > 1:
        return pj(BAM, f"{wildcards['sample']}.merged.bam")
    else:
        rg_id = rgs[0]['info']['ID']
        return pj(BAM, f"{wildcards['sample']}.{rg_id}.sorted.bam")

rule markdup:
    """Mark duplicates using samtools markdup."""
    input:
        bam=get_markdup_input_bam
    output:
        mdbams=temp(pj(BAM,"{sample}.markdup.bam")),
        mdbams_bai=temp(pj(BAM,"{sample}.markdup.bam.bai")),
        MD_stat=temp(pj(STAT,"{sample}.markdup.stat"))
    priority: 20
    params:
        machine=2500,
    # machine error rate, default is 2500
    # NovaSeq uses 100
        no_dedup =lambda wildcards: 1 if SAMPLEINFO[wildcards['sample']]['no_dedup'] else 0
    log:
        samtools_markdup=pj(LOG,"Aligner","{sample}.markdup.log")
    resources:
        n="1",
        mem_mb=get_mem_mb_markdup,
        temp_loc=lambda wildcards: pj(f"markdup_temporary_{wildcards.sample}")
    conda: CONDA_MAIN
    #write index is buggy in samtools 1.17, 2/110 invalid index, probably race condition due to multithreading.
    #switching to single thread
    shell:
        """
            if [ {params.no_dedup} -eq 1 ]; then
                cp {input.bam} {output.mdbams}
                samtools index {output.mdbams}
                touch {output.MD_stat}
            else
                TMP_SSD="/scratch-node/${{USER}}.${{SLURM_JOB_ID}}"
                if [ ! -d "$TMP_SSD" ] || [ ! -w "$TMP_SSD" ]; then CAND=$(ls -1dt /scratch-node/${{USER}}.* 2>/dev/null | head -n1); if [ -n "$CAND" ] && [ -d "$CAND" ] && [ -w "$CAND" ]; then TMP_SSD="$CAND"; fi; fi
                if [ -d "$TMP_SSD" ] && [ -w "$TMP_SSD" ]; then TMPDIR_USE="$TMP_SSD"; elif [ -n "$SLURM_TMPDIR" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then TMPDIR_USE="$SLURM_TMPDIR"; else TMPDIR_USE="{resources.temp_loc}"; fi
                JOB_ID="${{SLURM_JOB_ID}}"; if [ -z "$JOB_ID" ]; then JOB_ID="${{SLURM_JOBID}}"; fi; if [ -z "$JOB_ID" ]; then JOB_ID="$$"; fi
                MDTMP="$TMPDIR_USE/markdup/$JOB_ID/{wildcards.sample}"
                mkdir -p "$(dirname "$MDTMP")"
                samtools markdup -T "$MDTMP" -f {output.MD_stat} -S -d {params.machine} {input.bam} --write-index {output.mdbams}##idx##{output.mdbams_bai} 2> {log.samtools_markdup}
            fi
        """

rule mCRAM:
    """Convert bam to mapped cram."""
    input:
        bam=rules.markdup.output.mdbams,
        bai=rules.markdup.output.mdbams_bai
    output:
        cram=temp(pj(CRAM,"{sample}.mapped_hg38.cram")),
        crai=temp(pj(CRAM,"{sample}.mapped_hg38.cram.crai"))
    resources:
        n="2",
        mem_mb=1500
    priority: 30
    conda: CONDA_MAIN
    log:
        pj(LOG,"Aligner","{sample}.mCRAM.log")
    shell:
        "samtools view --output-fmt cram,version=3.1,archive --reference {REF} -@ {resources.n} --write-index -o {output.cram}##idx##{output.crai} {input.bam} 2> {log}"




