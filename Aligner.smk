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
import tempfile
from collections.abc import Mapping

import read_samples
from common import *
import utils

FAILED_JOBS_LOG = pj(LOG, "failed_jobs.jsonl")


def _config_bool(value, name):
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {'1', 'true', 'yes', 'on'}:
        return True
    if normalized in {'0', 'false', 'no', 'off', ''}:
        return False
    raise ValueError(f"{name} must be true or false, got {value!r}")


FUSE_EXTERNAL_ADAPTER = _config_bool(
    config.get('fuse_external_adapter', True), 'fuse_external_adapter'
)
EXTERNAL_ADAPTER_LEASE_MODE = str(
    config.get('external_adapter_lease_mode', 'required')
).strip().lower()
if EXTERNAL_ADAPTER_LEASE_MODE not in {'required', 'optional', 'disabled'}:
    raise ValueError(
        "external_adapter_lease_mode must be required, optional, or disabled"
    )


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

    External inputs are materialized by the routed ``start_sample`` rule.  Its
    protocol-independent marker is written only after validation succeeds.
    """

    sinfo = SAMPLEINFO[wildcards['sample']]
    prefixpath = sinfo['prefix']
    files = []
    route_ready_added = False
    for f in itertools.chain(sinfo['file1'],sinfo['file2']):
        if not f:
            continue
        f = append_prefix(prefixpath,f)

        if f.startswith('archive:') or f.startswith('dcache:'):
            if not route_ready_added:
                files.append(ancient(pj(
                    SOURCEDIR, wildcards['sample'] + '.route_ready'
                )))
                files.append(ancient(external_data_dir(
                    wildcards['sample'], sinfo
                )))
                route_ready_added = True
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
        time = get_time('get_readgroups'),
        n="1",
        mem_mb=256
    params:
        sample=lambda wildcards: SAMPLEINFO[wildcards['sample']],
        prefixpath=lambda wildcards: (
            external_data_dir(
                wildcards['sample'], SAMPLEINFO[wildcards['sample']]
            )
            if SAMPLEINFO[wildcards['sample']]['from_external']
            else SAMPLEINFO[wildcards['sample']]['prefix']
        ),
        pipeline_root=os.path.dirname(srcdir('read_samples.py')),
        warningfile=lambda wildcards: pj(
            SAMPLEINFODIR, wildcards['sample'] + '.warnings'
        )
    conda: CONDA_MAIN
    script:
        "scripts/get_readgroups.py"

rule archive_get:
    """Stage a batch of files from archive.

    The batch is defined in SAMPLEFILE_TO_BATCHES.
    Once the batch is retrieved, an indicator file is written to indicate that the batch is available.
    """
    output:
        temp(pj(FETCHDIR,'{samplefile}.archive_{batchnr}.retrieved'))
    resources:
        time = get_time('archive_get'),
        arch_use_add=lambda wildcards:
        SAMPLEFILE_TO_BATCHES[wildcards['samplefile']]['archive'][int(wildcards['batchnr'])]['size'],
        partition="archive",
        n="0.1",
        mem_mb=256
    run:
        dname = os.path.dirname(str(output))        
        batch = SAMPLEFILE_TO_BATCHES[wildcards['samplefile']]['archive'][int(wildcards['batchnr'])]
        excluded_map = SAMPLEFILE_TO_EXCLUDED_SAMPLES.get(wildcards['samplefile'], {})
        files = []
        for sample in batch['samples']:
            sinfo = SAMPLEINFO.get(sample)
            if sinfo is None:
                if excluded_map.get(sample, {}).get('info') is not None:
                    print(f"[archive_get] Skipping excluded sample {sample}", flush=True)
                    continue
                raise KeyError(sample)
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
            last_status = {}
            poll_count = 0
            pending_report_every = int(os.environ.get("ARCHIVE_GET_PENDING_REPORT_EVERY", "10"))
            pending_report_max = int(os.environ.get("ARCHIVE_GET_PENDING_REPORT_MAX", "25"))
            while True:
                poll_count += 1
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
                        last_status[filepath] = status
                        if status in ("DUL", "REG"):
                            done.add(filepath)
                print(f"[archive_get] Poll: {len(done)}/{len(check_paths)} online (DUL/REG)", flush=True)
                if pending_report_every > 0 and (poll_count == 1 or poll_count % pending_report_every == 0):
                    pending = [p for p in check_paths if p not in done]
                    print(f"[archive_get] Pending: {len(pending)} file(s) not online yet", flush=True)
                    for p in pending[:pending_report_max]:
                        print(f"[archive_get]   ({last_status.get(p, 'UNK')}) {p}", flush=True)
                    if len(pending) > pending_report_max:
                        print(f"[archive_get]   ... and {len(pending) - pending_report_max} more", flush=True)
                time.sleep(30)

        # Mark batch as retrieved only when staged
        print(f"[archive_get] Writing retrieved flag: {str(output)}", flush=True)
        try:
            os.makedirs(os.path.dirname(str(output)), exist_ok=True)
        except Exception:
            pass
        with open(str(output), "w") as _f:
            _f.write("")


def _dcache_source_file(sample, filename):
    source_uri = append_prefix(sample['prefix'], filename)
    endpoint = read_samples.parse_dcache_uri(source_uri)
    if endpoint is None:
        raise ValueError(
            f"Expected dCache source for sample {sample['sample']}, got {source_uri!r}"
        )
    remote, remote_path = endpoint
    expected_remote = sample.get('source_remote')
    if expected_remote and remote != expected_remote:
        raise ValueError(
            f"Mixed dCache remotes for sample {sample['sample']}: "
            f"{expected_remote!r} and {remote!r}"
        )
    return remote, remote_path


def _dcache_local_destination(sample, filename, destination_root):
    relative = os.path.normpath(str(filename).lstrip('/'))
    if relative in ('', '.') or relative == '..' or relative.startswith('../'):
        raise ValueError(
            f"Unsafe dCache destination path for sample {sample['sample']}: {filename!r}"
        )
    return os.path.join(str(destination_root), relative)


rule dcache_get:
    """Stage one stable sample batch directly from Snellius."""
    output:
        temp(pj(FETCHDIR, '{samplefile}.dcache_{batchnr}.retrieved'))
    resources:
        time=get_time('dcache_get'),
        dcache_use_add=lambda wildcards:
        SAMPLEFILE_TO_BATCHES[wildcards['samplefile']]['dcache'][int(wildcards['batchnr'])]['size'],
        n="0.2",
        mem_mb=512
    params:
        transfer_script=srcdir('scripts/dcache_transfer.py')
    run:
        batch = SAMPLEFILE_TO_BATCHES[wildcards['samplefile']]['dcache'][int(wildcards['batchnr'])]
        excluded_map = SAMPLEFILE_TO_EXCLUDED_SAMPLES.get(wildcards['samplefile'], {})
        remote_paths = []
        remotes = set()
        configs = set()

        for sample_name in batch['samples']:
            sinfo = SAMPLEINFO.get(sample_name)
            if sinfo is None:
                if excluded_map.get(sample_name, {}).get('info') is not None:
                    print(f"[dcache_get] Skipping excluded sample {sample_name}", flush=True)
                    continue
                raise KeyError(sample_name)
            if not sinfo['need_retrieval']:
                continue

            remotes.add(sinfo.get('source_remote'))
            configs.add(sinfo.get('source_config'))
            for filename in itertools.chain(sinfo['file1'], sinfo['file2']):
                if not filename:
                    continue
                remote, remote_path = _dcache_source_file(sinfo, filename)
                remotes.add(remote)
                remote_paths.append(remote_path)

        remotes.discard(None)
        configs.discard(None)
        if len(remotes) != 1 or len(configs) != 1:
            raise ValueError(
                f"dCache batch {wildcards.samplefile}#{wildcards.batchnr} must use "
                f"one remote/config; got remotes={sorted(remotes)} configs={sorted(configs)}"
            )

        os.makedirs(FETCHDIR, exist_ok=True)
        fd, list_path = tempfile.mkstemp(
            prefix=f".{wildcards.samplefile}.dcache-stage-",
            suffix=".txt",
            dir=FETCHDIR,
            text=True,
        )
        try:
            with os.fdopen(fd, 'w', encoding='utf-8') as handle:
                for remote_path in remote_paths:
                    handle.write(remote_path + '\n')

            if remote_paths:
                subprocess.run(
                    [
                        sys.executable,
                        str(params.transfer_script),
                        'stage',
                        '--config',
                        next(iter(configs)),
                        '--remote',
                        next(iter(remotes)),
                        '--file-list',
                        list_path,
                        '--lifetime',
                        str(config.get('dcache_stage_lifetime', '7D')),
                        '--poll-seconds',
                        str(config.get('dcache_stage_poll_seconds', 60)),
                        '--stage-timeout',
                        str(config.get('dcache_stage_timeout', 86400)),
                    ],
                    check=True,
                )
        finally:
            try:
                os.unlink(list_path)
            except FileNotFoundError:
                pass

        os.makedirs(os.path.dirname(str(output[0])), exist_ok=True)
        with open(str(output[0]), 'w', encoding='utf-8'):
            pass


def retrieve_batch(wildcards):  #{{{
    """Return the batch-stage marker needed before routed materialization."""

    sample = SAMPLEINFO[wildcards['sample']]
    route = _start_sample_route(wildcards)
    if route == 'active':
        return []

    route_ready = pj(SOURCEDIR, wildcards['sample'] + '.route_ready')
    legacy_ready = pj(SOURCEDIR, wildcards['sample'] + f'.{route}_retrieved')
    destination = external_data_dir(wildcards['sample'], sample)
    if os.path.exists(route_ready) or (
        os.path.exists(legacy_ready) and os.path.isdir(destination)
    ):
        # Adopt a complete old-layout materialization without restaging its
        # whole stable batch merely to create the new universal marker.
        return []

    batch = SAMPLE_TO_BATCH[wildcards['sample']]
    # batch has format <protocol>_<batchnr> (for example archive_0)
    if batch is None:
        return ancient(pj(
            FETCHDIR,
            wildcards['sample'] +
            ".finished_samples_not_assigned_to_retrieval_batch",
        ))
    else:
        return ancient(pj(
            FETCHDIR,
            os.path.basename(sample['samplefile']) + f".{batch}.retrieved",
        ))


#}}}

def _start_sample_route(wildcards):
    route = SAMPLEINFO[wildcards['sample']].get('from_external')
    return str(route).lower() if route else 'active'


def _start_sample_time(wildcards, attempt=1):
    rule_name = {
        'active': 'start_sample',
        'archive': 'archive_to_active',
        'dcache': 'dcache_to_active',
    }[_start_sample_route(wildcards)]
    return get_time(rule_name)(wildcards, attempt)


def _start_sample_partition(wildcards):
    return 'archive' if _start_sample_route(wildcards) == 'archive' else 'compute'


def _start_sample_cores(wildcards):
    return {'active': 0.1, 'archive': 0.6, 'dcache': 1.0}[
        _start_sample_route(wildcards)
    ]


def _start_sample_mem_mb(wildcards):
    return 1024 if _start_sample_route(wildcards) == 'dcache' else 256


def _start_sample_active_add(wildcards):
    started = pj(SOURCEDIR, wildcards['sample'] + '.started')
    route_ready = pj(SOURCEDIR, wildcards['sample'] + '.route_ready')
    if os.path.exists(started) and not os.path.exists(route_ready):
        # One-time adoption of the old two-step layout: its start job already
        # took the old lifecycle reservation. The old active-input formula did
        # not include the source bytes, so add exactly that migration delta for
        # an unfinished active sample. External routes already included it.
        sample = SAMPLEINFO[wildcards['sample']]
        finished = pj(SOURCEDIR, wildcards['sample'] + '.finished')
        if not sample.get('from_external') and not os.path.exists(finished):
            return sample['filesize']
        return 0
    return active_use_gb(wildcards)


def _start_sample_tier_remove(wildcards, route):
    if _start_sample_route(wildcards) != route:
        return 0
    legacy = pj(SOURCEDIR, wildcards['sample'] + f'.{route}_retrieved')
    route_ready = pj(SOURCEDIR, wildcards['sample'] + '.route_ready')
    finished = pj(SOURCEDIR, wildcards['sample'] + '.finished')
    if (
        os.path.exists(legacy)
        or os.path.exists(route_ready)
        or os.path.exists(finished)
    ):
        # Do not release durable accounting twice while adopting an old run.
        return 0
    return SAMPLEINFO[wildcards['sample']]['filesize']


def _start_sample_pattern(*routes):
    samples = sorted(
        re.escape(sample)
        for sample in SAMPLEINFO
        if _start_sample_route({'sample': sample}) in routes
    )
    return '(?:' + '|'.join(samples) + ')' if samples else r'(?!)'


START_SAMPLE_ACTIVE_PATTERN = _start_sample_pattern('active')
START_SAMPLE_ARCHIVE_PATTERN = _start_sample_pattern('archive')
START_SAMPLE_DCACHE_PATTERN = _start_sample_pattern('dcache')
START_SAMPLE_EXTERNAL_PATTERN = _start_sample_pattern('archive', 'dcache')


def _run_start_sample(wildcards, output, params):
    helper_dir = os.path.dirname(str(params.job_helper))
    if helper_dir not in sys.path:
        sys.path.insert(0, helper_dir)
    from start_sample_job import run_start_sample_job

    run_start_sample_job(
        sample=SAMPLEINFO[wildcards['sample']],
        sample_name=str(wildcards.sample),
        expected_route=str(params.expected_route),
        destination=(
            str(output.materialized)
            if hasattr(output, 'materialized')
            else None
        ),
        started=str(output.started),
        route_ready=str(output.route_ready),
        append_prefix=append_prefix,
        dcache_source_file=_dcache_source_file,
        transfer_script=str(params.transfer_script),
        source_dir=SOURCEDIR,
        dcache_download_workers=max(
            1, int(config.get('dcache_download_workers', 2))
        ),
        dcache_download_lock_slots=max(1, int(config.get(
            'dcache_download_lock_slots',
            config.get('dcache_transfer_slots', 4),
        ))),
    )


#}}}

rule start_sample_active:
    """Reserve the sample lifecycle after validating active source input."""
    input:
        retrieve_batch
    output:
        started=pj(SOURCEDIR,"{sample}.started"),
        route_ready=pj(SOURCEDIR,"{sample}.route_ready")
    wildcard_constraints:
        sample=START_SAMPLE_ACTIVE_PATTERN
    resources:
        time=_start_sample_time,
        partition=_start_sample_partition,
        active_use_add=_start_sample_active_add,
        arch_use_remove=0,
        dcache_use_remove=0,
        dcache_download_slots=0,
        mem_mb=_start_sample_mem_mb,
        n=_start_sample_cores
    params:
        expected_route='active',
        job_helper=srcdir('scripts/start_sample_job.py'),
        transfer_script=srcdir('scripts/dcache_transfer.py')
    run:
        _run_start_sample(wildcards, output, params)


rule start_sample_archive:
    """Reserve and materialize one archive sample in a tracked temp directory."""
    input:
        retrieve_batch
    output:
        started=pj(SOURCEDIR,"{sample}.started"),
        route_ready=temp(pj(SOURCEDIR,"{sample}.route_ready")),
        materialized=temp(directory(pj(SOURCEDIR,"{sample}.data")))
    wildcard_constraints:
        sample=START_SAMPLE_ARCHIVE_PATTERN
    resources:
        time=_start_sample_time,
        partition=_start_sample_partition,
        active_use_add=_start_sample_active_add,
        arch_use_remove=lambda wildcards: _start_sample_tier_remove(
            wildcards, 'archive'
        ),
        dcache_use_remove=0,
        dcache_download_slots=0,
        mem_mb=_start_sample_mem_mb,
        n=_start_sample_cores
    params:
        expected_route='archive',
        job_helper=srcdir('scripts/start_sample_job.py'),
        transfer_script=srcdir('scripts/dcache_transfer.py')
    run:
        _run_start_sample(wildcards, output, params)


rule start_sample_dcache:
    """Reserve and download one dCache sample in a tracked temp directory."""
    input:
        retrieve_batch
    output:
        started=pj(SOURCEDIR,"{sample}.started"),
        route_ready=temp(pj(SOURCEDIR,"{sample}.route_ready")),
        materialized=temp(directory(pj(SOURCEDIR,"{sample}.dcache_data")))
    wildcard_constraints:
        sample=START_SAMPLE_DCACHE_PATTERN
    resources:
        time=_start_sample_time,
        partition=_start_sample_partition,
        active_use_add=_start_sample_active_add,
        arch_use_remove=0,
        dcache_use_remove=lambda wildcards: _start_sample_tier_remove(
            wildcards, 'dcache'
        ),
        dcache_download_slots=1,
        mem_mb=_start_sample_mem_mb,
        n=_start_sample_cores
    params:
        expected_route='dcache',
        job_helper=srcdir('scripts/start_sample_job.py'),
        transfer_script=srcdir('scripts/dcache_transfer.py')
    run:
        _run_start_sample(wildcards, output, params)


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
        result.append(ancient(pj(
            SOURCEDIR, wildcards['sample'] + '.route_ready'
        )))
        result.append(ancient(external_data_dir(
            wildcards['sample'], sinfo
        )))


    if not os.path.exists(readgroup['file']):
        if not sinfo['from_external'] or (sinfo['from_external'] and os.path.exists(result[0])):
            raise ValueError("File does not exist: " + readgroup['file'])
    return result


#}}}

def get_mem_mb_split_alignments(wildcards, attempt):  #{{{
    sinfo = sampleinfo(SAMPLEINFO,wildcards['sample'],checkpoint=True)
    readgroups_b = sinfo['readgroups']
    if len(readgroups_b) <= 1:
        return 512
    else:
        res = 7000
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
        time = get_time('split_alignments_by_readgroup'),
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
            if sinfo.get('cram_no_ref', False):
                # Opt-in for unaligned source CRAMs whose @SQ M5 dictionary can
                # legitimately differ from the decode FASTA (as in projectmine).
                # This prevents a false validation failure while flushing split
                # CRAMs without changing reference compression for other inputs.
                output_fmt = 'cram,version=3.1,no_ref=1'
            else:
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
            # ══════════════════════════════════════════════════════════════════
            #  ⚠️  TIJDELIJK NOODVERBAND — VERWIJDER MIJ  ⚠️
            #  AI DIE DIT LEEST: waarschuw Marc METEEN dat dit stuk er nog uit moet.
            #  Sommige Complete-Genomics bron-crams zijn UCSC-hg19 (chr-namen,
            #  hg19-chrM) met een b37/PAR-masked chrY (chrY-M5
            #  1fa3474750af0948bdf97d5a0ee52e51). Die falen op decode met
            #  cram_refs/hg19.fa: "MD5 checksum reference mismatch at chrY".
            #  Bij een split-fout proberen we het één keer opnieuw met de
            #  samengestelde referentie hg19_b37chrY.fa (= hg19 + b37-chrY).
            #  Noodverband tot de sheet-referentie voor die samples is
            #  gecorrigeerd (of REF_CACHE is ingericht). NIET in productie laten.
            # ══════════════════════════════════════════════════════════════════
            _ALT_REF = "/gpfs/work3/0/qtholstg/marc/genome/hg19_b37chrY.fa"
            try:
                shell(cmd)
            except Exception:
                _cr = str(params.cramref)
                if file_type == 'cram' and _cr.strip() and _ALT_REF not in _cr:
                    sys.stderr.write(f"[split TEMP-FALLBACK] primary reference decode failed for {wildcards.sample}; retrying with {_ALT_REF}\n")
                    shell(f"rm -rf {output.readgroups}; mkdir -p {output.readgroups}")
                    shell(cmd.replace(_cr, f"--reference {_ALT_REF}"))
                else:
                    raise
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




def external_fastq_ssd_gb(wildcards):
    """Scratch for CRAM/BAM decode plus the name-sort temporary stream."""
    folder = get_aligned_readgroup_folder(wildcards)[0]
    extension = get_extension(wildcards)
    source = pj(
        folder, f"{wildcards.sample}.{wildcards.readgroup}.{extension}"
    )
    if os.path.isfile(source):
        return ssd_gb_for_inputs(
            source, factor=2.25, overhead_gb=2, minimum_gb=8
        )

    # During initial DAG construction the checkpoint-generated split BAM/CRAM
    # does not exist yet.  Use the whole sample source size as a conservative
    # upper bound for one read group; once present, retries use the exact file.
    source_gb_upper_bound = float(
        sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)['filesize']
    )
    return max(8, int(math.ceil(source_gb_upper_bound * 2.25 + 2)))


rule external_alignments_to_fastq:
    """Convert a sample bam/cram file to fastq files.
    """
    input: get_aligned_readgroup_folder
    output:
        fq1=temp(FQ + "/{sample}.{readgroup}_R1.fastq.gz"),
        fq2=temp(FQ + "/{sample}.{readgroup}_R2.fastq.gz"),
        singletons=temp(FQ + "/{sample}.{readgroup}.extracted_singletons.fq.gz"),
    resources:
        time = get_time('external_alignments_to_fastq'),
        n="1.5",
        mem_mb=lambda wildcards, attempt: (attempt - 1) * 14250 * 0.5 + 14250,
        tmpdir=tmpdir,
        ssd_use="required",
        ssd_gb=external_fastq_ssd_gb
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
            if [ ! -d "$TMP_SSD" ] || [ ! -w "$TMP_SSD" ]; then CAND=$(ls -1dt /scratch-node/${{USER}}.* 2>/dev/null | head -n1 || true); if [ -n "${{CAND:-}}" ] && [ -d "$CAND" ] && [ -w "$CAND" ]; then TMP_SSD="$CAND"; fi; fi
            if [ -d "$TMP_SSD" ] && [ -w "$TMP_SSD" ]; then TMPDIR_USE="$TMP_SSD"; elif [ -n "${{SLURM_TMPDIR:-}}" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then TMPDIR_USE="$SLURM_TMPDIR"; else TMPDIR_USE="{resources.tmpdir}"; fi
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
            files.append(ancient(pj(
                SOURCEDIR, wildcards['sample'] + '.route_ready'
            )))
            files.append(ancient(external_data_dir(
                wildcards['sample'], sinfo
            )))

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
        time = get_time('adapter_removal'),
        n="5",
        mem_mb=512,
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


def external_adapter_ssd_gb(wildcards):
    folder = get_aligned_readgroup_folder(wildcards)[0]
    extension = get_extension(wildcards)
    source = pj(folder, f"{wildcards.sample}.{wildcards.readgroup}.{extension}")
    if os.path.isfile(source):
        return ssd_gb_for_inputs(
            source, factor=4.5, overhead_gb=8, minimum_gb=32
        )
    source_gb_upper_bound = float(
        sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)['filesize']
    )
    return max(32, int(math.ceil(source_gb_upper_bound * 4.5 + 8)))


def external_alignment_path(wildcards):
    folder = get_aligned_readgroup_folder(wildcards)[0]
    extension = get_extension(wildcards)
    return pj(folder, f"{wildcards.sample}.{wildcards.readgroup}.{extension}")


EXTERNAL_ALIGNMENT_SAMPLE_PATTERN = '(?:' + '|'.join(
    re.escape(sample)
    for sample, sinfo in SAMPLEINFO.items()
    if sinfo.get('file_type') not in {'fastq', 'fastq_paired'}
) + ')'
if EXTERNAL_ALIGNMENT_SAMPLE_PATTERN == '(?:)':
    EXTERNAL_ALIGNMENT_SAMPLE_PATTERN = r'(?!)'


if FUSE_EXTERNAL_ADAPTER:
    rule external_adapter_fused:
        """Extract BAM/CRAM FASTQs once and remove adapters on assigned SSD."""
        input:
            aligned=get_aligned_readgroup_folder,
            started=ancient(pj(SOURCEDIR,"{sample}.started"))
        output:
            raw_fq1=temp(pj(FQ,"{sample}.{readgroup}_R1.fastq.gz")),
            raw_fq2=temp(pj(FQ,"{sample}.{readgroup}_R2.fastq.gz")),
            singletons=temp(pj(FQ,"{sample}.{readgroup}.extracted_singletons.fq.gz")),
            for_f=temp(pj(FQ,"{sample}.{readgroup}.fastq.cut_1.fq.gz")),
            rev_f=temp(pj(FQ,"{sample}.{readgroup}.fastq.cut_2.fq.gz")),
            adapter_removal=ensure(pj(STAT,"{sample}.{readgroup}.adapter_removal.log"), non_empty=True),
            fastq_stats=pj(STAT,"{sample}.{readgroup}.fastq.stats.tsv"),
            adapters=pj(STAT,"{sample}.{readgroup}.fastq.adapters")
        wildcard_constraints:
            sample=EXTERNAL_ALIGNMENT_SAMPLE_PATTERN
        log:
            runner=pj(LOG,"Aligner","{sample}.{readgroup}.external_adapter_fused.log"),
            io_profile=pj(LOG,"Aligner","{sample}.{readgroup}.external_adapter_fused.io.json")
        params:
            runner=srcdir('scripts/run_fused_external_adapter.py'),
            alignment=external_alignment_path,
            cram_options=get_cram_ref,
            adapters=ADAPTERS,
            fastq_stats=srcdir('scripts/fastq_stats.py'),
            rmdups=srcdir('scripts/remove_interleaved_duplicates.py'),
            rescuer=srcdir('scripts/fastq_pair_rescue.py'),
            remove_duplicated_reads=lambda wc: int(SAMPLEINFO[wc.sample].get('remove_duplicated_reads', False)),
            error_file=lambda wc: str(SAMPLEINFO[wc.sample]['samplefile']) + '.errors',
            lease_mode=EXTERNAL_ADAPTER_LEASE_MODE,
            lease_command=zslurm_lease_command(config)
        conda: CONDA_MAIN
        priority: 10
        resources:
            time=get_time('external_adapter_fused'),
            n="5",
            mem_mb=lambda wildcards, attempt: (
                (attempt - 1) * 14250 * 0.5 + 14250
            ),
            attempt=lambda wildcards, attempt: attempt,
            ssd_use="required",
            ssd_gb=external_adapter_ssd_gb
        shell:
            """
            python {params.runner:q} \
                --input-alignment {params.alignment:q} \
                --cram-options {params.cram_options:q} \
                --sample {wildcards.sample:q} \
                --readgroup {wildcards.readgroup:q} \
                --adapter-list {params.adapters:q} \
                --fastq-stats-script {params.fastq_stats:q} \
                --remove-duplicates-script {params.rmdups:q} \
                --pair-rescue-script {params.rescuer:q} \
                --remove-duplicated-reads {params.remove_duplicated_reads} \
                --attempt {resources.attempt} \
                --error-file {params.error_file:q} \
                --output-raw-forward {output.raw_fq1:q} \
                --output-raw-reverse {output.raw_fq2:q} \
                --output-singletons {output.singletons:q} \
                --output-forward {output.for_f:q} \
                --output-reverse {output.rev_f:q} \
                --output-adapter-log {output.adapter_removal:q} \
                --output-fastq-stats {output.fastq_stats:q} \
                --output-adapters {output.adapters:q} \
                --metrics {log.io_profile:q} \
                --initial-cores {resources.n} \
                --initial-memory-mb {resources.mem_mb} \
                --adapter-cores 5 \
                --adapter-memory-mb 1024 \
                --lease-mode {params.lease_mode:q} \
                --lease-command {params.lease_command:q} \
                --ssd-gb {resources.ssd_gb} \
                2> {log.runner:q}
            """

    ruleorder: external_adapter_fused > external_alignments_to_fastq > adapter_removal

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


FUSE_KMER_SEX = _config_bool(
    config.get('fuse_kmer_sex', True), 'fuse_kmer_sex'
)
KMER_SEX_LEASE_MODE = str(
    config.get('kmer_sex_lease_mode', 'required')
).strip().lower()
if KMER_SEX_LEASE_MODE not in {'required', 'optional', 'disabled'}:
    raise ValueError(
        "kmer_sex_lease_mode must be required, optional, or disabled"
    )

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
        time = get_time('kmer_reads'),
        n="2",
        mem_mb=lambda wildcards, attempt: (attempt - 1) * 0.5 * int(36000) + int(36000),
        ssd_use="required",
        ssd_gb=lambda wildcards, input: ssd_gb_for_inputs(input.fastq, factor=1.0, overhead_gb=3, minimum_gb=8)
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
        time = get_time('get_validated_sex'),
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


if FUSE_KMER_SEX:
    rule kmer_sex_fused:
        """Build the KMC database and validate sex without GPFS intermediates."""
        input:
            fastq=get_all_prepared_fastq
        output:
            yaml=temp(pj(KMER,"{sample}.result.yaml")),
            chry=temp(pj(KMER,"{sample}.chry.tsv")),
            chrx=temp(pj(KMER,"{sample}.chrx.tsv")),
            chrm=temp(pj(KMER,"{sample}.chrm.tsv")),
            auto=temp(pj(KMER,"{sample}.auto.tsv"))
        log:
            runner=pj(LOG,"Aligner","{sample}.kmer_sex_fused.log"),
            io_profile=pj(LOG,"Aligner","{sample}.kmer_sex_fused.io.json")
        params:
            runner=srcdir('scripts/run_fused_kmer_sex.py'),
            process_sex=srcdir('scripts/process_sex.py'),
            kmer_chry=KMER_CHRY,
            kmer_chrx=KMER_CHRX,
            kmer_chrm=KMER_CHRM,
            kmer_auto=KMER_AUTO,
            lease_mode=KMER_SEX_LEASE_MODE,
            lease_command=zslurm_lease_command(config)
        conda: CONDA_KMC
        priority: 15
        resources:
            time=get_time('kmer_sex_fused'),
            n="2",
            mem_mb=lambda wildcards, attempt: (
                (attempt - 1) * 0.5 * 42000 + 42000
            ),
            ssd_use="required",
            ssd_gb=lambda wildcards, input: ssd_gb_for_inputs(
                input.fastq, factor=3.0, overhead_gb=8, minimum_gb=32
            )
        shell:
            """
            python {params.runner:q} \
                --fastq {input.fastq:q} \
                --sample {wildcards.sample:q} \
                --output-yaml {output.yaml:q} \
                --output-chry {output.chry:q} \
                --output-chrx {output.chrx:q} \
                --output-chrm {output.chrm:q} \
                --output-auto {output.auto:q} \
                --kmer-chry {params.kmer_chry:q} \
                --kmer-chrx {params.kmer_chrx:q} \
                --kmer-chrm {params.kmer_chrm:q} \
                --kmer-auto {params.kmer_auto:q} \
                --process-sex {params.process_sex:q} \
                --metrics {log.io_profile:q} \
                --initial-cores {resources.n} \
                --initial-memory-mb {resources.mem_mb} \
                --low-cores 0.5 \
                --low-memory-mb 3000 \
                --lease-mode {params.lease_mode:q} \
                --lease-command {params.lease_command:q} \
                --ssd-gb {resources.ssd_gb} \
                2> {log.runner:q}
            """

    ruleorder: kmer_sex_fused > get_validated_sex


# rule to align reads from cutted fq on hg38 ref
# use dragmap aligner
# samtools fixmate for future step with samtools mark duplicates


def get_prepared_fastq(wildcards):  #{{{
    """Utility function to get the path to the (adapter-removed) fastq files for a sample."""
    file1 = pj(FQ,wildcards['sample'] + '.' + wildcards['readgroup'] + '.fastq.cut_1.fq.gz')
    file2 = pj(FQ,wildcards['sample'] + '.' + wildcards['readgroup'] + '.fastq.cut_2.fq.gz')
    return [file1, file2]


#}}}


FUSE_ALIGNMENT_PHASES = _config_bool(
    config.get('fuse_alignment_phases', True), 'fuse_alignment_phases'
)
ALIGNMENT_LEASE_MODE = str(
    config.get('alignment_lease_mode', 'required')
).strip().lower()
if ALIGNMENT_LEASE_MODE not in {'required', 'optional', 'disabled'}:
    raise ValueError(
        "alignment_lease_mode must be required, optional, or disabled"
    )


def _fused_low_memory_mb(wildcards):
    # A 179-GB production readgroup made bam_merge reach about 15.3 GB RSS;
    # coordinate sort has also peaked around 13.5 GB. Keep enough headroom for
    # the entire merge/dechimer/sort tail at one monotonic target so the job
    # never needs to reacquire memory after releasing the 40-GB align lease.
    return 20000


def _fused_ignore_qual_flag(wildcards):
    return (
        '--ignore-qual-checksum-diff'
        if bool(SAMPLEINFO[wildcards['sample']].get('erf_correct', False))
        else ''
    )

rule align_reads:
    """Align reads to reference genome."""
    input:
        fastq=get_prepared_fastq,
        validated_sex=pj(KMER,"{sample}.result.yaml")
    output:
        bam=temp(pj(BAM,"{sample}.{readgroup}.aligned.bam")),
        dragmap_log=pj(STAT,"{sample}.{readgroup}.dragmap.log")            
    params:
        ref_dir=get_refdir_by_validated_sex,
        rg_params=get_readgroup_params
    conda:  CONDA_DRAGMAP
    priority: 15
    resources:
        time = get_time('align_reads'),
        n="22.75",#reducing thread count, as first part of dragmap is single threaded
        use_threads=24,
        mem_mb=lambda wildcards, attempt: (attempt - 1) * 0.25 * int(40000) + int(40000),
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
        fastq_stats=pj(STAT,"{sample}.{readgroup}.fastq.stats.tsv")
    output:
        bam=temp(pj(BAM,"{sample}.{readgroup}.dechimer.bam")),
        stats=pj(STAT,"{sample}.{readgroup}.dechimer_stats.tsv"),
        badmap_fastq1=temp(pj(FQ_BADMAP,"{sample}.{readgroup}.badmap_R1.fastq.gz")),
        badmap_fastq2=temp(pj(FQ_BADMAP,"{sample}.{readgroup}.badmap_R2.fastq.gz")),
        merge_stats=ensure(pj(STAT,"{sample}.{readgroup}.merge_stats.tsv"),non_empty=True),
        checked=temp(pj(BAM,"{sample}.{readgroup}.bam_checked")),
        check_stats=pj(STAT,"{sample}.{readgroup}.bam_check_stats.tsv")
    priority: 16
    params:
        bam_merge=srcdir(BAMMERGE),
        dechimer=srcdir(DECHIMER),
        bam_stats_compare_hts=srcdir('scripts/bam_stats_compare_hts.py')
    resources:
        time = get_time('merge_bam_alignment_dechimer'),
        n="1.5",
        mem_mb=lambda wildcards, attempt: attempt * 5500 if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else attempt * 4500,
        ssd_use="required",
        ssd_gb=lambda wildcards, input: ssd_gb_for_inputs(input.bam, factor=2.25, overhead_gb=2, minimum_gb=8)
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
        ignore_qual_flag = "--ignore-qual-checksum-diff" if bool(SAMPLEINFO[wildcards['sample']].get('erf_correct', False)) else ""

        cmd_stage1 = (
            "set -o pipefail; "
            f"samtools view -h --threads 2 {bami_q} "
            f"| {bam_merge} -a {fq1} -b {fq2} -ua {ua} -ub {ub} -s {merge_stats} "
            f"| tee >(python3 {bam_stats} -i - --threads 2 --fastq-stats {fastq_stats} {ignore_qual_flag} -s {check_stats} -c {checked} > /dev/null) "
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
                f"| tee >(python3 {bam_stats} -i - --threads 2 --fastq-stats {fastq_stats} {ignore_qual_flag} -s {check_stats} -c {checked} > /dev/null) "
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


if FUSE_ALIGNMENT_PHASES:
    rule align_reads_fused:
        """Align, merge, conditionally dechimer, check, and sort one read group.

        DRAGMAP and its BAM are node-local.  After alignment the job returns
        CPU and memory through an absolute ZSlurm lease target. Coordinate
        sort runs before job completion at that same low-resource target.
        """
        input:
            prepared_fastq=get_prepared_fastq,
            validated_sex=pj(KMER,"{sample}.result.yaml"),
            source_fastq=get_fastqpaired,
            fastq_stats=pj(STAT,"{sample}.{readgroup}.fastq.stats.tsv")
        output:
            bam=temp(pj(BAM,"{sample}.{readgroup}.sorted.bam")),
            bai=temp(pj(BAM,"{sample}.{readgroup}.sorted.bam.bai")),
            dragmap_log=pj(STAT,"{sample}.{readgroup}.dragmap.log"),
            stats=pj(STAT,"{sample}.{readgroup}.dechimer_stats.tsv"),
            badmap_fastq1=temp(pj(FQ_BADMAP,"{sample}.{readgroup}.badmap_R1.fastq.gz")),
            badmap_fastq2=temp(pj(FQ_BADMAP,"{sample}.{readgroup}.badmap_R2.fastq.gz")),
            merge_stats=ensure(
                pj(STAT,"{sample}.{readgroup}.merge_stats.tsv"),
                non_empty=True,
            ),
            checked=temp(pj(BAM,"{sample}.{readgroup}.bam_checked")),
            check_stats=pj(STAT,"{sample}.{readgroup}.bam_check_stats.tsv")
        log:
            runner=pj(LOG,"Aligner","{sample}.{readgroup}.align_fused.log"),
            io_profile=pj(
                LOG,"Aligner","{sample}.{readgroup}.align_fused.io.json"
            )
        params:
            runner=srcdir('scripts/run_fused_alignment.py'),
            ref_dir=get_refdir_by_validated_sex,
            bam_merge=srcdir(BAMMERGE),
            dechimer=srcdir(DECHIMER),
            bam_stats=srcdir('scripts/bam_stats_compare_hts.py'),
            lease_mode=ALIGNMENT_LEASE_MODE,
            lease_command=zslurm_lease_command(config),
            low_memory_mb=_fused_low_memory_mb,
            ignore_qual_flag=_fused_ignore_qual_flag
        conda: CONDA_ALIGN_FUSED
        priority: 16
        resources:
            time=get_time('align_reads_fused'),
            n="22.75",
            use_threads=24,
            mem_mb=lambda wildcards, attempt: (
                (attempt - 1) * 0.25 * 40000 + 40000
            ),
            ssd_use="required",
            # First estimate: the sort tail adds input, output, and spill data
            # to the earlier aligned/merged/dechimer peak. Calibrate from the
            # per-phase align_fused.io.json measurements.
            ssd_gb=lambda wildcards, input: ssd_gb_for_inputs(
                input.prepared_fastq,
                factor=4.5,
                overhead_gb=6,
                minimum_gb=24,
            )
        shell:
            """
            python {params.runner:q} \
                --prepared-fastq1 {input.prepared_fastq[0]:q} \
                --prepared-fastq2 {input.prepared_fastq[1]:q} \
                --source-fastq1 {input.source_fastq[0]:q} \
                --source-fastq2 {input.source_fastq[1]:q} \
                --fastq-stats {input.fastq_stats:q} \
                --reference-dir {params.ref_dir:q} \
                --sample {wildcards.sample:q} \
                --readgroup {wildcards.readgroup:q} \
                --output-bam {output.bam:q} \
                --output-bai {output.bai:q} \
                --dragmap-log {output.dragmap_log:q} \
                --dechimer-stats {output.stats:q} \
                --badmap-fastq1 {output.badmap_fastq1:q} \
                --badmap-fastq2 {output.badmap_fastq2:q} \
                --merge-stats {output.merge_stats:q} \
                --checked {output.checked:q} \
                --check-stats {output.check_stats:q} \
                --metrics {log.io_profile:q} \
                --bam-merge {params.bam_merge:q} \
                --dechimer {params.dechimer:q} \
                --bam-stats {params.bam_stats:q} \
                --dechimer-threshold {DECHIMER_THRESHOLD} \
                --align-threads {resources.use_threads} \
                --initial-cores {resources.n} \
                --initial-memory-mb {resources.mem_mb} \
                --low-cores 6 \
                --low-memory-mb {params.low_memory_mb} \
                --sort-threads 2 \
                --sort-memory-mb 6000 \
                --sort-compression-level 1 \
                --lease-mode {params.lease_mode:q} \
                --lease-command {params.lease_command:q} \
                --ssd-gb {resources.ssd_gb} \
                {params.ignore_qual_flag} \
                2> {log.runner:q}
            """

    # The fused rule deliberately retains the legacy provenance/stat outputs.
    # Prefer it for every overlapping product when both rule definitions are
    # present in the DAG.
    ruleorder: align_reads_fused > merge_bam_alignment_dechimer > align_reads

if not FUSE_ALIGNMENT_PHASES:
    rule sort_bam_alignment:
        """Sort bam alignment by chromosome and position."""
        input:
            in_bam=pj(BAM,"{sample}.{readgroup}.dechimer.bam")
        output:
            bam=temp(pj(BAM,"{sample}.{readgroup}.sorted.bam")),
            bai=temp(pj(BAM,"{sample}.{readgroup}.sorted.bam.bai"))
        conda: CONDA_MAIN
        log:
            samtools_sort=pj(LOG,"Aligner","{sample}.{readgroup}.samtools_sort.log"),
            io_profile=pj(LOG,"Aligner","{sample}.{readgroup}.sort_bam_alignment.io.json"),
        priority: 17
        resources:
            time = get_time('sort_bam_alignment'),
            tmpdir=tmpdir,
            n="1.3",
            mem_mb=13000,
            ssd_use="required",
            ssd_gb=lambda wildcards, input: ssd_gb_for_inputs(input.in_bam, factor=1.15, overhead_gb=2, minimum_gb=6)
        params:
            sort_runner=srcdir("scripts/run_samtools_sort_ssd.py"),
            memory_per_core=6000
        shell:
            """
                python {params.sort_runner:q} \
                    --input {input.in_bam:q} \
                    --output-bam {output.bam:q} \
                    --output-bai {output.bai:q} \
                    --metrics {log.io_profile:q} \
                    --threads 2 \
                    --memory-mb {params.memory_per_core} \
                    --compression-level 1 \
                    --ssd-gb {resources.ssd_gb} \
                    2> {log.samtools_sort:q}
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
        time = get_time('merge_rgs'),
        n="1",
        mem_mb=1250
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
        time = get_time('merge_rgs_badmap'),
        n="1",
        mem_mb=150
    shell:
        """
        zcat {input.fastq} | bgzip > {output.fastq} 
        """


def get_mem_mb_markdup(wildcards, attempt):  #{{{
    # Intentionally size for representative use rather than the long tail.
    # zslurm_chief keeps node-level memory headroom, while a rare failure can
    # use Snakemake's attempt-based escalation below.
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
        time = get_time('markdup'),
        n="1",
        mem_mb=get_mem_mb_markdup,
        # fastqs + intermediate bams are gone once markdup runs -> hand that share of
        # the start_sample reservation back now (see active_release_markdup).
        active_use_remove=active_release_markdup,
        temp_loc=lambda wildcards: pj(f"markdup_temporary_{wildcards.sample}"),
        ssd_use="required",
        ssd_gb=4
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
                if [ ! -d "$TMP_SSD" ] || [ ! -w "$TMP_SSD" ]; then CAND=$(ls -1dt /scratch-node/${{USER}}.* 2>/dev/null | head -n1 || true); if [ -n "${{CAND:-}}" ] && [ -d "$CAND" ] && [ -w "$CAND" ]; then TMP_SSD="$CAND"; fi; fi
                MDROOT=""
                if [ -d "$TMP_SSD" ] && [ -w "$TMP_SSD" ]; then TMPDIR_USE="$TMP_SSD"; elif [ -n "${{SLURM_TMPDIR:-}}" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then TMPDIR_USE="$SLURM_TMPDIR"; else TMPDIR_USE="{resources.temp_loc}"; MDROOT="{resources.temp_loc}"; fi
                JOB_ID="${{SLURM_JOB_ID}}"; if [ -z "$JOB_ID" ]; then JOB_ID="${{SLURM_JOBID}}"; fi; if [ -z "$JOB_ID" ]; then JOB_ID="$$"; fi
                MDTMP="$TMPDIR_USE/markdup/$JOB_ID/{wildcards.sample}"
                mkdir -p "$(dirname "$MDTMP")"
                # markdup had NO temp cleanup (unlike aligner_sort/aligner_fastq), so its
                # working-dir fallback left markdup_temporary_<sample> dirs behind. Clean
                # our own job subtree on exit; rmdir shared parents only if empty; and
                # remove the per-sample fallback root (MDROOT) -- never $TMPDIR_USE itself
                # when it is shared scratch (/scratch-node or $SLURM_TMPDIR).
                trap 'rm -rf "$TMPDIR_USE/markdup/$JOB_ID" 2>/dev/null || true; rmdir "$TMPDIR_USE/markdup" 2>/dev/null || true; [ -n "${{MDROOT:-}}" ] && rmdir "$MDROOT" 2>/dev/null || true' EXIT INT TERM
                samtools markdup -T "$MDTMP" -f {output.MD_stat} -S -d {params.machine} {input.bam} --write-index {output.mdbams}##idx##{output.mdbams_bai} 2> {log.samtools_markdup}
            fi
        """


localrules: release_materialized_source

rule release_materialized_source:
    """Keep external source data through alignment, then let temp GC reclaim it."""
    input:
        ready=pj(SOURCEDIR, "{sample}.route_ready"),
        materialized=lambda wildcards: external_data_dir(
            wildcards['sample'], SAMPLEINFO[wildcards['sample']]
        ),
        bam=pj(BAM, "{sample}.markdup.bam")
    output:
        marker=temp(touch(pj(SOURCEDIR, "{sample}.materialized_consumed")))
    wildcard_constraints:
        sample=START_SAMPLE_EXTERNAL_PATTERN
    shell:
        "touch {output.marker:q}"


rule mCRAM:
    """Convert bam to mapped cram."""
    input:
        bam=rules.markdup.output.mdbams,
        bai=rules.markdup.output.mdbams_bai
    output:
        cram=temp(pj(CRAM,"{sample}.mapped_hg38.cram")),
        crai=temp(pj(CRAM,"{sample}.mapped_hg38.cram.crai"))
    resources:
        time = get_time('mCRAM'),
        n="2",
        # Full-depth Knight WGS CRAM conversion was observed at about 1.44 GB,
        # leaving virtually no headroom with the previous 1.5 GB request.
        mem_mb=2500
    priority: 30
    conda: CONDA_MAIN
    log:
        pj(LOG,"Aligner","{sample}.mCRAM.log")
    shell:
        "samtools view --output-fmt cram,version=3.1,archive --reference {REF} -@ {resources.n} --write-index -o {output.cram}##idx##{output.crai} {input.bam} 2> {log}"
