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
import readgroup_checkpoint as _readgroup_checkpoint

import read_samples
from common import *
import utils

FAILED_JOBS_LOG = pj(LOG, "failed_jobs.jsonl")


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
    sample=PROCESSING_SAMPLE_PATTERN,
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
            filename = _readgroup_checkpoint.output(sample)
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

        if any(f.startswith(route + ':') for route in ('archive', 'dcache', 's3')):
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
        # CRAM creation, encryption and upload are fused in Encrypt.smk. The
        # durable receipt replaces the former plaintext-CRAM aggregate target,
        # so Align cannot pin CRAM payloads on active GPFS.
        expand("{cram}/{sample}.mapped_hg38.cram.copied", sample=sample_names, cram=CRAM)

checkpoint get_readgroups:
    """Get the readgroup info for a sample.

    Once the checkpoint rule is executed, the readgroup info is available in the file SAMPLEINFODIR/<sample>.dat.

    Once readgroup info is available, snakemake will recalculate the DAG.
    """
    input:
        get_source_files,
        ancient(pj(SOURCEDIR,"{sample}.started"))
    output:
        pj(SAMPLEINFODIR,"{sample}.dat")
    resources:
        time = get_time('get_readgroups'),
        n="0.5",
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

# Snakemake isolates the checkpoint namespace of every workflow module.  Make
# this canonical checkpoint available to the Stat and Kraken modules after it
# has been registered by the checkpoint directive above.
_readgroup_checkpoint.register(checkpoints.get_readgroups)

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
                    res = subprocess.run([DAGET, "-av", *chunk], capture_output=True, text=True)
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
                        res = subprocess.run([DALS, "-l", *chunk], capture_output=True, text=True)
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


def _dcache_stage_timeout_seconds(resource_time_seconds):
    """Keep tape-stage polling inside this attempt's ZSlurm time budget."""

    configured_timeout = config.get('dcache_stage_timeout')
    if configured_timeout is not None:
        return int(configured_timeout)
    return max(60, int(resource_time_seconds) - 600)


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
                stage_timeout_seconds = _dcache_stage_timeout_seconds(
                    resources.time
                )
                print(
                    f"[dcache_get] Stage timeout {stage_timeout_seconds}s "
                    f"within {int(resources.time)}s job budget",
                    flush=True,
                )
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
                        str(stage_timeout_seconds),
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

    if wildcards['sample'] in REUSED_SAMPLES:
        # Keep stable batch numbering, but no graph edge from an adopted
        # completed sample back to a shared (possibly missing) stage receipt.
        return []
    sample = SAMPLEINFO[wildcards['sample']]
    route = _start_sample_route(wildcards)
    if route in ('active', 's3'):
        # S3 is object storage and needs no batch-level tape staging.  The
        # per-sample start rule performs the complete requester-pays download.
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
        's3': 's3_to_active',
    }[_start_sample_route(wildcards)]
    return get_time(rule_name)(wildcards, attempt)


def _start_sample_partition(wildcards):
    return 'archive' if _start_sample_route(wildcards) == 'archive' else 'compute'


def _start_sample_cores(wildcards):
    # Snakemake rounds float resources to integers during expansion. Keep
    # fractional ZSlurm core reservations as strings so the native executor
    # receives their exact values.
    return {'active': '0.1', 'archive': '0.6', 'dcache': '1.0', 's3': '0.5'}[
        _start_sample_route(wildcards)
    ]


def _start_sample_mem_mb(wildcards):
    return {'dcache': 1024, 's3': 512}.get(
        _start_sample_route(wildcards), 256
    )


def _start_sample_active_add(wildcards):
    if wildcards['sample'] in REBUILD_SAMPLES and os.path.exists(
        pj(SOURCEDIR, wildcards['sample'] + '.finished')
    ):
        # Only a completed lifecycle needs a fresh reservation on rebuild.
        # Selecting an unfinished sample must not double its existing budget.
        return active_use_gb(wildcards)
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
        if sample not in REUSED_SAMPLES
        and _start_sample_route({'sample': sample}) in routes
    )
    return '(?:' + '|'.join(samples) + ')' if samples else r'(?!)'


START_SAMPLE_ACTIVE_PATTERN = _start_sample_pattern('active')
START_SAMPLE_ARCHIVE_PATTERN = _start_sample_pattern('archive')
START_SAMPLE_DCACHE_PATTERN = _start_sample_pattern('dcache')
START_SAMPLE_S3_PATTERN = _start_sample_pattern('s3')
START_SAMPLE_EXTERNAL_PATTERN = _start_sample_pattern('archive', 'dcache', 's3')


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
        aws_cli=str(config.get('aws_cli', 'aws')),
        s3_max_attempts=max(1, int(config.get('s3_max_attempts', 6))),
        s3_initial_backoff_seconds=max(
            0, float(config.get('s3_initial_backoff_seconds', 30))
        ),
        s3_max_backoff_seconds=max(
            0, float(config.get('s3_max_backoff_seconds', 300))
        ),
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


rule start_sample_s3:
    """Download one requester-pays S3 sample into tracked active storage."""
    input:
        retrieve_batch
    output:
        started=pj(SOURCEDIR,"{sample}.started"),
        route_ready=temp(pj(SOURCEDIR,"{sample}.route_ready")),
        materialized=temp(directory(pj(SOURCEDIR,"{sample}.s3_data")))
    wildcard_constraints:
        sample=START_SAMPLE_S3_PATTERN
    resources:
        time=_start_sample_time,
        partition=_start_sample_partition,
        active_use_add=_start_sample_active_add,
        arch_use_remove=0,
        dcache_use_remove=0,
        # S3 transfers use their own manager-wide concurrency pool and do not
        # consume dCache download capacity.
        s3_download_slots=1,
        mem_mb=_start_sample_mem_mb,
        n=_start_sample_cores
    params:
        expected_route='s3',
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
        res = 7500
    return attempt * res


#}}}

def get_n_split_alignments(wildcards):  #{{{
    """Reserve measured split throughput only when the rule does real work."""
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    readgroups = [
        rg for rg in sinfo['readgroups']
        if wildcards['filename'] in rg['file']
    ]
    if len(readgroups) > 1 or sinfo.get('erf_correct', False):
        # A real-data sweep measured 2.87 cores on average with samtools -@ 3.
        return "3.0"
    # The single-RG path only creates a link and completion marker.
    return "0.9"


#}}}

def get_threads_split_alignments(wildcards):  #{{{
    """Keep workers submitted by an older chief within their CPU lease."""
    lease_max_cores = os.environ.get('ZSLURM_LEASE_MAX_CORES')
    if lease_max_cores is not None and float(lease_max_cores) < 2.0:
        # Before this change split reserved 0.9 core and effectively passed
        # zero extra HTSlib threads via atoi("0.9").  Such a worker rereads
        # this file, so preserve that behavior until its old chief is gone.
        return 0
    return 3


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
        n=get_n_split_alignments,
        use_threads=get_threads_split_alignments,
        mem_mb=get_mem_mb_split_alignments
    conda: CONDA_MAIN
    priority: 99
    params:
        fixer=srcdir('scripts/fix_bam_rg_pairs'),
        reference_selector=srcdir('scripts/select_cram_reference.py'),
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
            selected_reference = shell(
                "python {selector:q} --cram {cram:q} "
                "--primary-reference {primary:q} "
                "--hg19-reference {hg19:q} "
                "--hg19-b37-chry-reference {hg19_b37_chry:q} "
                "--hg38-reference {hg38:q}",
                selector=str(params.reference_selector),
                cram=readfile,
                primary=reference_file,
                hg19=HG19_REFERENCE,
                hg19_b37_chry=HG19_B37_CHRY_REFERENCE,
                hg38=HG38_CRAM_REFERENCE,
                read=True,
            ).strip()
            cramref = f"--reference {shlex.quote(selected_reference)}"
            rflag = f"-r {shlex.quote(selected_reference)}"
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
            cramref = ""
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
                    {params.fixer} -i {readfile} -o {output.readgroups}/{wildcards.sample}.{readgroup_id}.{extension_in} {rflag} --threads {resources.use_threads}
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
                pre = f"{params.fixer} -i {readfile} -o {sanitized} {rflag} --threads {resources.use_threads}\n                "
                inpath = sanitized
            else:
                pre = ""
                inpath = readfile
            cmd = f"""
                set -euo pipefail
                mkdir -p {output.readgroups}
                {pre}samtools split -@ {resources.use_threads} --output-fmt {output_fmt} {cramref} {inpath} -f "{output.readgroups}/{wildcards.sample}.%!.{extension}"
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

NATIVE_FASTQ_SAMPLE_PATTERN = '(?:' + '|'.join(
    re.escape(sample)
    for sample, sinfo in SAMPLEINFO.items()
    if sample not in REUSED_SAMPLES
    and sinfo.get('file_type') in {'fastq', 'fastq_paired'}
) + ')'
if NATIVE_FASTQ_SAMPLE_PATTERN == '(?:)':
    NATIVE_FASTQ_SAMPLE_PATTERN = r'(?!)'


rule adapter_removal:
    """Native FASTQ entry point to the shared adapter-processing implementation."""
    input:
        get_fastqpaired,
        ancient(pj(SOURCEDIR,"{sample}.started"))
    output:
        for_f=temp(pj(FQ,"{sample}.{readgroup}.fastq.cut_1.fq.gz")),
        rev_f=temp(pj(FQ,"{sample}.{readgroup}.fastq.cut_2.fq.gz")),
        adapter_removal=ensure(pj(STAT,"{sample}.{readgroup}.adapter_removal.log"), non_empty=True),
        fastq_stats=pj(STAT,"{sample}.{readgroup}.fastq.stats.tsv"),
        adapters=pj(STAT,"{sample}.{readgroup}.fastq.adapters"),
    wildcard_constraints:
        sample=NATIVE_FASTQ_SAMPLE_PATTERN
    log:
        runner=pj(LOG,"Aligner","{sample}.{readgroup}.adapter_removal.log"),
        io_profile=pj(LOG,"Aligner","{sample}.{readgroup}.adapter_removal.io.json")
    priority: 10
    conda: CONDA_MAIN
    params:
        runner=srcdir('scripts/run_fastq_adapter.py'),
        adapters=ADAPTERS,
        fastq_stats=srcdir('scripts/fastq_stats.py'),
        rmdups=srcdir('scripts/remove_interleaved_duplicates.py'),
        rescuer=srcdir('scripts/fastq_pair_rescue.py'),
        remove_duplicated_reads=lambda wc: int(SAMPLEINFO[wc.sample].get('remove_duplicated_reads', False)),
        error_file=lambda wc: str(SAMPLEINFO[wc.sample]['samplefile']) + '.errors'
    resources:
        time=get_time('adapter_removal'),
        n="5",
        mem_mb=512,
        attempt=lambda wildcards, attempt: attempt
    shell:
        """
        python {params.runner:q} \
            --input-forward {input[0]:q} --input-reverse {input[1]:q} \
            --sample {wildcards.sample:q} --readgroup {wildcards.readgroup:q} \
            --adapter-list {params.adapters:q} \
            --fastq-stats-script {params.fastq_stats:q} \
            --remove-duplicates-script {params.rmdups:q} \
            --pair-rescue-script {params.rescuer:q} \
            --remove-duplicated-reads {params.remove_duplicated_reads} \
            --attempt {resources.attempt} --error-file {params.error_file:q} \
            --output-forward {output.for_f:q} --output-reverse {output.rev_f:q} \
            --output-adapter-log {output.adapter_removal:q} \
            --output-fastq-stats {output.fastq_stats:q} \
            --output-adapters {output.adapters:q} \
            --adapter-memory-mb {resources.mem_mb} \
            --metrics {log.io_profile:q} 2> {log.runner:q}
        """


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
    if sample not in REUSED_SAMPLES
    and sinfo.get('file_type') not in {'fastq', 'fastq_paired'}
) + ')'
if EXTERNAL_ALIGNMENT_SAMPLE_PATTERN == '(?:)':
    EXTERNAL_ALIGNMENT_SAMPLE_PATTERN = r'(?!)'


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
        hg19_reference=HG19_REFERENCE,
        hg19_b37_chry_reference=HG19_B37_CHRY_REFERENCE,
        hg38_reference=HG38_CRAM_REFERENCE,
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
        # AdapterRemoval uses all five cores in production. Keep them for
        # the whole fused job so no phase has to reacquire released cores.
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
            --hg19-reference {params.hg19_reference:q} \
            --hg19-b37-chry-reference {params.hg19_b37_chry_reference:q} \
            --hg38-reference {params.hg38_reference:q} \
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
            --adapter-memory-mb 768 \
            --lease-mode {params.lease_mode:q} \
            --lease-command {params.lease_command:q} \
            --ssd-gb {resources.ssd_gb} \
            2> {log.runner:q}
        """


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


# KMC stage 1 can create 12 FASTQ readers, 12 splitters, one bin writer and
# one bin reader with the command-line settings below.  Reserve for observed
# average CPU use rather than the theoretical thread maximum: 2 cores yielded
# about 67% CPU reservation efficiency over 2,778 successful production jobs.
KMC_RESERVED_CORES = 2
KMER_SEX_LEASE_MODE = str(
    config.get('kmer_sex_lease_mode', 'required')
).strip().lower()
if KMER_SEX_LEASE_MODE not in {'required', 'optional', 'disabled'}:
    raise ValueError(
        "kmer_sex_lease_mode must be required, optional, or disabled"
    )

rule kmer_sex_fused:
    """Build the KMC database and validate sex without GPFS intermediates."""
    input:
        fastq=get_all_prepared_fastq
    output:
        yaml=pj(KMER,"{sample}.result.yaml"),
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
        # This is an average scheduler reservation; the patched KMC
        # invocation itself is unchanged.
        n="1.6",
        use_threads=KMC_RESERVED_CORES,
        # Estimated standalone peak RSS is about 33.8 GiB. Reserve ten
        # percent less than the former 38-GB baseline; retries retain the
        # existing 50% step-up.
        mem_mb=lambda wildcards, attempt: (
            (attempt - 1) * 0.5 * 34200 + 34200
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
            --kmc-threads {resources.use_threads} \
            --low-cores 1.0 \
            --low-memory-mb 3000 \
            --lease-mode {params.lease_mode:q} \
            --lease-command {params.lease_command:q} \
            --ssd-gb {resources.ssd_gb} \
            2> {log.runner:q}
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


ALIGNMENT_LEASE_MODE = str(
    config.get('alignment_lease_mode', 'required')
).strip().lower()
if ALIGNMENT_LEASE_MODE not in {'required', 'optional', 'disabled'}:
    raise ValueError(
        "alignment_lease_mode must be required, optional, or disabled"
    )


def _fused_low_memory_mb(wildcards):
    # Phase-level PSS/RSS reconstruction puts the alignment-tail P95 near
    # 13.5 GB. Reserve ten percent below the former 15-GB target.
    return 13500


def _fused_ignore_qual_flag(wildcards):
    return (
        '--ignore-qual-checksum-diff'
        if bool(SAMPLEINFO[wildcards['sample']].get('erf_correct', False))
        else ''
    )

# --enable-sampling true used for (unmapped) bam input. It prevents bugs when in output bam information about whicj read is 1st or 2nd in pair.
#--preserve-map-align-order 1 was tested, so that unaligned and aligned bam have sam read order (requires thread synchronization). But reduces performance by 1/3.  Better to let mergebam job deal with the issue.


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
        # Initial alignment peak RSS is about 37 GiB. Start five percent
        # below the former 40-GB baseline while retaining retry scaling.
        mem_mb=lambda wildcards, attempt: (
            (attempt - 1) * 0.25 * 38000 + 38000
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
            --low-cores 2 \
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

# The fused rule retains the provenance/stat outputs consumed downstream.

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
        n="0.7",
        mem_mb=200
    shell:
        """
        zcat {input.fastq} | bgzip > {output.fastq} 
        """


def get_mem_mb_markdup(wildcards, attempt):  #{{{
    # Intentionally size for representative use rather than the long tail.
    # zslurm_chief keeps node-level memory headroom, while a rare failure can
    # use Snakemake's attempt-based escalation below.
    res = 3000 if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else 150
    #large range of memory usage for markdup
    return (attempt - 1) * res * 3 + res


#}}}

def get_n_merge_markdup(wildcards):
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    # Preserve the measured reservations of the former separate phases: the
    # merge phase averages 2.20 cores, while single-RG markdup averages 0.93.
    return "2.3" if len(sinfo['readgroups']) > 1 else "0.95"


def get_mem_mb_merge_markdup(wildcards, attempt):
    markdup_mem = get_mem_mb_markdup(wildcards, attempt)
    sinfo = sampleinfo(SAMPLEINFO, wildcards['sample'], checkpoint=True)
    # A multi-RG exome still has to accommodate the former 384-MB merge phase.
    return max(markdup_mem, 384 if len(sinfo['readgroups']) > 1 else 0)


def get_ssd_gb_merge_markdup(wildcards, input):
    # Multi-RG peak: local merged BAM + markdup spill + final BAM/index.
    # Single-RG jobs omit the local merged copy. These are deliberately
    # conservative first estimates and are recorded by the runner for tuning.
    factor = 4.0 if len(input.bam) > 1 else 3.0
    return ssd_gb_for_inputs(
        input.bam, factor=factor, overhead_gb=6, minimum_gb=12
    )


MERGE_MARKDUP_LEASE_MODE = str(
    config.get('merge_markdup_lease_mode', 'required')
).strip().lower()
if MERGE_MARKDUP_LEASE_MODE not in {'required', 'optional', 'disabled'}:
    raise ValueError(
        "merge_markdup_lease_mode must be required, optional, or disabled"
    )

rule markdup:
    """Merge read groups and mark duplicates without a GPFS merged BAM."""
    input:
        bam=get_readgroups_bam,
        bai=get_readgroups_bai,
        checks=get_readgroup_checks
    output:
        mdbams=temp(pj(BAM,"{sample}.markdup.bam")),
        mdbams_bai=temp(pj(BAM,"{sample}.markdup.bam.bai")),
        MD_stat=pj(STAT,"{sample}.markdup.stat")
    priority: 20
    params:
        runner=srcdir('scripts/run_fused_merge_markdup.py'),
        machine=2500,
        # machine error rate, default is 2500; NovaSeq uses 100
        no_dedup=lambda wildcards: 1 if SAMPLEINFO[wildcards['sample']]['no_dedup'] else 0,
        lease_mode=MERGE_MARKDUP_LEASE_MODE,
        lease_command=zslurm_lease_command(config)
    log:
        runner=pj(LOG,"Aligner","{sample}.merge_markdup_fused.log"),
        merge=pj(LOG,"Aligner","{sample}.mergereadgroups.log"),
        samtools_markdup=pj(LOG,"Aligner","{sample}.markdup.log"),
        io_profile=pj(LOG,"Aligner","{sample}.merge_markdup_fused.io.json")
    resources:
        time=get_time('merge_markdup_fused'),
        n=get_n_merge_markdup,
        use_threads=3,
        mem_mb=get_mem_mb_merge_markdup,
        # Do not require a late storage increase: it can block completed
        # readgroup inputs from reaching the phase that releases their budget.
        active_use_remove=active_release_markdup,
        ssd_use="required",
        ssd_gb=get_ssd_gb_merge_markdup
    conda: CONDA_MAIN
    shell:
        """
        python {params.runner:q} \
            --input-bam {input.bam:q} \
            --input-bai {input.bai:q} \
            --check-marker {input.checks:q} \
            --sample {wildcards.sample:q} \
            --output-bam {output.mdbams:q} \
            --output-bai {output.mdbams_bai:q} \
            --output-stat {output.MD_stat:q} \
            --merge-log {log.merge:q} \
            --markdup-log {log.samtools_markdup:q} \
            --metrics {log.io_profile:q} \
            --no-dedup {params.no_dedup} \
            --optical-distance {params.machine} \
            --merge-threads {resources.use_threads} \
            --initial-cores {resources.n} \
            --initial-memory-mb {resources.mem_mb} \
            --markdup-cores 0.95 \
            --markdup-memory-mb {resources.mem_mb} \
            --lease-mode {params.lease_mode:q} \
            --lease-command {params.lease_command:q} \
            --ssd-gb {resources.ssd_gb} \
            2> {log.runner:q}
        """


localrules: release_materialized_source

rule release_materialized_source:
    """Reclaim external source data once every readgroup alignment is valid."""
    input:
        ready=pj(SOURCEDIR, "{sample}.route_ready"),
        materialized=lambda wildcards: external_data_dir(
            wildcards['sample'], SAMPLEINFO[wildcards['sample']]
        ),
        bam=get_readgroups_bam,
        bai=get_readgroups_bai,
        checks=get_readgroup_checks
    output:
        marker=touch(pj(SOURCEDIR, "{sample}.materialized_consumed"))
    wildcard_constraints:
        sample=START_SAMPLE_EXTERNAL_PATTERN
    shell:
        "touch {output.marker:q}"
