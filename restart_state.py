"""Reuse completed samples without reconstructing their processing DAG.

The root invocation freezes the reuse decision in a shared manifest. Its path
is inherited by nested worker invocations, so completing a sample during a run
cannot change which rules a later worker sees.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import tempfile
import uuid

MANIFEST_ENV = 'SHORT_READ_RESTART_MANIFEST'
SCHEMA = 1


def sample_pattern(samples):
    names = sorted(re.escape(s) for s in samples)
    return '(?:' + '|'.join(names) + ')' if names else r'(?!)'


def parse_bool(value):
    if isinstance(value, bool):
        return value
    if str(value).lower() in ('1', 'true', 'yes', 'on'):
        return True
    if str(value).lower() in ('0', 'false', 'no', 'off'):
        return False
    raise ValueError(f'Expected boolean, got {value!r}')


def sample_digest(samples):
    return hashlib.sha256('\n'.join(sorted(samples)).encode()).hexdigest()


def product_settings(config):
    return {key:config.get(key, default) for key, default in
            (('END_POINT', 'gVCF'), ('caller', 'Deepvariant'), ('chrM', 'Yes'))}


def load_manifest(path, root, samples):
    data = json.loads(Path(path).read_text())
    if data.get('schema_version') != SCHEMA:
        raise ValueError('Unsupported completed-sample manifest schema')
    if data.get('workdir') != str(Path(root).resolve()):
        raise ValueError('Completed-sample manifest belongs to another workdir')
    if data.get('sample_digest') != sample_digest(samples):
        raise ValueError('Sample list changed; prepare a new restart manifest')
    reused = frozenset(data['reused_samples'])
    if not reused <= set(samples):
        raise ValueError('Completed-sample manifest contains unknown samples')
    rebuild = set(data.get('rebuild_samples', []))
    if rebuild - set(samples) or rebuild & reused:
        raise ValueError('Completed-sample manifest contains an invalid rebuild selection')
    return data


def allowed_stat_member(name, sample):
    path = PurePosixPath(name)
    return (bool(path.parts) and not path.is_absolute() and '..' not in path.parts
            and path.parts[0] in ('stats', 'kmer')
            and path.name.startswith(sample + '.'))


def restore_stats(root, sample, replace_newer=False, backup_root=None):
    """Restore archived results atomically, retaining the original timestamps.

    An interrupted redundant rerun can overwrite old stats. Explicit recovery
    may restore those newer files too, retaining a backup before replacement.
    Normal startup only restores missing files.
    """
    root = Path(root)
    archive = root / 'stats' / f'{sample}.stats.tar.gz'
    finished_time = (root / 'source' / f'{sample}.finished').stat().st_mtime
    restored = []
    with tarfile.open(archive, 'r:gz') as handle:
        for member in handle:
            if not member.isfile() or not allowed_stat_member(member.name, sample):
                raise ValueError(f'Unexpected stats archive member: {member.name!r}')
            dest = root / member.name
            if not dest.resolve().is_relative_to(root.resolve()):
                raise ValueError(f'Archive destination escapes workdir: {dest}')
            exists = dest.exists()
            if exists and not (replace_newer and dest.stat().st_mtime > finished_time):
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            if exists:
                if backup_root is None:
                    raise ValueError('Replacing stats requires a backup directory')
                backup = Path(backup_root) / member.name
                backup.parent.mkdir(parents=True, exist_ok=True)
                if backup.exists():
                    raise ValueError(f'Recovery backup already exists: {backup}')
                shutil.copy2(dest, backup)
            fd, tmp = tempfile.mkstemp(prefix='.' + dest.name + '.', dir=dest.parent)
            try:
                with os.fdopen(fd, 'wb') as output, handle.extractfile(member) as source:
                    shutil.copyfileobj(source, output)
                os.chmod(tmp, member.mode & 0o777)
                os.utime(tmp, (member.mtime, member.mtime))
                os.replace(tmp, dest)
            finally:
                if os.path.exists(tmp):
                    os.unlink(tmp)
            restored.append(member.name)
    return restored


def completion_products(sample, sinfo, config, levels):
    """Files needed by the supported gVCF endpoint's cohort consumers."""
    products = [f'source/{sample}.finished', f'cram/{sample}.mapped_hg38.cram.copied',
                f'stats/{sample}.stats.tar.gz', f'sampleinfo/{sample}.dat',
                f'kmer/{sample}.result.yaml']
    for suffix in ('samtools.stat', 'samtools.exome.stat', 'bam_all.tsv', 'bam_exome.tsv',
                   'pre_adapter_summary_metrics', 'bait_bias_summary_metrics',
                   'pre_adapter_detail_metrics', 'bait_bias_detail_metrics',
                   'hs_metrics', 'phase_stats.tsv'):
        products.append(f'stats/{sample}.{suffix}')
    products.extend([f'stats/contam/{sample}.verifybamid.pca2.selfSM',
                     f'stats/cov/{sample}.regions.bed.gz'])
    for suffix in ('report.tsv', 'bracken_report.tsv', 'kraken_summary.tsv',
                   'report_bracken_species.tsv', 'read_classification.tsv.gz'):
        products.append(f'kraken/{sample}.{suffix}')
    if config.get('chrM', 'Yes') == 'Yes':
        products += [f'chrM_analysis/variants/gvcf/{sample}.chrM_merged_BP_annotated.g.vcf.gz',
                     f'chrM_analysis/variants/gvcf/{sample}.chrM_merged_BP_annotated.g.vcf.gz.tbi',
                     f'chrM_analysis/{sample}.done',
                     f'stats/{sample}.chrM_read_stats.tsv', f'stats/{sample}.numt_read_stats.tsv']
        # NUMT gVCFs only feed the per-sample .done marker, not cohort uploads.
        # Older runs legitimately collected those temp files after completion.
    caller = config.get('caller', 'Deepvariant')
    if caller in ('Deepvariant', 'BOTH'):
        wgs = 'wgs' in sinfo['sample_type'].lower()
        regions = levels[1] if wgs else levels[0]
        for region in regions:
            paths = [f'deepvariant/gVCF/{region}/{sample}.{region}.wg.vcf.gz']
            if wgs:
                paths.append(f'deepvariant/gVCF/exome_extract/{region}/{sample}.{region}.wg.vcf.gz')
            for path in paths:
                products.extend([path, path + '.tbi'])
            products.append(f'stats/deepvariant_bcftools/{sample}.{region}.summary.tsv')
    if caller in ('HaplotypeCaller', 'BOTH'):
        wgs = 'wgs' in sinfo['sample_type'].lower()
        for region in levels[1] if wgs else levels[0]:
            folder = 'exome_extract' if wgs else 'reblock'
            products.append(f'gvcf_conv/{folder}/{region}/{sample}.{region}.wg.vcf.gz')
    return products


def refresh_retrieval(samples, batches, reused, rebuild):
    """Filter stage contents, without changing stable batch IDs or membership."""
    for sample in reused:
        samples[sample]['need_retrieval'] = False
    for sample in rebuild:
        if samples[sample].get('from_external'):
            samples[sample]['need_retrieval'] = True
    for routes in batches.values():
        for route_batches in routes.values():
            for batch in route_batches:
                batch['size'] = sum(samples[s]['filesize'] for s in batch['samples']
                                    if s in samples and samples[s].get('need_retrieval'))


def verify_rule_filters(rules, reused):
    """Reject per-rule overrides that reopen completed samples before any DAG work."""
    if not reused:
        return
    from snakemake.io import WILDCARD_REGEX
    allowed = {'tar_badmap_fastqs', 'copy_badmap_to_dcache'}
    checked = 0
    for rule in rules:
        if rule.name in allowed:
            continue
        constraints = {m.group('constraint') for output in rule.output
                       for m in WILDCARD_REGEX.finditer(str(output))
                       if m.group('name') == 'sample'}
        for constraint in constraints:
            pattern = re.compile(constraint or '.+')
            example = next((s for s in reused if pattern.fullmatch(s)), None)
            if example is not None:
                raise ValueError(f'Rule {rule.name} still accepts reused sample {example}; '
                                 'its sample constraint must also exclude REUSED_SAMPLES')
        checked += bool(constraints)
    print(f'[restart] Verified finished-sample exclusion in {checked} processing rules', flush=True)


def prepare(root, samples, config, levels, *, replace_newer=False, backup_root=None):
    root = Path(root).resolve()
    rebuild = config.get('restart_rebuild_samples', [])
    if isinstance(rebuild, str):
        rebuild = [s for s in rebuild.split(',') if s]
    rebuild = set(rebuild)
    if rebuild - set(samples):
        raise ValueError(f'Unknown rebuild samples: {sorted(rebuild - set(samples))}')
    finished = {p.name[:-9] for p in (root/'source').glob('*.finished')} & set(samples)
    reused = finished - rebuild
    restored = {}
    # Directory inventories avoid hundreds of thousands of serial GPFS stats.
    inventories = {}
    def exists(relative):
        parent, name = os.path.split(relative)
        if parent not in inventories:
            directory = root/parent
            inventories[parent] = set(os.listdir(directory)) if directory.is_dir() else set()
        return name in inventories[parent]

    # Recover only when needed. The full DAG validates other missing inputs;
    # expensive producers for reused samples are never available as a fallback.
    recovery = [s for s in sorted(reused) if replace_newer
                or not exists(f'kmer/{s}.result.yaml')
                or any(not exists(p) for p in completion_products(s, samples[s], config, levels)
                       if p.startswith('stats/'))]
    if recovery:
        print(f'[restart] Recovering archived statistics for {len(recovery)} completed samples', flush=True)
        with ThreadPoolExecutor(max_workers=4) as pool:
            tasks = {pool.submit(restore_stats, root, s, replace_newer, backup_root):s for s in recovery}
            for number, future in enumerate(as_completed(tasks), 1):
                restored[tasks[future]] = future.result()
                if number % 100 == 0 or number == len(tasks):
                    print(f'[restart] Recovered {number}/{len(tasks)} stats archives', flush=True)
    inventories.clear()  # Recovery may have published previously missing files.
    missing = {}
    for s in sorted(reused):
        paths = completion_products(s, samples[s], config, levels)
        absent = [p for p in paths if not exists(p)]
        if not exists(f'fq_badmap/{s}.badmap.tar.copied'):
            if not exists(f'fq_badmap/{s}.badmap.fastqs.tar.gz'):
                if not all(exists(f'fq_badmap/{s}.badmap.{r}.fastq.gz') for r in ('R1','R2')):
                    absent.append(f'fq_badmap/{s}.badmap.tar.copied (or uploadable badmap data)')
        if absent:
            missing[s] = absent
    if missing:
        error_path = root/'.snakemake/restart_missing_products.json'
        error_path.parent.mkdir(exist_ok=True)
        error_path.write_text(json.dumps(missing, indent=2))
        raise ValueError(f'{len(missing)} finished samples lack required products; see {error_path}. '
                         'Recover them or explicitly select restart_rebuild_samples. No samples were silently reprocessed.')
    data = {'schema_version':SCHEMA, 'workdir':str(root), 'sample_digest':sample_digest(samples),
            'product_settings':product_settings(config),
            'reused_samples':sorted(reused), 'rebuild_samples':sorted(rebuild),
            'restored_files':restored}
    directory = root/'.snakemake/restart_manifests'
    directory.mkdir(parents=True, exist_ok=True)
    path = directory/f'{uuid.uuid4().hex}.json'
    path.write_text(json.dumps(data, indent=2))
    error_path = root/'.snakemake/restart_missing_products.json'
    if error_path.exists():
        error_path.replace(directory/f'{path.stem}.resolved_missing_products.json')
    return path, data


def configure(config, samples, levels, *, is_main_process=True):
    """Freeze on the controller, load on workers, never rescan on a worker.

    A worker submitted by a pre-upgrade controller has no manifest. Keep its
    original (unfiltered) DAG semantics during a rolling source deployment.
    New controllers export the manifest, including when the reuse set is empty.
    An explicitly supplied but invalid manifest still fails validation.
    """
    enabled = parse_bool(config.get('reuse_finished_samples', True))
    if not enabled or config.get('END_POINT', 'gVCF') != 'gVCF':
        return frozenset(), frozenset(), r'[\w\d_\-@]+'
    path = config.get('restart_manifest') or os.environ.get(MANIFEST_ENV)
    if not is_main_process and not path:
        print('[restart] Legacy worker without a frozen manifest; '
              'keeping controller selection (no completed-sample scan)', flush=True)
        return frozenset(), frozenset(), r'[\w\d_\-@]+'
    if path:
        data = load_manifest(path, Path.cwd(), samples)
        if data.get('product_settings', product_settings({})) != product_settings(config):
            raise ValueError('Endpoint/caller/chrM settings changed; prepare a new restart manifest')
        explicit = config.get('restart_rebuild_samples')
        if explicit:
            names = explicit.split(',') if isinstance(explicit, str) else explicit
            if set(names) != set(data.get('rebuild_samples', [])):
                raise ValueError('Rebuild selection differs from the frozen restart manifest')
    else:
        path, data = prepare(Path.cwd(), samples, config, levels)
    os.environ[MANIFEST_ENV] = str(Path(path).resolve())
    reused = frozenset(data['reused_samples'])
    rebuild = frozenset(data.get('rebuild_samples', []))
    print(f'[restart] Reusing {len(reused)} completed samples; '
          f'{len(samples)-len(reused)} samples eligible for processing; manifest={path}', flush=True)
    pattern = sample_pattern(set(samples)-reused) if reused else r'[\w\d_\-@]+'
    return reused, rebuild, pattern
