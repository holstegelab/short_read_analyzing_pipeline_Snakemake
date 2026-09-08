"""Restart and rolling-worker checks using the real root workflow in scratch.

Never point these tests at a production run directory.
"""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import restart_state as rs

pytestmark = pytest.mark.skipif(
    not Path('/gpfs/work3/0/qtholstg/hg38_res_v2/databases/Adapters_illumina.txt').is_file(),
    reason='Real workflow smoke tests require the configured reference bundle',
)


def fixture(root):
    for mate in (1, 2):
        with gzip.open(root / f'R{mate}.fq.gz', 'wt') as handle:
            handle.write(f'@pair/{mate}\n' + 'ACGT' * 15 + '\n+\n!' + 'I' * 59 + '\n')
    (root / 'cohort.tsv').write_text(''.join(
        f'TEST\t{s}\tfastq_paired\tillumina_wgs\tWGS\tM\tR1.fq.gz\tR2.fq.gz\n'
        for s in ('TEST_A', 'TEST_B')
    ))
    (root / 'source').mkdir()
    # Deliberately incomplete products: a legacy worker must not run preflight.
    (root / 'source/TEST_A.started').touch()
    (root / 'source/TEST_A.finished').touch()


def frozen(root, caller='Deepvariant', chrm='Yes'):
    path = root / 'frozen.json'
    path.write_text(json.dumps({
        'schema_version': 1, 'workdir': str(root.resolve()),
        'sample_digest': rs.sample_digest(['TEST_A', 'TEST_B']),
        'product_settings': rs.product_settings({'caller': caller, 'chrM': chrm}),
        'reused_samples': ['TEST_B'], 'rebuild_samples': [],
    }))
    return path


def command(root, *args, manifest=None):
    env = dict(os.environ, PYTHONPATH=str(REPO), SKIP_FASTQ_VALIDATION='1')
    env.pop('SNAKEMAKE_PROFILE', None)
    env.pop(rs.MANIFEST_ENV, None)
    if manifest:
        env[rs.MANIFEST_ENV] = str(manifest)
    result = subprocess.run(
        [sys.executable, '-m', 'snakemake', '--snakefile', str(REPO / 'Snakefile'),
         '--cores', '1', '--nolock', *args],
        cwd=root, env=env, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, timeout=90,
    )
    assert result.returncode == 0, result.stdout
    return result.stdout


@pytest.mark.parametrize('mode', ['subprocess', 'remote'])
def test_pre_upgrade_worker_parses_without_adopting_newly_finished_sample(tmp_path, mode):
    fixture(tmp_path)
    output = command(tmp_path, '--list-rules', '--mode', mode)
    assert 'Legacy worker without a frozen manifest' in output
    assert 'external_adapter_fused' in output
    assert not (tmp_path / '.snakemake/restart_manifests').exists()
    assert not (tmp_path / '.snakemake/restart_missing_products.json').exists()


@pytest.mark.parametrize('caller,chrm', [
    ('Deepvariant', 'Yes'), ('HaplotypeCaller', 'No'), ('BOTH', 'Yes'),
])
def test_protected_root_and_worker_enforce_filters_in_all_imported_modules(tmp_path, caller, chrm):
    fixture(tmp_path)
    path = frozen(tmp_path, caller, chrm)
    for mode in ('default', 'remote'):
        output = command(tmp_path, '--list-rules', '--mode', mode,
                         '--config', f'caller={caller}', f'chrM={chrm}', manifest=path)
        assert 'Reusing 1 completed samples; 1 samples eligible' in output
        assert 'Verified finished-sample exclusion' in output
        assert 'Legacy worker' not in output


@pytest.mark.skipif(not shutil.which('AdapterRemoval') or not shutil.which('pigz'),
                    reason='Real AdapterRemoval and pigz required')
@pytest.mark.parametrize('protected', [False, True])
def test_adapter_worker_executes_with_legacy_or_frozen_controller_selection(tmp_path, protected):
    # Production has this extension installed. Build before entering the legacy
    # tee pipeline: setuptools' stdout is not FASTQ data. This test intentionally
    # does not change the adapter implementation as part of the restart upgrade.
    subprocess.run(
        ['python3', '-c', 'from fastcheck_loader import ensure_fastcheck; '
         'assert ensure_fastcheck()[0]'],
        env=dict(os.environ, PYTHONPATH=str(REPO / 'scripts')),
        check=True, capture_output=True, text=True, timeout=90,
    )
    fixture(tmp_path)
    path = frozen(tmp_path) if protected else None
    target = 'fq/TEST_A.TEST_A_rg0.fastq.cut_1.fq.gz'
    output = command(
        tmp_path, target, '--allowed-rules', 'adapter_removal',
        '--mode', 'remote', '--force-use-threads', '--target-files-omit-workdir-adjustment',
        manifest=path,
    )
    assert 'localrule adapter_removal:' in output
    assert gzip.open(tmp_path / target, 'rt').read().startswith('@pair/1')
    assert (tmp_path / 'stats/TEST_A.TEST_A_rg0.fastq.stats.tsv').is_file()
    assert not (tmp_path / '.snakemake/restart_manifests').exists()
