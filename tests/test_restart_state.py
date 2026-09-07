import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import restart_state as rs


def make_archive(root, members):
    (root/'stats').mkdir(exist_ok=True)
    with tarfile.open(root/'stats/A.stats.tar.gz', 'w:gz') as archive:
        for name, value in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(value)
            info.mtime = 100
            info.mode = 0o640
            archive.addfile(info, io.BytesIO(value))
    (root/'source').mkdir(exist_ok=True)
    (root/'source/A.finished').touch()
    os.utime(root/'source/A.finished', (200,200))


def test_restore_keeps_mtime_and_does_not_overwrite_existing(tmp_path):
    make_archive(tmp_path, {'kmer/A.result.yaml':b'sex: M\n', 'stats/A.qc':b'old'})
    (tmp_path/'stats/A.qc').write_bytes(b'current')
    restored = rs.restore_stats(tmp_path, 'A')
    assert restored == ['kmer/A.result.yaml']
    assert (tmp_path/'kmer/A.result.yaml').read_bytes() == b'sex: M\n'
    assert (tmp_path/'kmer/A.result.yaml').stat().st_mtime == 100
    assert (tmp_path/'stats/A.qc').read_bytes() == b'current'


def test_recover_overwritten_statistics_keeps_backup(tmp_path):
    make_archive(tmp_path, {'stats/A.qc':b'original'})
    (tmp_path/'stats/A.qc').write_bytes(b'partial rerun')
    rs.restore_stats(tmp_path, 'A', True, tmp_path/'backup')
    assert (tmp_path/'stats/A.qc').read_bytes() == b'original'
    assert (tmp_path/'backup/stats/A.qc').read_bytes() == b'partial rerun'


@pytest.mark.parametrize('name', ['../outside', '/absolute', 'stats/B.qc', 'source/A.finished', '.'])
def test_recovery_rejects_unrelated_or_escaping_members(tmp_path, name):
    make_archive(tmp_path, {name:b'bad'})
    with pytest.raises(ValueError, match='Unexpected'):
        rs.restore_stats(tmp_path, 'A')


def manifest(root, samples=('A','B')):
    path = root/'manifest.json'
    path.write_text(json.dumps({'schema_version':1, 'workdir':str(root),
        'sample_digest':rs.sample_digest(samples), 'reused_samples':['A'],
        'rebuild_samples':[]}))
    return path


def test_frozen_manifest_is_shared_with_workers(tmp_path, monkeypatch):
    path = manifest(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(rs.MANIFEST_ENV, str(path))
    before = rs.configure({}, {'A':{},'B':{}}, ([],[]))
    (tmp_path/'source').mkdir()
    (tmp_path/'source/B.finished').touch()
    after = rs.configure({}, {'A':{},'B':{}}, ([],[]), is_main_process=False)
    assert before == after
    assert before[0] == {'A'}
    assert before[2] == '(?:B)'
    with pytest.raises(ValueError, match='Sample list changed'):
        rs.load_manifest(path, tmp_path, ['A','B','C'])
    with pytest.raises(ValueError, match='settings changed'):
        rs.configure({'caller':'HaplotypeCaller'}, {'A':{},'B':{}}, ([],[]))


def test_legacy_worker_does_not_scan_or_restore_completed_samples(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv(rs.MANIFEST_ENV, raising=False)
    (tmp_path/'source').mkdir()
    (tmp_path/'source/A.finished').touch()
    def forbidden(*args, **kwargs):
        pytest.fail('A legacy worker must not run restart preflight')
    monkeypatch.setattr(rs, 'prepare', forbidden)
    reused, rebuild, pattern = rs.configure({}, {'A':{}}, ([],[]), is_main_process=False)
    assert not reused and not rebuild
    import re
    assert re.fullmatch(pattern, 'A')
    assert rs.MANIFEST_ENV not in os.environ
    assert not (tmp_path/'.snakemake').exists()


def test_new_root_exports_even_an_empty_manifest(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv(rs.MANIFEST_ENV, raising=False)
    before = rs.configure({}, {'A':{}}, ([],[]))
    path = Path(os.environ[rs.MANIFEST_ENV])
    assert json.loads(path.read_text())['reused_samples'] == []
    (tmp_path/'source').mkdir()
    (tmp_path/'source/A.finished').touch()
    assert rs.configure({}, {'A':{}}, ([],[]), is_main_process=False) == before


def test_worker_with_missing_explicit_manifest_never_recomputes(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(rs.MANIFEST_ENV, str(tmp_path/'missing.json'))
    with pytest.raises(FileNotFoundError):
        rs.configure({}, {'A':{}}, ([],[]), is_main_process=False)


def test_stage_filter_keeps_batch_membership_and_ids():
    samples = {'A':{'filesize':10, 'need_retrieval':True, 'from_external':'dcache'},
               'B':{'filesize':20, 'need_retrieval':False, 'from_external':'dcache'},
               'C':{'filesize':30, 'need_retrieval':True, 'from_external':'dcache'}}
    batches = {'cohort':{'dcache':[{'samples':['A','B'], 'size':10},
                                  {'samples':['C'], 'size':30}], 'archive':[]}}
    rs.refresh_retrieval(samples, batches, {'A'}, {'B'})
    assert batches['cohort']['dcache'] == [
        {'samples':['A','B'], 'size':20}, {'samples':['C'], 'size':30}]
    assert samples['A']['need_retrieval'] is False
    assert samples['B']['need_retrieval'] is True


def test_rebuild_reserves_again_only_for_completed_lifecycles():
    from types import SimpleNamespace
    aligner = (REPO/'Aligner.smk').read_text()
    source = 'def _start_sample_active_add' + aligner.split('def _start_sample_active_add',1)[1].split('\ndef _start_sample_tier_remove',1)[0]
    present = {'source/A.started'}
    ns = {'os':SimpleNamespace(path=SimpleNamespace(exists=lambda p:p in present)),
          'pj':os.path.join, 'SOURCEDIR':'source', 'REBUILD_SAMPLES':{'A'},
          'SAMPLEINFO':{'A':{'from_external':'dcache'}},
          'active_use_gb':lambda wc:123}
    exec(source, ns)
    assert ns['_start_sample_active_add']({'sample':'A'}) == 0
    present.add('source/A.finished')
    assert ns['_start_sample_active_add']({'sample':'A'}) == 123


def test_external_adapter_specific_filter_also_excludes_completed_samples():
    import re
    aligner = (REPO/'Aligner.smk').read_text()
    source = 'EXTERNAL_ALIGNMENT_SAMPLE_PATTERN =' + aligner.split('EXTERNAL_ALIGNMENT_SAMPLE_PATTERN =',1)[1].split('\nrule external_adapter_fused:',1)[0]
    ns = {'re':re, 'SAMPLEINFO':{'A':{'file_type':'cram'},'B':{'file_type':'cram'},
                                 'C':{'file_type':'fastq'}}, 'REUSED_SAMPLES':{'A'}}
    exec(source, ns)
    pattern = ns['EXTERNAL_ALIGNMENT_SAMPLE_PATTERN']
    assert re.fullmatch(pattern, 'B')
    assert not re.fullmatch(pattern, 'A')
    assert not re.fullmatch(pattern, 'C')


def test_parsed_rule_guard_rejects_local_constraint_overrides():
    from types import SimpleNamespace
    safe = SimpleNamespace(name='compute', output=['stats/{sample,(?:B)}.qc'])
    allowed = SimpleNamespace(name='copy_badmap_to_dcache', output=['badmap/{sample}.copied'])
    rs.verify_rule_filters([safe,allowed], {'A'})
    unsafe = SimpleNamespace(name='external_adapter_fused', output=['fq/{sample,(?:A|B)}.fq'])
    with pytest.raises(ValueError, match='external_adapter_fused still accepts reused sample A'):
        rs.verify_rule_filters([safe,unsafe], {'A'})


def test_native_adapter_filter_excludes_finished_fastq_samples():
    import re
    text = (REPO / 'Aligner.smk').read_text()
    source = 'NATIVE_FASTQ_SAMPLE_PATTERN =' + text.split('NATIVE_FASTQ_SAMPLE_PATTERN =', 1)[1].split('\nrule adapter_removal:', 1)[0]
    ns = {'re': re, 'REUSED_SAMPLES': {'A'}, 'SAMPLEINFO': {
        'A': {'file_type': 'fastq_paired'}, 'B': {'file_type': 'fastq_paired'},
        'C': {'file_type': 'cram'},
    }}
    exec(source, ns)
    pattern = ns['NATIVE_FASTQ_SAMPLE_PATTERN']
    assert re.fullmatch(pattern, 'B')
    assert not re.fullmatch(pattern, 'A')
    assert not re.fullmatch(pattern, 'C')


def test_preflight_recovers_before_freezing_and_reports_irrecoverable_files(tmp_path, monkeypatch):
    make_archive(tmp_path, {'kmer/A.result.yaml':b'sex: M\n', 'stats/A.qc':b'good'})
    (tmp_path/'fq_badmap').mkdir()
    (tmp_path/'fq_badmap/A.badmap.tar.copied').touch()
    required = ['source/A.finished', 'stats/A.qc', 'kmer/A.result.yaml']
    monkeypatch.setattr(rs, 'completion_products', lambda *a: list(required))
    path, state = rs.prepare(tmp_path, {'A':{}}, {}, ([],[]))
    assert state['reused_samples'] == ['A']
    assert (tmp_path/'stats/A.qc').read_bytes() == b'good'
    assert rs.load_manifest(path, tmp_path, ['A'])['reused_samples'] == ['A']
    required.append('gvcf/A.g.vcf.gz')
    with pytest.raises(ValueError, match='No samples were silently reprocessed'):
        rs.prepare(tmp_path, {'A':{}}, {}, ([],[]))
    report = tmp_path/'.snakemake/restart_missing_products.json'
    assert json.loads(report.read_text()) == {'A':['gvcf/A.g.vcf.gz']}
    required.pop()
    new, _ = rs.prepare(tmp_path, {'A':{}}, {}, ([],[]))
    assert not report.exists()
    assert new.with_name(new.stem+'.resolved_missing_products.json').exists()


@pytest.mark.parametrize('missing_qc', [False, True])
def test_restart_dag_cannot_reopen_finished_sample_via_shared_stage(tmp_path, missing_qc):
    path = manifest(tmp_path)
    for name in ['source/A.started','source/A.finished','cram/A.copied'] + ([] if missing_qc else ['stats/A.qc']):
        dest=tmp_path/name
        dest.parent.mkdir(exist_ok=True)
        dest.touch()
    # The old start's incomplete record must not force its producer back in.
    import base64
    incomplete = tmp_path/'.snakemake/incomplete'
    incomplete.mkdir(parents=True)
    (incomplete/base64.urlsafe_b64encode(b'source/A.started').decode()).write_text('{"external_jobid":"old"}')
    snake = tmp_path/'Snakefile'
    snake.write_text('''
from restart_state import configure
reused, rebuild, processing = configure({}, {'A':{},'B':{}}, ([],[]),
                                      is_main_process=workflow.is_main_process)
rule all:
    input: 'source/A.finished', 'source/B.finished', 'cohort.tsv'
rule stage:
    output: temp('fetch/batch')
    shell: 'touch {output}'
rule start:
    input: lambda wc: [] if wc.sample in reused else ancient('fetch/batch')
    output: 'source/{sample}.started', temp('source/{sample}.data')
    wildcard_constraints: sample=processing
    shell: 'touch {output}'
rule compute:
    input: 'source/{sample}.data'
    output: 'stats/{sample}.qc', temp('bams/{sample}.bam')
    wildcard_constraints: sample=processing
    shell: 'touch {output}'
rule upload:
    input: 'bams/{sample}.bam'
    output: 'cram/{sample}.copied'
    wildcard_constraints: sample=processing
    shell: 'touch {output}'
rule finish:
    input: 'cram/{sample}.copied', 'stats/{sample}.qc'
    output: 'source/{sample}.finished'
    wildcard_constraints: sample=processing
    shell: 'touch {output}'
rule cohort:
    input: 'stats/A.qc', 'stats/B.qc'
    output: 'cohort.tsv'
    shell: 'touch {output}'
''')
    env = dict(os.environ, PYTHONPATH=str(REPO), **{rs.MANIFEST_ENV:str(path)})
    env.pop('SNAKEMAKE_PROFILE',None)
    result=subprocess.run([sys.executable,'-m','snakemake','--cores','1','--dry-run','--rerun-incomplete'],
                          cwd=tmp_path,env=env,capture_output=True,text=True,timeout=45)
    text=result.stdout+result.stderr
    assert 'wildcards: sample=A' not in text
    if missing_qc:
        assert result.returncode != 0
        assert 'MissingInputException' in text
        assert 'stats/A.qc' in text
    else:
        assert result.returncode == 0, text
        assert 'wildcards: sample=B' in text
        # Simulate the nested worker's one-rule invocation after another job
        # creates B.finished. The worker must still see B's original rules.
        (tmp_path/'source/B.data').touch()
        (tmp_path/'source/B.finished').touch()
        worker = subprocess.run([
            sys.executable, '-m', 'snakemake', '--cores', '1', '--dry-run',
            '--target-jobs', 'compute:sample=B', '--allowed-rules', 'compute',
            '--force', '--ignore-incomplete', '--mode', 'remote'],
            cwd=tmp_path, env=env, capture_output=True, text=True, timeout=45)
        worker_text = worker.stdout + worker.stderr
        assert worker.returncode == 0, worker_text
        assert 'rule compute:' in worker_text
        assert 'wildcards: sample=B' in worker_text
        assert 'wildcards: sample=A' not in worker_text
