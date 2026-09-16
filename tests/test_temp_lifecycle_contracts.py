from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def rule_body(text: str, name: str) -> str:
    marker = f"rule {name}:"
    start = text.index(marker)
    next_rule = text.find("\nrule ", start + len(marker))
    return text[start:] if next_rule < 0 else text[start:next_rule]


def test_external_start_routes_publish_tracked_temp_directories():
    aligner = (REPO / "Aligner.smk").read_text()
    snakefile = (REPO / "Snakefile").read_text()
    archive = rule_body(aligner, "start_sample_archive")
    dcache = rule_body(aligner, "start_sample_dcache")
    s3 = rule_body(aligner, "start_sample_s3")
    release = rule_body(aligner, "release_materialized_source")

    assert 'materialized=temp(directory(pj(SOURCEDIR,"{sample}.data")))' in archive
    assert (
        'materialized=temp(directory(pj(SOURCEDIR,"{sample}.dcache_data")))'
        in dcache
    )
    assert "sample=START_SAMPLE_ARCHIVE_PATTERN" in archive
    assert "sample=START_SAMPLE_DCACHE_PATTERN" in dcache
    assert 'materialized=temp(directory(pj(SOURCEDIR,"{sample}.s3_data")))' in s3
    assert "sample=START_SAMPLE_S3_PATTERN" in s3
    assert "s3_download_slots=1" in s3
    assert "dcache_download_slots=1" not in s3
    assert "_start_sample_route(wildcards) == 'archive'" in aligner
    assert "files.append(ancient(external_data_dir(" in aligner
    assert 'bam=pj(BAM, "{sample}.markdup.bam")' in release
    assert 'materialized=lambda wildcards: external_data_dir(' in release
    assert 'f"{wildcards.sample}.materialized_consumed"' in snakefile


def test_raw_intermediates_are_temporary_but_completed_results_are_durable():
    aligner = (REPO / "Aligner.smk").read_text()
    snakefile = (REPO / "Snakefile").read_text()
    stats = (REPO / "Stat.smk").read_text()

    assert aligner.count(
        'badmap_fastq1=temp(pj(FQ_BADMAP,"{sample}.{readgroup}.badmap_R1.fastq.gz"))'
    ) == 1
    assert aligner.count(
        'badmap_fastq2=temp(pj(FQ_BADMAP,"{sample}.{readgroup}.badmap_R2.fastq.gz"))'
    ) == 1
    assert 'tar=pj(STAT,"{sample}.stats.tar.gz")' in stats
    assert "coverage_regions=pj(STAT, 'cov', '{sample}.regions.bed.gz')" in stats
    assert 'os.path.join(FQ_BADMAP, se + ".*.badmap_*.fastq.gz")' in snakefile
    assert 'os.path.join(STAT, "cov", se + ".*")' not in snakefile
    assert 'os.path.join(STAT, "*.stats_bundle.tar.gz")' in snakefile


def test_markdup_ssd_reservation_includes_merge_spill_and_final_bam():
    aligner = (REPO / "Aligner.smk").read_text()
    block = rule_body(aligner, "markdup")

    assert 'ssd_use="required"' in block
    assert "ssd_gb=get_ssd_gb_merge_markdup" in block
    assert "factor = 4.0 if len(input.bam) > 1 else 3.0" in aligner
    assert "overhead_gb=6, minimum_gb=12" in aligner


def test_cram_encryption_and_upload_have_no_active_storage_payload():
    encrypt = (REPO / "Encrypt.smk").read_text()
    block = rule_body(encrypt, "cram_encrypt_fused")

    assert 'bam=pj(BAM,"{sample}.markdup.bam")' in block
    assert 'copied=pj(CRAM,"{sample}.mapped_hg38.cram.copied")' in block
    assert 'sum=pj(CRAM,"{sample}.mapped_hg38.cram.ADLER32")' in block
    assert 'mapped_hg38.cram.c4gh"))' not in block
    assert 'mapped_hg38.cram.crai"))' not in block
    assert 'pj(CRAM,"{sample}.mapped_hg38.cram")' not in block
    assert "ssd_use=\"required\"" in block
    assert "factor=1.5, overhead_gb=6, minimum_gb=12" in block
    assert 'dcache_upload_slots="0.05"' in block
    assert "--upload-script" in block
    assert "rule copy_to_dcache:" not in encrypt


def test_cram_transfer_share_is_removed_from_active_reservation():
    common = (REPO / "common.py").read_text()
    assert "ACTIVE_CRAM_TRANSFER_FRAC_REMOVED = 0.15" in common
    assert "return ACTIVE_RESERVATION_FRAC * res" in common
    assert "def active_release_upload" not in common
    assert "0.30 / ACTIVE_RESERVATION_FRAC" in common
