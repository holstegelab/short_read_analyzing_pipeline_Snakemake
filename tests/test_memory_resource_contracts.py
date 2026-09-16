from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def rule_block(filename, rule_name):
    text = (REPO / filename).read_text()
    marker = f"rule {rule_name}:"
    start = text.index(marker)
    next_rule = text.find("\nrule ", start + len(marker))
    return text[start : next_rule if next_rule >= 0 else None]


def test_large_fused_memory_targets_follow_phase_peak_estimates():
    aligner = (REPO / "Aligner.smk").read_text()
    assert "(attempt - 1) * 0.5 * 34200 + 34200" in rule_block(
        "Aligner.smk", "kmer_sex_fused"
    )
    assert "(attempt - 1) * 0.25 * 38000 + 38000" in rule_block(
        "Aligner.smk", "align_reads_fused"
    )
    assert "return 13500" in aligner
    assert "--low-memory-mb {params.low_memory_mb}" in rule_block(
        "Aligner.smk", "align_reads_fused"
    )
    assert "--low-memory-mb 4000" in rule_block(
        "Deepvariant.smk", "deepvariant_phasing_fused"
    )


def test_clear_overreservations_are_reduced():
    expected = {
        ("Stat.smk", "tar_badmap_fastqs"): "mem_mb=384",
        ("Stat.smk", "copy_badmap_to_dcache"): "mem_mb=512",
        ("chrM_analysis.smk", "chrm_extract_align_fused"): "mem_mb=2000",
    }
    for (filename, rule_name), target in expected.items():
        assert target in rule_block(filename, rule_name)


def test_clear_underreservations_are_raised():
    expected = {
        ("Aligner.smk", "merge_rgs_badmap"): "mem_mb=200",
        ("Stat.smk", "tar_stats_per_sample"): "mem_mb=200",
        ("chrM_analysis.smk", "chrm_mutect_tail_fused"): "mem_mb=2500",
    }
    for (filename, rule_name), target in expected.items():
        assert target in rule_block(filename, rule_name)
    snakefile = (REPO / "Snakefile").read_text()
    assert snakefile.count("mem_mb=200") >= 3


def test_split_alignments_base_memory_is_seven_and_a_half_gb():
    aligner = (REPO / "Aligner.smk").read_text()
    start = aligner.index("def get_mem_mb_split_alignments")
    end = aligner.index("rule split_alignments_by_readgroup", start)

    assert "res = 7500" in aligner[start:end]


def test_markdup_wgs_base_memory_is_three_gb():
    aligner = (REPO / "Aligner.smk").read_text()
    start = aligner.index("def get_mem_mb_markdup")
    end = aligner.index("def get_n_merge_markdup", start)

    assert "res = 3000 if 'wgs'" in aligner[start:end]
    assert "else 150" in aligner[start:end]


def test_storage_fusions_retain_phase_memory_peaks_and_only_shrink():
    markdup = rule_block("Aligner.smk", "markdup")
    cram = rule_block("Encrypt.smk", "cram_encrypt_fused")
    assert "mem_mb=get_mem_mb_merge_markdup" in markdup
    assert "--markdup-memory-mb {resources.mem_mb}" in markdup
    assert "mem_mb=1800" in cram
    assert "--encrypt-memory-mb 512" in cram


def test_runner_defaults_match_pipeline_lease_targets():
    deepvariant = (REPO / "scripts/run_fused_deepvariant_phasing.py").read_text()
    extract = (REPO / "scripts/run_fused_chrm_extract_align.py").read_text()
    tail = (REPO / "scripts/run_fused_chrm_tail.py").read_text()
    assert '"--low-memory-mb", type=float, default=4000' in deepvariant
    assert '"--memory-mb", type=int, default=2000' in extract
    assert '"--memory-mb", type=int, default=2500' in tail
