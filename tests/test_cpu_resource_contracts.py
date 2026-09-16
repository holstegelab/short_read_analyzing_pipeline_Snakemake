from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def rule_block(filename, rule_name):
    text = (REPO / filename).read_text()
    marker = f"rule {rule_name}:"
    start = text.index(marker)
    next_rule = text.find("\nrule ", start + len(marker))
    return text[start : next_rule if next_rule >= 0 else None]


def checkpoint_block(filename, checkpoint_name):
    text = (REPO / filename).read_text()
    marker = f"checkpoint {checkpoint_name}:"
    start = text.index(marker)
    next_rule = text.find("\nrule ", start + len(marker))
    return text[start : next_rule if next_rule >= 0 else None]


def test_input_independent_cpu_reservations_follow_observed_average():
    markdup = rule_block("Aligner.smk", "markdup")
    cram = rule_block("Encrypt.smk", "cram_encrypt_fused")
    assert "n=get_n_merge_markdup" in markdup
    assert "--markdup-cores 0.95" in markdup
    assert 'n="1.85"' in cram
    assert "--encrypt-cores 0.45" in cram
    assert 'dcache_upload_slots="0.05"' in cram


def test_fused_cram_keeps_two_samtools_threads_separate_from_scheduler_reservation():
    block = rule_block("Encrypt.smk", "cram_encrypt_fused")
    assert "use_threads=2" in block
    assert "--cram-threads {resources.use_threads}" in block


def test_split_cpu_reservation_only_grows_for_real_split_or_sanitization():
    aligner = (REPO / "Aligner.smk").read_text()
    start = aligner.index("def get_n_split_alignments")
    end = aligner.index("rule split_alignments_by_readgroup", start)
    helper = aligner[start:end]

    assert "wildcards['filename'] in rg['file']" in helper
    assert "len(readgroups) > 1 or sinfo.get('erf_correct', False)" in helper
    assert 'return "3.0"' in helper
    assert 'return "0.9"' in helper

    start = aligner.index("def get_threads_split_alignments")
    end = aligner.index("rule split_alignments_by_readgroup", start)
    helper = aligner[start:end]
    assert "ZSLURM_LEASE_MAX_CORES" in helper
    assert "float(lease_max_cores) < 2.0" in helper
    assert "return 0" in helper
    assert "return 3" in helper


def test_simple_cpu_reservations_follow_observed_average():
    expected = {
        ("Aligner.smk", "merge_rgs_badmap"): 'n="0.7"',
        ("Stat.smk", "tar_stats_per_sample"): 'n="0.55"',
        ("Stat.smk", "chrM_and_numt_read_stats"): 'n="0.6"',
        ("Stat.smk", "whatsHap_phase_stats"): 'n="0.5"',
        ("Stat.smk", "tar_badmap_fastqs"): 'n="0.65"',
        ("Stat.smk", "copy_badmap_to_dcache"): 'n="0.5"',
    }
    for (filename, rule_name), reservation in expected.items():
        assert reservation in rule_block(filename, rule_name)
    assert 'n="0.5"' in checkpoint_block("Aligner.smk", "get_readgroups")


def test_fused_cpu_reservations_keep_tool_parallelism_separate():
    external = rule_block("Aligner.smk", "external_adapter_fused")
    assert 'n="5"' in external
    assert "--extract-cores" not in external
    assert "--adapter-cores 5" in external

    kmer = rule_block("Aligner.smk", "kmer_sex_fused")
    assert 'n="1.6"' in kmer
    assert "use_threads=KMC_RESERVED_CORES" in kmer
    assert "--kmc-threads {resources.use_threads}" in kmer
    assert "--low-cores 1.0" in kmer

    extract = rule_block("chrM_analysis.smk", "chrm_extract_align_fused")
    assert 'n="2.0"' in extract
    assert "threads_per_tool=2" in extract

    tail = rule_block("chrM_analysis.smk", "chrm_mutect_tail_fused")
    assert 'n="1.1"' in tail


def test_fractional_reservations_do_not_reduce_integer_tool_threads():
    merge_markdup = rule_block("Aligner.smk", "markdup")
    assert "use_threads=3" in merge_markdup
    assert "--merge-threads {resources.use_threads}" in merge_markdup

    split = rule_block("Aligner.smk", "split_alignments_by_readgroup")
    assert "n=get_n_split_alignments" in split
    assert "use_threads=get_threads_split_alignments" in split
    assert "-@ {resources.use_threads}" in split
    assert "--threads {resources.use_threads}" in split
    assert "-@ {resources.n}" not in split
    assert "--threads {resources.n}" not in split

    chrm_stats = rule_block("Stat.smk", "chrM_and_numt_read_stats")
    assert "use_threads=1" in chrm_stats
    assert "n = resources.use_threads" in chrm_stats


def test_route_and_completion_reservations_are_updated():
    aligner = (REPO / "Aligner.smk").read_text()
    assert "'s3': '0.5'" in aligner
    snakefile = (REPO / "Snakefile").read_text()
    assert snakefile.count('n="0.5"') >= 3
