from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def rule_block(filename, rule_name):
    text = (REPO / filename).read_text()
    marker = f"rule {rule_name}:"
    start = text.index(marker)
    next_rule = text.find("\nrule ", start + len(marker))
    return text[start : next_rule if next_rule >= 0 else None]


def test_low_io_rules_allow_normal_compute_nodes():
    optional_rules = {
        "Deepvariant.smk": (
            "deepvariant_phasing_fused",
        ),
        "Deepvariant_apptainer.smk": ("deepvariant_apptainer",),
        "chrM_analysis.smk": (
            "chrm_extract_align_fused",
            "chrm_mutect_tail_fused",
        ),
    }
    for filename, rules in optional_rules.items():
        for rule in rules:
            assert 'ssd_use="possible"' in rule_block(filename, rule)


def test_fused_optional_rules_pass_a_shared_scratch_fallback():
    for filename, rule in (
        ("Deepvariant.smk", "deepvariant_phasing_fused"),
        ("chrM_analysis.smk", "chrm_extract_align_fused"),
        ("chrM_analysis.smk", "chrm_mutect_tail_fused"),
    ):
        block = rule_block(filename, rule)
        assert "tmpdir=tmpdir" in block
        assert "--shared-scratch-base {resources.tmpdir:q}" in block


def test_high_io_fused_bam_qc_stays_ssd_required():
    assert 'ssd_use="required"' in rule_block("Stat.smk", "bam_qc_fused")
