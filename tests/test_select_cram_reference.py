import importlib.util
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "select_cram_reference.py"
SPEC = importlib.util.spec_from_file_location("select_cram_reference", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def sq_header(chr1_m5, chry_m5):
    return (
        f"@SQ\tSN:chr1\tLN:1\tM5:{chr1_m5}\n"
        f"@SQ\tSN:chrY\tLN:1\tM5:{chry_m5}\n"
    )


def choose(header):
    return MODULE.choose_reference(
        header,
        primary_reference="primary.fa",
        hg19_reference="hg19.fa",
        hg19_b37_chry_reference="hg19_b37chrY.fa",
        hg38_reference="hg38.fa",
    )


def test_selects_hg38_from_header_m5():
    selected, reason = choose(
        sq_header(MODULE.HG38_CHR1_M5, MODULE.HG38_CHRY_M5)
    )
    assert selected == "hg38.fa"
    assert reason == "hg38 M5 signature"


def test_selects_standard_hg19_from_header_m5():
    selected, reason = choose(
        sq_header(MODULE.HG19_CHR1_M5, MODULE.HG19_CHRY_M5)
    )
    assert selected == "hg19.fa"
    assert reason == "hg19 M5 signature"


def test_selects_hybrid_hg19_b37_chry_from_header_m5():
    selected, reason = choose(
        sq_header(MODULE.HG19_CHR1_M5, MODULE.HG19_B37_CHRY_M5)
    )
    assert selected == "hg19_b37chrY.fa"
    assert reason == "hg19+b37-chrY M5 signature"


def test_unknown_signature_keeps_configured_reference():
    selected, reason = choose(sq_header("unknown", "unknown"))
    assert selected == "primary.fa"
    assert reason == "unknown M5 signature; keeping configured reference"


def test_aligner_uses_selector_instead_of_exception_fallback():
    aligner = (REPO / "Aligner.smk").read_text()
    block_start = aligner.index("rule split_alignments_by_readgroup:")
    block_end = aligner.index("def get_aligned_readgroup_folder", block_start)
    block = aligner[block_start:block_end]

    assert "scripts/select_cram_reference.py" in block
    assert "[split TEMP-FALLBACK]" not in block
    assert "except Exception" not in block
