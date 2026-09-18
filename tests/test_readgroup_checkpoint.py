import importlib.util
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "readgroup_checkpoint_under_test", REPO / "readgroup_checkpoint.py"
)
readgroup_checkpoint = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(readgroup_checkpoint)


class FakeCheckpoint:
    def __init__(self):
        self.samples = []

    def get(self, *, sample):
        self.samples.append(sample)
        return SimpleNamespace(output=[f"sampleinfo/{sample}.dat"])


def test_shared_checkpoint_returns_canonical_output():
    checkpoint = FakeCheckpoint()
    readgroup_checkpoint.register(checkpoint)

    assert readgroup_checkpoint.output("SAMPLE_1") == "sampleinfo/SAMPLE_1.dat"
    assert checkpoint.samples == ["SAMPLE_1"]


def test_readgroup_dependent_modules_use_shared_checkpoint():
    for name in ("Aligner.smk", "Stat.smk", "Kraken.smk"):
        source = (REPO / name).read_text()
        assert "_readgroup_checkpoint.output(sample)" in source
        assert "checkpoints.get_readgroups.get" not in source


def test_aligner_is_imported_once_before_consumers():
    source = (REPO / "Snakefile").read_text()
    assert source.count("use rule * from Aligner") == 1
    assert source.index("use rule * from Aligner") < source.index(
        "use rule * from Stat"
    )
    assert source.index("use rule * from Aligner") < source.index(
        "use rule * from Kraken"
    )
