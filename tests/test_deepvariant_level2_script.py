import os
import runpy
import tarfile
from pathlib import Path
from types import SimpleNamespace


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "extract_and_tar_deepvariant_level2.py"
)


def test_extract_and_tar_uses_activated_bcftools_and_preserves_inputs(
    tmp_path, monkeypatch
):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    bcftools = bin_dir / "bcftools"
    bcftools.write_text(
        """#!/usr/bin/env python3
import pathlib
import sys

if sys.argv[1] == "view":
    output = pathlib.Path(sys.argv[sys.argv.index("-o") + 1])
    output.write_bytes(b"fake-bgzf")
elif sys.argv[1] == "index":
    pathlib.Path(sys.argv[-1] + ".tbi").write_bytes(b"fake-index")
else:
    raise SystemExit(2)
""",
        encoding="utf-8",
    )
    bcftools.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")

    interval = tmp_path / "region.bed"
    interval.write_text("chr1\t0\t10\n", encoding="utf-8")
    sources = []
    for sample in ("sample_a", "sample_b"):
        source = tmp_path / f"{sample}.vcf.gz"
        index = tmp_path / f"{sample}.vcf.gz.tbi"
        source.write_bytes(b"source")
        index.write_bytes(b"index")
        sources.extend([source, index])

    output = tmp_path / "out" / "cohort.R0.dv.wes.gvcf.tar"
    fake_snakemake = SimpleNamespace(
        params=SimpleNamespace(
            samples=["sample_a", "sample_b"],
            region="R0",
            samplefile="cohort",
            dataset="wes",
            interval=str(interval),
        ),
        input=SimpleNamespace(gvcfs=[str(path) for path in sources]),
        output=SimpleNamespace(tar=str(output)),
    )

    runpy.run_path(str(SCRIPT), init_globals={"snakemake": fake_snakemake})

    assert output.is_file()
    with tarfile.open(output, "r") as handle:
        assert set(handle.getnames()) == {
            "sample_a.R0.dv.wes.g.vcf.gz",
            "sample_a.R0.dv.wes.g.vcf.gz.tbi",
            "sample_b.R0.dv.wes.g.vcf.gz",
            "sample_b.R0.dv.wes.g.vcf.gz.tbi",
        }
    assert all(path.is_file() for path in sources)
