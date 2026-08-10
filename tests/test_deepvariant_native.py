import importlib.util
import io
from pathlib import Path
import tarfile
import zipfile

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts/prepare_deepvariant_native.py"
DEEPVARIANT_RULES = Path(__file__).parents[1] / "Deepvariant.smk"
SPEC = importlib.util.spec_from_file_location("prepare_deepvariant_native", SCRIPT)
NATIVE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NATIVE)


def add_tar_file(archive, name, contents=b""):
    member = tarfile.TarInfo(name)
    member.size = len(contents)
    archive.addfile(member, io.BytesIO(contents))


def test_safe_layer_member_rejects_paths_outside_rootfs():
    with pytest.raises(ValueError):
        NATIVE.safe_layer_member("../../outside")
    with pytest.raises(ValueError):
        NATIVE.safe_layer_member("/absolute")


def test_apply_layer_honours_whiteouts(tmp_path):
    rootfs = tmp_path / "rootfs"
    (rootfs / "opaque").mkdir(parents=True)
    (rootfs / "obsolete").write_text("old")
    (rootfs / "opaque/old").write_text("old")
    layer = tmp_path / "layer.tar"
    with tarfile.open(layer, "w") as archive:
        add_tar_file(archive, ".wh.obsolete")
        add_tar_file(archive, "opaque/.wh..wh..opq")
        add_tar_file(archive, "opaque/new", b"new")

    NATIVE.apply_layer(layer, rootfs)

    assert not (rootfs / "obsolete").exists()
    assert not (rootfs / "opaque/old").exists()
    assert (rootfs / "opaque/new").read_text() == "new"
    assert not list(rootfs.rglob(".wh.*"))


def test_patch_bazel_zip_relocates_second_stage_python(tmp_path):
    archive_path = tmp_path / "binary.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("__main__.py", "PYTHON_BINARY = '/usr/bin/python3'\n")
        archive.writestr("runfiles/payload", b"payload")

    launcher = tmp_path / "runtime/bin/python"
    NATIVE.patch_bazel_zip(archive_path, launcher)

    with zipfile.ZipFile(archive_path) as archive:
        main = archive.read("__main__.py").decode()
        assert f"PYTHON_BINARY = {str(launcher)!r}" in main
        assert archive.read("runfiles/payload") == b"payload"


def test_launcher_reproduces_image_environment_without_leaking_library_path(tmp_path):
    rootfs = tmp_path / "rootfs"
    loader = rootfs / "lib/x86_64-linux-gnu/ld-linux-x86-64.so.2"
    python = rootfs / "usr/bin/python3.10"
    loader.parent.mkdir(parents=True)
    python.parent.mkdir(parents=True)
    loader.write_bytes(b"loader")
    python.write_bytes(b"python")
    loader.chmod(0o755)
    python.chmod(0o755)

    prefix = Path("/runtime/deepvariant-1.9.0")
    launcher = NATIVE.launcher_text(
        rootfs,
        prefix,
        prefix / "rootfs/opt/deepvariant/bin/run_deepvariant_native.py",
    )

    assert "export TF_ENABLE_ONEDNN_OPTS=1" in launcher
    assert "export VERSION=1.9.0" in launcher
    assert "unset PYTHONPATH" in launcher
    assert "LD_LIBRARY_PATH" in launcher
    assert "export LD_LIBRARY_PATH" not in launcher


def snakemake_rule_body(name):
    text = DEEPVARIANT_RULES.read_text()
    return text.split(f"rule {name}:", 1)[1].split("\nrule ", 1)[0]


def test_production_deepvariant_rule_is_namespace_free():
    rule = snakemake_rule_body("deepvariant")
    assert "get_deepvariant_native_runner" in rule
    assert "DEEPVARIANT" in rule
    assert "container:" not in rule


def test_apptainer_deepvariant_is_an_isolated_fallback():
    rule = snakemake_rule_body("deepvariant_apptainer")
    assert "docker://google/deepvariant:1.9.0" in rule
    assert "DEEPVARIANT_APPTAINER" in rule
