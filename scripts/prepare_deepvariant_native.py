#!/usr/bin/env python3
"""Prepare a namespace-free DeepVariant runtime from the official image.

The resulting runtime executes the Python binaries from the image with the
image's own dynamic loader and shared libraries.  It does not use Apptainer,
Singularity, Docker, Podman, chroot, or a user namespace at run time.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path, PurePosixPath
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile


DEFAULT_IMAGE = "docker://google/deepvariant:1.9.0"
RUNTIME_VERSION = 1
PYTHON_WRAPPERS = {
    "make_examples": "make_examples.zip",
    "call_variants": "call_variants.zip",
    "postprocess_variants": "postprocess_variants.zip",
    "vcf_stats_report": "vcf_stats_report.zip",
    "show_examples": "show_examples.zip",
    "runtime_by_region_vis": "runtime_by_region_vis.zip",
    "multisample_make_examples": "multisample_make_examples.zip",
    "labeled_examples_to_vcf": "labeled_examples_to_vcf.zip",
    "convert_to_saved_model": "convert_to_saved_model.zip",
    "make_examples_somatic": "make_examples_somatic.zip",
    "train": "train.zip",
}


def run(command: list[str], announce: bool = True) -> None:
    if announce:
        print("+ " + shlex.join(command), flush=True)
    subprocess.run(command, check=True)


def safe_layer_member(name: str) -> PurePosixPath:
    member = PurePosixPath(name)
    if member.is_absolute() or ".." in member.parts:
        raise ValueError(f"Unsafe path in image layer: {name!r}")
    return member


def remove_path(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def apply_layer(layer: Path, rootfs: Path) -> None:
    """Apply one Docker layer, including OCI whiteout semantics."""
    whiteouts: list[tuple[Path, str]] = []
    with tarfile.open(layer, "r:*") as layer_tar:
        for member in layer_tar:
            relative = safe_layer_member(member.name)
            if member.islnk():
                safe_layer_member(member.linkname)
            if relative.name.startswith(".wh."):
                whiteouts.append((rootfs.joinpath(*relative.parent.parts), relative.name))

    # Whiteouts describe removals from lower layers and must be applied before
    # this layer's ordinary members are extracted.
    for parent, name in whiteouts:
        if name == ".wh..wh..opq":
            if parent.is_dir():
                for child in parent.iterdir():
                    remove_path(child)
        else:
            target = parent / name.removeprefix(".wh.")
            if target.exists() or target.is_symlink():
                remove_path(target)

    run(
        ["tar", "-xf", str(layer), "-C", str(rootfs), "--no-same-owner"],
        announce=False,
    )

    # GNU tar materializes regular-file whiteouts.  They are metadata, not
    # files in the merged root filesystem.
    for parent, name in whiteouts:
        marker = parent / name
        if marker.exists() or marker.is_symlink():
            remove_path(marker)


def unpack_docker_archive(archive: Path, build_dir: Path) -> Path:
    archive_dir = build_dir / "docker-archive"
    rootfs = build_dir / "rootfs"
    archive_dir.mkdir()
    rootfs.mkdir()
    run(["tar", "-xf", str(archive), "-C", str(archive_dir), "--no-same-owner"])

    manifest_path = archive_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if len(manifest) != 1:
        raise ValueError(f"Expected one image in {archive}, found {len(manifest)}")
    for layer_name in manifest[0]["Layers"]:
        relative = safe_layer_member(layer_name)
        layer = archive_dir.joinpath(*relative.parts)
        if not layer.is_file():
            raise FileNotFoundError(f"Image layer is missing: {layer}")
        apply_layer(layer, rootfs)

    shutil.rmtree(archive_dir)
    return rootfs


def unpack_sif(sif: Path, build_dir: Path) -> Path:
    listing = subprocess.run(
        ["apptainer", "sif", "list", str(sif)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    filesystem_id = None
    for line in listing.splitlines():
        fields = [field.strip() for field in line.split("|")]
        if len(fields) >= 5 and "FS (Squashfs" in fields[4]:
            filesystem_id = fields[0]
            break
    if filesystem_id is None:
        raise ValueError(f"No Squashfs filesystem partition found in {sif}")

    squashfs = build_dir / "rootfs.squashfs"
    print(f"+ apptainer sif dump {filesystem_id} {sif} > {squashfs}", flush=True)
    with squashfs.open("wb") as output:
        subprocess.run(
            ["apptainer", "sif", "dump", filesystem_id, str(sif)],
            check=True,
            stdout=output,
        )
    rootfs = build_dir / "rootfs"
    run(["unsquashfs", "-no-xattrs", "-d", str(rootfs), str(squashfs)])
    squashfs.unlink()
    return rootfs


def find_image_python(rootfs: Path) -> Path:
    candidates = sorted((rootfs / "usr/bin").glob("python3.*"), reverse=True)
    for candidate in candidates:
        version = candidate.name.removeprefix("python3.")
        if version.isdigit() and candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    raise FileNotFoundError("No executable /usr/bin/python3.X found in image")


def find_loader(rootfs: Path) -> Path:
    candidates = [
        rootfs / "lib/x86_64-linux-gnu/ld-linux-x86-64.so.2",
        rootfs / "lib64/ld-linux-x86-64.so.2",
    ]
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    raise FileNotFoundError("No x86_64 glibc dynamic loader found in image")


def runtime_library_dirs(rootfs: Path) -> list[Path]:
    candidates = [
        rootfs / "lib/x86_64-linux-gnu",
        rootfs / "usr/lib/x86_64-linux-gnu",
        rootfs / "usr/local/lib",
        rootfs / "opt/conda/lib",
        rootfs / "opt/conda/envs/bio/lib",
    ]
    return [path for path in candidates if path.is_dir()]


def launcher_text(
    build_rootfs: Path,
    prefix: Path,
    target: Path | None,
    unbuffered: bool = False,
    preserve_pythonpath: bool = False,
) -> str:
    rootfs = prefix / "rootfs"
    loader = rootfs / find_loader(build_rootfs).relative_to(build_rootfs)
    python = rootfs / find_image_python(build_rootfs).relative_to(build_rootfs)
    libraries = ":".join(
        str(rootfs / path.relative_to(build_rootfs))
        for path in runtime_library_dirs(build_rootfs)
    )
    unset_names = (
        "PYTHONSTARTUP PYTHONUSERBASE VIRTUAL_ENV CONDA_PREFIX "
        "CONDA_DEFAULT_ENV LD_LIBRARY_PATH"
    )
    if not preserve_pythonpath:
        unset_names = "PYTHONPATH " + unset_names
    python_flag = " -u" if unbuffered else ""
    target_arg = "" if target is None else " " + shlex.quote(str(target))
    return f"""#!/usr/bin/env bash
set -euo pipefail
unset {unset_names}
export PYTHONHOME={shlex.quote(str(rootfs / 'usr'))}
export PYTHONNOUSERSITE=1
export TF_ENABLE_ONEDNN_OPTS=1
export DV_GPU_BUILD=0
export VERSION=1.9.0
export PYTHON_VERSION=3.10
export PATH={shlex.quote(str(prefix / 'bin') + ':/usr/bin:/bin')}
exec {shlex.quote(str(loader))} --library-path {shlex.quote(libraries)} {shlex.quote(str(python))}{python_flag}{target_arg} "$@"
"""


def write_executable(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(0o755)


def patch_bazel_zip(path: Path, python_launcher: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with zipfile.ZipFile(path, "r") as source, zipfile.ZipFile(
        temporary, "w"
    ) as destination:
        for member in source.infolist():
            contents = source.read(member)
            if member.filename == "__main__.py":
                text = contents.decode()
                relocated = text.replace(
                    "PYTHON_BINARY = '/usr/bin/python3'",
                    f"PYTHON_BINARY = {str(python_launcher)!r}",
                )
                if relocated == text:
                    raise ValueError(f"No Bazel PYTHON_BINARY found in {path}")
                contents = relocated.encode()
            destination.writestr(member, contents)
    os.replace(temporary, path)


def prepare_launchers(build_dir: Path, final_prefix: Path) -> None:
    build_rootfs = build_dir / "rootfs"
    final_rootfs = final_prefix / "rootfs"
    deepvariant_bin = build_rootfs / "opt/deepvariant/bin"
    if not deepvariant_bin.is_dir():
        raise FileNotFoundError("Image has no /opt/deepvariant/bin directory")

    original_entrypoint = deepvariant_bin / "run_deepvariant.py"
    patched_entrypoint = deepvariant_bin / "run_deepvariant_native.py"
    entrypoint_text = original_entrypoint.read_text()
    # The upstream orchestrator uses absolute /opt paths for its binaries and
    # bundled models. Point those paths at the extracted root filesystem.
    entrypoint_text = entrypoint_text.replace(
        "/opt/deepvariant/bin", str(final_prefix / "bin")
    )
    entrypoint_text = entrypoint_text.replace(
        "/opt/models", str(final_rootfs / "opt/models")
    )
    entrypoint_text = entrypoint_text.replace(
        "/opt/smallmodels", str(final_rootfs / "opt/smallmodels")
    )
    patched_entrypoint.write_text(entrypoint_text)

    bin_dir = build_dir / "bin"
    bin_dir.mkdir()
    python_launcher = final_prefix / "bin/python"
    write_executable(
        bin_dir / "python",
        launcher_text(
            build_rootfs,
            final_prefix,
            target=None,
            preserve_pythonpath=True,
        ),
    )
    write_executable(
        bin_dir / "run_deepvariant",
        launcher_text(
            build_rootfs,
            final_prefix,
            final_rootfs / "opt/deepvariant/bin/run_deepvariant_native.py",
        ),
    )
    for command, zip_name in PYTHON_WRAPPERS.items():
        patch_bazel_zip(deepvariant_bin / zip_name, python_launcher)
        write_executable(
            bin_dir / command,
            launcher_text(
                build_rootfs,
                final_prefix,
                final_rootfs / "opt/deepvariant/bin" / zip_name,
                unbuffered=command in {"labeled_examples_to_vcf", "make_examples_somatic"},
            ),
        )

    parallel = build_rootfs / "usr/bin/parallel"
    if not parallel.is_file():
        raise FileNotFoundError("Image has no GNU Parallel executable")
    shutil.copy2(parallel, bin_dir / "parallel")
    (bin_dir / "parallel").chmod(0o755)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True, help="New runtime directory")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--image", help=f"OCI image (default: {DEFAULT_IMAGE})")
    source.add_argument("--docker-archive", type=Path, help="Existing skopeo Docker archive")
    source.add_argument("--sif", type=Path, help="Existing SIF image (useful on Spider)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    prefix = args.prefix.expanduser().resolve()
    if any(character.isspace() for character in str(prefix)):
        raise ValueError("The runtime prefix may not contain whitespace")
    if prefix.exists():
        raise FileExistsError(f"Refusing to replace existing runtime: {prefix}")

    prefix.parent.mkdir(parents=True, exist_ok=True)
    build_dir = Path(tempfile.mkdtemp(prefix=f".{prefix.name}.build-", dir=prefix.parent))
    source_description: str
    try:
        if args.sif:
            source_path = args.sif.expanduser().resolve(strict=True)
            source_description = f"sif:{source_path}"
            unpack_sif(source_path, build_dir)
        else:
            if args.docker_archive:
                archive = args.docker_archive.expanduser().resolve(strict=True)
                source_description = f"docker-archive:{archive}"
            else:
                image = args.image or DEFAULT_IMAGE
                archive = build_dir / "deepvariant-image.tar"
                source_description = image
                run(
                    [
                        "skopeo",
                        "copy",
                        "--retry-times",
                        "3",
                        image,
                        f"docker-archive:{archive}:deepvariant-native:1.9.0",
                    ]
                )
            unpack_docker_archive(archive, build_dir)
            if archive.parent == build_dir and archive.exists():
                archive.unlink()

        prepare_launchers(build_dir, prefix)
        marker = {
            "runtime_format": RUNTIME_VERSION,
            "deepvariant_version": "1.9.0",
            "source": source_description,
            "entrypoint": "bin/run_deepvariant",
        }
        (build_dir / ".deepvariant-native.json").write_text(
            json.dumps(marker, indent=2, sort_keys=True) + "\n"
        )
        build_dir.chmod(0o2750)
        os.replace(build_dir, prefix)
        run([str(prefix / "bin/run_deepvariant"), "--version"])
    except Exception:
        if build_dir.exists():
            shutil.rmtree(build_dir)
        # A failed smoke test can only leave the prefix created by this run.
        if prefix.exists() and not (prefix / ".deepvariant-native.ready").exists():
            shutil.rmtree(prefix)
        raise

    (prefix / ".deepvariant-native.ready").write_text("DeepVariant 1.9.0\n")
    print(f"Native DeepVariant runtime ready: {prefix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
