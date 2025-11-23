from setuptools import setup, Extension
import os
import shlex
import subprocess


def _pkg_config_tokens(pkg: str, flag: str):
    try:
        out = subprocess.check_output(["pkg-config", flag, pkg], universal_newlines=True)
        return shlex.split(out.strip())
    except Exception:
        return []


def _parse_cflags(tokens):
    include_dirs = []
    extra_compile_args = []
    define_macros = []
    for t in tokens:
        if t.startswith("-I"):
            include_dirs.append(t[2:])
        elif t.startswith("-D"):
            if "=" in t[2:]:
                k, v = t[2:].split("=", 1)
                define_macros.append((k, v))
            else:
                define_macros.append((t[2:], None))
        else:
            extra_compile_args.append(t)
    return include_dirs, extra_compile_args, define_macros


def _parse_libs(tokens):
    library_dirs = []
    libraries = []
    extra_link_args = []
    for t in tokens:
        if t.startswith("-L"):
            library_dirs.append(t[2:])
        elif t.startswith("-l"):
            libraries.append(t[2:])
        else:
            extra_link_args.append(t)
    return library_dirs, libraries, extra_link_args


cflags = _pkg_config_tokens("htslib", "--cflags")
libs = _pkg_config_tokens("htslib", "--libs")

inc, extra_c, defs = _parse_cflags(cflags)
lib_dirs, libs_list, extra_l = _parse_libs(libs)

# Fallbacks if pkg-config is unavailable/missing
HTSLIB_INCLUDE = os.environ.get("HTSLIB_INCLUDE")
HTSLIB_LIBDIR = os.environ.get("HTSLIB_LIBDIR")
if HTSLIB_INCLUDE and HTSLIB_INCLUDE not in inc:
    inc.append(HTSLIB_INCLUDE)
if HTSLIB_LIBDIR and HTSLIB_LIBDIR not in lib_dirs:
    lib_dirs.append(HTSLIB_LIBDIR)

# Common conda layout fallback: use CONDA_PREFIX/include and CONDA_PREFIX/lib
CONDA_PREFIX = os.environ.get("CONDA_PREFIX") or os.environ.get("VIRTUAL_ENV")
if CONDA_PREFIX:
    cp_inc = os.path.join(CONDA_PREFIX, "include")
    cp_lib = os.path.join(CONDA_PREFIX, "lib")
    if os.path.isdir(cp_inc) and cp_inc not in inc:
        inc.append(cp_inc)
    if os.path.isdir(cp_lib) and cp_lib not in lib_dirs:
        lib_dirs.append(cp_lib)

if not libs_list:
    # Minimal set; conda-built htslib should pull deps transitively, but add common ones just in case
    libs_list = ["hts", "z", "bz2", "lzma"]

ext = Extension(
    name="fastcheck_hts",
    sources=["fastcheck_hts.c"],
    include_dirs=inc,
    library_dirs=lib_dirs,
    libraries=libs_list,
    extra_compile_args=["-O3", "-std=c11", "-Wall", "-Wextra"] + extra_c,
    extra_link_args=extra_l,
    define_macros=defs,
)

setup(
    name="fastcheck_hts",
    version="0.1.0",
    description="htslib-based direct BAM/CRAM checksum and stats",
    ext_modules=[ext],
)
