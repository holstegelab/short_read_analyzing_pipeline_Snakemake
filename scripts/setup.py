from setuptools import setup, Extension

fastcheck_ext = Extension(
    name="fastcheck",
    sources=["fastcheck.c"],
    extra_compile_args=["-O3", "-std=c11", "-Wall", "-Wextra"],
)

setup(
    name="fastcheck",
    version="0.1.0",
    description="C extension to accelerate checksum comparisons between BAM/SAM and FASTQ pairs",
    ext_modules=[fastcheck_ext],
)
