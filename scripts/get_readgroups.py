#!/usr/bin/env python3
"""Resolve sample read groups inside the rule's Conda environment."""

import sys
from pathlib import Path


pipeline_root = Path(str(snakemake.params.pipeline_root)).resolve()
sys.path.insert(0, str(pipeline_root))

import read_samples  # noqa: E402
import utils  # noqa: E402


sample, warnings = read_samples.get_readgroups(
    snakemake.params.sample,
    str(snakemake.params.prefixpath),
)

output_path = Path(str(snakemake.output[0]))
output_path.parent.mkdir(parents=True, exist_ok=True)

if warnings:
    warning_path = Path(str(snakemake.params.warningfile))
    warning_path.parent.mkdir(parents=True, exist_ok=True)
    with warning_path.open('w', encoding='utf-8') as handle:
        for warning in warnings:
            print(f"WARNING: {warning}")
            handle.write(f"{warning}\n")

utils.save(sample, str(output_path))
