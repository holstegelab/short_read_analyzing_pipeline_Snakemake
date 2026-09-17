# Release baseline

This release freezes the fused-only workflow on the upstream main line and
makes installation-specific paths explicit. It was merged to `main` after the
isolated Snellius and Spider canaries. A merge does not replace code already
loaded by running Snakemake controllers or ZSlurm managers; deploying it still
requires an explicit checkout/update and, for ZSlurm, a planned new manager.

The machine-readable component pins and resource archive hashes are in
[`deployment/component-lock.yaml`](deployment/component-lock.yaml). The
immutable tag `codex-spider-portability-2026-09-17.2` identifies the merged,
tested portability release; use the exact component commits in the lock
rather than the moving heads of their repositories.

## Scope of this baseline

- Pre-merge `origin/main` at `dc6cd25`, retaining the resource audit anchored
  at `c386a7d`, plus the active-storage fusion and portability commits listed
  in this release history.
- The reviewed fused-only cleanup, rebased as `8720c46`; predecessor rule
  implementations remain removed.
- Site configuration schema 1, loaded before `common` and `constants` in both
  controller and worker reparsing.
- The ZSlurm executor environment fix at `9ed793e`, so site/restart/storage
  variables are transported in the XML-RPC environment rather than a shell
  prefix.
- Existing resource releases kept separate from the encrypted Crypt4GH key
  archive. No tokens, key contents or passwords are committed here.

The SHA-256 file under `deployment/` pins the environment specifications,
post-deploy scripts and native build recipe. These are reproducibility inputs,
not solver-complete Conda lockfiles: the YAML files still contain several
unpinned dependencies. Preserve the created Conda package records for a
production release or generate platform locks before long-term archival.

## Validation boundary

The original baseline pipeline run passed **169 tests with no failures or
skips** while
comparing representative complete DAGs and real adapter outputs against a
fresh `c386a7d` archive. A separate explicit-Snellius-site run passed the real
root-workflow/worker-reparse checks. `git diff --check` passed.

The Git-ignored `fastcheck`, `fastcheck_hts`, `bam_merge`, `bam_dechimer` and
`fix_bam_rg_pairs` artifacts were compiled in this fresh worktree and both
extensions were imported under the preprocessing Python/htslib ABI. The native
executor suite passed **37 tests**, including literal spaces, quotes and
newlines in the shell-free worker environment.

The hg19/b37+chrY supplement was inspected after creation and uploaded with the
large-file checksum timeout. `dcache_cp` verified it against dCache; the remote
Adler-32 is recorded in the component lock.

These checks validate unchanged Snellius defaults, explicit Snellius site
configuration, Spider path resolution and worker environment transport. They
are not an end-to-end Spider certification.

The Spider portability release implements the allocation-scratch contract
described here. Its publication run passed **195 tests with 25 skips**, and
ZSlurm passed **162 tests**. Real Spider pilots subsequently exercised the
manager/chief RPC, dynamic lease resize, isolated scratch cleanup and a
two-rule paired-FASTQ preprocessing DAG. The same DAG completed on Snellius
with the preceding Snakemake fork commit `83c79c1d`, confirming that existing
installations do not require an immediate Snakemake update. Full alignment,
DeepVariant, dCache publication, restart and Apptainer canaries are still
required before Spider production.
Archive commands (`daget`, `dals`, `darelease`) remain Snellius-only unless a
target backend is installed and tested.

See [`DEPLOYMENT.md`](DEPLOYMENT.md) for the installation sequence and a
readiness checklist.
