# Release baseline

This checkout freezes the fused-only workflow on the current upstream main
line and makes installation-specific paths explicit. It is prepared in the
isolated branch `codex/spider-release-baseline-20260914`; it does not alter a
running Snellius workflow.

The machine-readable component pins and resource archive hashes are in
[`deployment/component-lock.yaml`](deployment/component-lock.yaml). The local
tag `codex-spider-portability-baseline-2026-09-14` identifies the tested final
pipeline commit. The tag and branch must be pushed deliberately before another
site can fetch them.

## Scope of this baseline

- Current `origin/main` at `c386a7d`, including the latest gCNV analysis and
  DeepVariant packaging changes.
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

The final pipeline run passed **169 tests with no failures or skips** while
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

The later Spider portability candidate implements the allocation-scratch
contract described here without changing this historical baseline tag. It has
unit coverage plus a real Spider allocation resolver canary, but still needs a
complete ZSlurm manager/worker and biological smoke run before production.
Archive commands (`daget`, `dals`, `darelease`) remain Snellius-only unless a
target backend is installed and tested.

See [`DEPLOYMENT.md`](DEPLOYMENT.md) for the installation sequence and a
readiness checklist.
