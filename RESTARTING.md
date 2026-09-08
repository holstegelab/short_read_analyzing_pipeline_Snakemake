# Restarting completed samples safely

For the `gVCF` endpoint, completed samples are reused by default. A `.finished`
marker selects a candidate, then startup validates the retained cohort inputs
and CRAM-upload receipt. Missing QC files are recovered from the local
`stats/<sample>.stats.tar.gz`, preserving archive timestamps. No source download
or tape staging is used for that recovery. Missing unarchiveable results produce
a report at `.snakemake/restart_missing_products.json` and stop startup.

The startup decision is frozen in `.snakemake/restart_manifests/<uuid>.json`.
`SHORT_READ_RESTART_MANIFEST` is inherited by worker Snakemake invocations and
declared in `envvars`; ZSlurm also transfers the submission environment. The
manifest checks the work directory and sample list. Do not reuse a manifest
after changing the sample list or the intended rebuild selection.

## How the dependency cut works

- Cohort input lists still contain **all** samples.
- Per-sample processing rules only match samples outside the frozen reuse set.
  Thus saved QC/gVCF files and completion receipts enter the DAG as existing
  inputs, without producers that could reconstruct the alignment pipeline.
- A reused sample's `retrieve_batch()` returns no dependency. Batch membership
  and numbering stay stable, while staging excludes reused samples.
- The inexpensive badmap archive/upload rules remain eligible for finished
  samples with uploadable local data. Missing badmap data cannot cause alignment
  to run again for an adopted sample.
- Old incomplete records cannot force excluded producers back into the graph.
  Incomplete records for unfinished samples retain their usual retry behavior.
- After loading all modules, a guard inspects every parsed per-sample output
  constraint. A rule-specific override that still accepts a reused sample stops
  startup before DAG construction. This includes the external-adapter override.

Small QC results, readgroup metadata, final cohort gVCFs, per-sample stats
archives, and upload/completion receipts are no longer temporary. Large raw
inputs, FASTQs, BAMs, intermediate variant files, and uploaded cohort tar files
remain temporary. Per-sample QC archives and coverage are also excluded from the
old `pipeline.done` startup cleanup, so a subsequent restart retains its recovery
inputs. This does not change output names or computational results.

## Existing interrupted runs

Stop the root controller and its jobs before repairing outputs. Do not run two
controllers against the same directory. Prepare explicitly if desired:

```bash
python /path/to/pipeline/scripts/prepare_restart.py --workdir /path/to/run
```

Pass the returned path with `--config restart_manifest=/absolute/manifest.json`
for a dry-run and the corresponding restart. Omit that option on a later restart
to freeze a fresh completed set. Do not export a stale manifest in a shell profile.

If a redundant rerun overwrote completed statistics, the recovery command can
also take `--replace-newer --backup-dir /absolute/new-backup-directory`. This
explicit mode restores archived files whose mtimes are later than `.finished`,
backing up each overwritten file. Normal startup restores only missing files.

`--config restart_rebuild_samples=sample1,sample2` exempts selected samples from
reuse and makes them eligible for normal processing; it does not force existing
outputs to be recomputed. Inspect a dry-run and select explicit targets/force
options for deliberate reprocessing. Stage sizes are refreshed without batch
renumbering. `reuse_finished_samples=False` opts out entirely, restoring legacy
DAG reconstruction; it is unsafe for the restart failure described above.

Prefer deploying before starting a run, or use a separate checkout for a stopped
workflow. This restart-protection upgrade explicitly supports workers from an
already-running, pre-upgrade controller: a worker without a manifest keeps the
legacy sample selection and never scans `.finished` markers or restores files.
The controller's existing DAG and temporary-file cleanup decisions remain in
effect until that run stops. Protection starts with the next root invocation;
editing the checkout does not retrofit protection into an in-memory controller.

Workers with a manifest always use its frozen selection, even if more samples
have since finished. A missing or invalid explicitly specified manifest is an
error, not a reason to compute a fresh selection on a worker. This compatibility
path is specific to the restart upgrade, not permission for arbitrary live rule
or input/output changes. Output filenames, computational commands and CPU/memory
requests are unchanged. An explicitly selected rebuild of a completed sample
does acquire a fresh active-storage lifecycle reservation; unfinished samples
do not acquire that reservation twice.

## Regression checks

`python -m pytest -q tests/test_restart_state.py tests/test_temp_lifecycle_contracts.py`
tests archived recovery and backups, unsafe archive rejection, frozen worker
selection, stable batch filtering, and a real Snakemake DAG with shared staging,
an old incomplete start, and both present/missing finished QC.
