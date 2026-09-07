# Fused-only checkout: changes and validation

Date: 2026-09-07. Branch: `codex/fused-only-20260907`.
Checkout: `/gpfs/home2/hulsmanm/projects/short_read_pipeline_fused_only_20260907`.
Baseline snapshot commit: `c02883b` (ProjectMine restart checkout contents,
including the previously uncommitted pipeline and restart fixes).

The main and ProjectMine source worktrees were not edited. No production
workflow directory, queued job, running manager or engine was changed. Commits
in this branch are local; no push or production launch is part of this task.

## Changes

- Removed the 20 superseded rule definitions across the seven fused families,
  their toggles/ruleorder chains, obsolete helpers/runtime table entries,
  standalone sort runner/pilot and three obsolete environment definitions.
- Kept native FASTQ `adapter_removal`, with disjoint sample constraints and a
  small external runner. Its processing is shared with
  `external_adapter_fused` through `scripts/adapter_processing.py`.
- Kept simultaneous adapter identification and trimming, quality detection,
  optional deduplication and retry rescue; no thread/memory retuning.
- Moved 11 unchanged runtime functions/classes (scratch, atomic publication,
  leases) from the alignment runner to `scripts/pipeline_runtime.py`.
- Moved optional Apptainer comparison rules into `Deepvariant_apptainer.smk`,
  still included for the existing explicit target. Production DeepVariant
  remains native. Updated DRAGMAP/DeepVariant provenance lookups accordingly.
- Preserved source/restart protections and retained cohort products. Split,
  sample-level merge/markdup, CRAM/encryption and the other non-replaced stages
  remain. Custom C algorithms and their Python reference/fallback code remain.
- Updated current documentation; marked old optimization notes historical.
  Spider changes are a **plan only** in [PORTABILITY_PLAN.md](PORTABILITY_PLAN.md).

## Checks performed

Final complete test run, with real tools and frozen-baseline comparisons:
**130 passed, zero failures/skips, in 72.01 seconds**. `git diff --check` passed.

1. Baseline tests before edits: **101 passed, 5 skipped** (samtools/pigz were
   not on the initial test PATH). The subsequent complete run uses real tools.
2. AST snapshot comparison: **33 input/output/resource sections unchanged**
   across the seven fused stages plus native adapters, split, markdup and
   sample-level merge. This includes `temp`/`ensure` annotations and resource
   expressions, not just filenames.
3. The 11 extracted runtime definitions have identical ASTs before/after;
   scheduler-lease behavior is unchanged by moving them. Existing lease,
   profiling, environment and runner tests remain in the suite.
4. Real Snakemake DAG construction for native FASTQ and BAM/CRAM, one/multiple
   readgroups, male/WGS and female/WES FASTQ cases and ERF-corrected BAM.
   No obsolete rule producer is selected; native and external adapter routes
   are mutually exclusive. One-readgroup samples retain the existing merge
   bypass; multiple readgroups retain `merge_rgs`.
5. Compared complete DAG output/producer lists to the baseline for FASTQ
   (one/two readgroups), BAM and CRAM through representative **QC, DeepVariant
   and chrM tail targets**: identical. Caller-rule loading also tested with
   DeepVariant, HaplotypeCaller and BOTH, and chrM configuration choices.
6. Executed a real native FASTQ worker invocation with root-workflow reparsing
   and `--allowed-rules adapter_removal --mode subprocess`.
7. Compared original inline native adapters with the new shared implementation
   using actual AdapterRemoval 2.3.3 and pigz: Phred33, Phred64, duplicate
   removal, and retry rescue with out-of-order pairs. Decoded trimmed FASTQs,
   FASTQ statistics, adapter-identification output, normalized trimming metrics
   and rescue `.errors` entries match.
8. Compared original/new external fused runners using tiny real **BAM and
   CRAM** inputs with samtools 1.17: extracted raw FASTQs, singletons, trimmed
   FASTQs, FASTQ statistics and adapter reports match; scratch is removed.
   These direct runner tests disable leases; protocol behavior is covered by
   the separate existing lease/mock-runner tests, not a live scheduler job.
9. Built `bam_merge`, `bam_dechimer`, `fix_bam_rg_pairs`, `fastcheck` and
   `fastcheck_hts` in this checkout; verified extension imports and dynamic
   library resolution against the preprocessing environment.
10. Kept the restart recovery/frozen-manifest/DAG/worker tests and added a
    completed-native-FASTQ constraint test. Finished samples cannot re-enter
    either adapter route; cohort input behavior is unchanged.

An initial real baseline comparison exposed a build prerequisite: automatic
`fastcheck` compilation printed messages into the streaming shell pipeline,
causing the *old* inline adapter run to fail. After prebuilding the extension,
the comparisons passed. Always build and validate tools before workflow jobs;
do not depend on an on-demand build inside a FASTQ stream.

## Reproduce the checks on this installation

The environment paths below describe this Snellius validation, not a portable
Spider installation recipe. The full integration comparisons need a frozen
copy of the baseline and access to the configured reference bundle. Without
`FUSED_CLEANUP_BASELINE`, comparison tests are skipped; without the reference
bundle, site-specific workflow tests are skipped. Core unit/contract tests
remain independent of those integration prerequisites.

```bash
cd /gpfs/home2/hulsmanm/projects/short_read_pipeline_fused_only_20260907
PIPELINE_TEST_PYTHON=/home/hulsmanm/all_data/research/conda/marc/miniconda3/envs/snakemake/bin/python
PIPELINE_PREPROCESS=/home/hulsmanm/.snakemake/d6b82cba48e3c9dc9255474984326c51_
PIPELINE_BASELINE=$(mktemp -d /tmp/fused-cleanup-baseline.XXXXXXXX)
git archive c02883b | tar -x -C "$PIPELINE_BASELINE"

for PIPELINE_BUILD in "$PWD" "$PIPELINE_BASELINE"; do
    HTSLIB_INCLUDE="$PIPELINE_PREPROCESS/include" \
    HTSLIB_LIBDIR="$PIPELINE_PREPROCESS/lib" \
    make -C "$PIPELINE_BUILD/scripts" all build-hts \
        PYTHON="$PIPELINE_PREPROCESS/bin/python" \
        HTSLIB_PREFIX="$PIPELINE_PREPROCESS"
done

PATH="$PIPELINE_PREPROCESS/bin:$PATH" \
FUSED_CLEANUP_BASELINE="$PIPELINE_BASELINE" \
"$PIPELINE_TEST_PYTHON" -m pytest -q
git diff --check
```

The compiled binaries in this checkout are local build artifacts, not committed
files. Rebuild after cloning elsewhere or changing Python/htslib environments.

## Limits and next use

These checks are not a full production-size sequencing run or biological
equivalence benchmark for every optional endpoint. DeepVariant, DRAGMAP,
GATK/chrM and BAM-QC algorithms were not replaced; their existing runner tests
mostly use fake tools. Real adapter comparisons use small synthetic inputs.
No Spider pilot, new cluster deployment or live ZSlurm submission was tested.

Use this branch for a separately scoped canary/new run. Do not switch the
source path of an already-running controller/worker population mid-run.
Before restarting an existing workdir, follow `RESTARTING.md`, use the frozen
manifest and inspect the intended DAG. Do not change workflow input/output
paths or staging manifests merely to point at this checkout.
