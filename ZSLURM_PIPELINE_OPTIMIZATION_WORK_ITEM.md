# Work item: short-read pipeline resource and I/O redesign

Historical design/progress notes. The current fused-only implementation and
deployment guidance are in [README.md](README.md) and
[PORTABILITY_PLAN.md](PORTABILITY_PLAN.md). Retired rule names, sort pilots and
`fuse_*` switches below describe earlier revisions, not current run options.

Status: **implementation in progress; chrM-tail and BAM-QC pilots complete**
Last updated: **2026-08-10**
Scope: `short_read_analyzing_pipeline_Snakemake`, ZSlurm, and the native
Snakemake ZSlurm executor

## Goal

Reduce shared/network-storage load and avoid CPU, memory, SSD, active-space, and
dCache-transfer deadlocks while preserving restartability and correct output.
Support several Snakemake pipelines concurrently, with explicit priority for
the run that should make progress first.

Do not implement executable pipeline changes while the current production run
is active. Documentation-only changes are safe.

## Current baseline

### Deployment

- ZSlurm source: `~/projects/cluster_manager`, with priority, dynamic leases,
  SSD accounting, and directional dCache pools implemented.
- Native executor: `~/repos/snakemake_executor`, with pipeline priority and
  directional dCache resources implemented.
- Persistent manager config contains:

  ```yaml
  dcache_download_slots: 4
  dcache_upload_slots: 4
  ```

- The manager started on 2026-07-26 still reports health schema 1. It therefore
  runs the old in-memory manager behavior until an authorized restart.
- The current pipeline worktree is heavily modified. Preserve all unrelated
  edits and inspect diffs before touching a rule.

### Measured `split_alignments_by_readgroup`

From `report-2026-07-26_13-36.tsv`, filtered to successful `projectmine`
(`mine_*`) jobs:

| Statistic | Runtime |
|---|---:|
| Successful jobs | 2,043 |
| Mean | 2:33:16 |
| Median | 2:32:57 |
| p90 | 2:56:11 |
| p95 | 3:05:39 |
| Maximum | 4:15:49 |
| Minimum | 1:42:14 |

The current configured runtime estimate is 8,700 seconds (2:25), which is below
the observed mean and well below p95. Re-estimate this rule before the next run;
a p95-based request plus scheduling buffer is safer than the current value.

The current report contains 12 failed `projectmine` split attempts. Diagnose
their causes separately before changing retry or memory behavior.

## Invariants

1. Every sample must reserve its active-space lifecycle footprint, including
   samples whose input already exists on active disk.
2. A job may not reserve destination space and then wait indefinitely for a
   separate transfer job that cannot dispatch because the destination is full.
3. A fused multi-phase job initially requests the maximum CPU and memory needed
   by any phase.
4. Resource changes inside a job are absolute and idempotent.
5. A later high-resource phase may start only after a successful lease
   reacquisition.
6. True storage offload means `/scratch-node`; `/scratch-local` is still GPFS.
7. `ssd_gb` represents peak simultaneous local bytes, not final output size.
8. Download and upload concurrency are separate from durable dCache GB.
9. Pipeline priority never preempts running work; it orders eligible waiting
   work.
10. Snakemake continues to own the DAG and retry semantics.

## Workstream 1: deploy and verify current ZSlurm changes

Perform this only after the active pipeline has drained and restart has been
authorized.

1. Record before-restart state:
   - instance name;
   - all active/assigned jobs;
   - engines and Slurm allocation IDs;
   - durable active/dCache/archive totals and in-use values;
   - legacy transfer total/in-use/pending;
   - relevant config;
   - manager and node reports.
2. Stop/drain the old manager without losing active jobs.
3. Start a new manager from the verified current source/install.
4. Start new chiefs; existing chiefs cannot gain lease support in place.
5. Verify manager health capabilities:

   ```text
   schema_version: 2
   submit_job_priority: true
   list_jobs_priority: true
   directional_dcache_slots: true
   ```

6. Verify `zsqueue` shows `PRIO` and JSON priority is numeric.
7. Verify raw status shows independent download and upload pools.
8. Verify durable accounting after restart before allowing new submissions.
9. Do not use `recompute-inuse` automatically: it can discard deliberately
   persistent producer accounting.

## Workstream 2: combine sample start and source routing

Replace the current three-rule control flow:

```text
start_sample
  -> archive_to_active
  -> dcache_to_active
```

with one simple `start_sample` job that selects and executes one route:

```text
input already active
  -> validate/claim, no copy

input in archive
  -> copy and, where required, convert to active

input in dCache
  -> checksum-verified download to active
```

Requirements:

- Take the same active lifecycle reservation for all three routes.
- The active route still reserves space even though it performs no copy.
- The archive route carries archive accounting appropriate to the bytes it
  stages/releases.
- Only the dCache route requests `dcache_download_slots=1`.
- Do not reserve a download slot while doing non-transfer preprocessing.
- Write the existing completion/indicator output only after validation succeeds.
- Make retries idempotent: recognize and validate a complete destination;
  clean or quarantine incomplete destinations.
- Preserve BAM, CRAM, and FASTQ input routing.
- Preserve CRAM reference and archive conversion requirements.

This removes the deadlock in which many `start_sample` jobs consume active
capacity while separate archive/dCache jobs needed to fill or eventually
release that capacity cannot start.

## Workstream 3: concrete phase combination and dynamic leases

### Keep `split_alignments_by_readgroup` before alignment

- BAM/CRAM input with multiple read groups continues through
  `split_alignments_by_readgroup`.
- FASTQ input bypasses that split.
- The split remains a separate job initially because it produces independently
  schedulable read-group inputs and currently takes roughly 2.5–3 hours.
- Raise its time estimate from 8,700 seconds using the measured distribution.

### First fusion candidate: per-read-group alignment chain

Evaluate and, if output/retry semantics remain sound, combine the linear
per-read-group segment:

```text
align_reads
  -> merge_bam_alignment
  -> dechimer
  -> check/prepare the aligned read-group output
```

Resource order:

1. run the highly threaded alignment phase at the job's maximum reservation;
2. perform any other high-thread work before releasing cores;
3. call `zslurm_lease set` to shrink for merge/dechimer/check phases;
4. finish without growing again.

Do not place a multithreaded sort after the shrink point. Choose one:

- execute the sort before releasing cores;
- keep `sort_bam_alignment` as a separate downstream job; or
- explicitly reacquire before sort and handle timeout, although this is less
  desirable than a monotonic high-to-low resource profile.

The initial implementation should prefer monotonic high-to-low behavior.

Example control flow:

```bash
run_alignment_with_24_threads
run_any_remaining_high_thread_step
zslurm_lease set --cores 2 --mem-gb <measured-low-phase-peak>
run_merge
run_dechimer
run_checks
```

The tool's own thread arguments must change with the lease; a lease alone does
not change threads or cgroups.

### Reacquisition behavior

When growth is genuinely required:

```bash
zslurm_lease set --cores 24 --mem-gb 48 --wait 3600
```

- The request waits FIFO.
- The local chief stops admitting new child jobs while the first acquire waits.
- Existing jobs are not preempted.
- CPU and memory are granted atomically.
- Timeout leaves the previous lease unchanged.
- Never start the high-resource tool before the acquire succeeds.

### Alternative to fusion

If one high-thread producer feeds several independent low-thread jobs, keep the
producer separate and let it complete/release normally; then let Snakemake
schedule the low-thread children in parallel. This retains better restart
granularity.

## Workstream 4: move repeated shared-storage I/O to node-local SSD

### Candidate pattern

For a large BAM/CRAM/FASTQ that is searched, scanned, sorted, or merged several
times:

1. request `ssd_use="required"`;
2. reserve a measured `ssd_gb`;
3. copy or stream one verified copy to `$TMPDIR` on `/scratch-node`;
4. perform repeated reads and intermediates locally;
5. copy only required final artifacts back;
6. clean local data on success and failure.

Do not describe `/scratch-local` as SSD offload: it remains GPFS.

### SSD estimate

For each rule or fused job, measure:

```text
peak local bytes =
    copied inputs
  + expanded/decompressed/converted copies
  + sort/merge/dechimer/spill temporaries
  + indexes
  + outputs waiting to be copied back
  + safety headroom
```

Use input-size-dependent functions with:

- multiplicative factor;
- fixed overhead;
- minimum reservation;
- upward rounding;
- explicit safety margin.

`ssd_use="required", ssd_gb=0` is not acceptable for storage-heavy work because
it constrains placement without preventing SSD overfill.

### Current rules to re-measure

The current dirty pipeline contains SSD reservations in alignment, GLnexus,
chrM, statistics, DeepVariant, and database-import rules. Treat the present
numbers as hypotheses. For each candidate log:

- `df` and `du` at phase boundaries;
- local input, temp, and output sizes;
- peak node SSD used and reserved;
- shared-filesystem read/write bytes;
- CPU use and I/O wait;
- sample input size/type.

Validate using large/worst-case samples, not only the median.

## Workstream 5: directional dCache transfers

### Manager

Persistent defaults:

```yaml
dcache_download_slots: 4
dcache_upload_slots: 4
```

Runtime control:

```bash
zscontrol transfer-slots --download N --upload N
```

TUI:

```text
DCache xfer D/U: download_inuse/download_total |
                  upload_inuse/upload_total
```

- key `7`: download total;
- key `8`: upload total.

Lowering a limit does not revoke running reservations; it blocks new matching
work until the pool drains.

### Pipeline/plugin

- Inbound transfer rules use `dcache_download_slots=1`.
- Outbound transfer rules use `dcache_upload_slots=1`.
- Never combine either with explicit legacy `dcache_transfer_slots` in one
  rule.
- Do not put directional slots in Snakemake's global resource-capacity list;
  ZSlurm is the cross-pipeline arbiter.
- Remove or explicitly justify any additional pipeline-local advisory lock.

The current uncommitted pipeline contains two directional download sites and
fourteen directional upload sites. Re-audit every actual dCache transfer after
the run.

Rolling compatibility:

- the new plugin sends directional metadata plus a conservative legacy fallback;
- the old manager throttles the fallback through one combined pool;
- the new manager ignores fallback when direction is present.

## Workstream 6: pipeline priority

Use one integer priority for every job submitted by one Snakemake invocation:

```bash
snakemake --executor zslurm --zslurm-priority 200 ...
```

or:

```yaml
executor: zslurm
zslurm-priority: 200
```

Properties:

- default `100`;
- higher values run first;
- zero and negative values are valid;
- no preemption;
- packing and FIFO/LIFO apply only inside equal-priority bands;
- ineligible high-priority work does not block eligible lower-priority work;
- priority also orders staging/download/upload/final-write jobs.

Use well-separated operational bands, for example:

- `300`: urgent/front pipeline;
- `200`: next production pipeline;
- `100`: normal/default;
- `0` or negative: background.

Verify the running manager supports priority before using a non-default value.

## Workstream 7: observability and resource retuning

Add phase-level instrumentation for candidate fused/SSD jobs:

- timestamps and phase name;
- requested/current/max lease;
- application thread count;
- process-tree PSS/RSS;
- local `df`/`du`;
- shared and local read/write counters where available;
- input/output size;
- return code and cleanup state.

After a representative run, aggregate:

- runtime mean/median/p90/p95/max by rule and input-size bin;
- peak memory vs reservation;
- used vs reserved core-hours;
- peak SSD vs reservation;
- shared-storage bytes and I/O wait;
- dCache transfer throughput and failure rate.

Tune from p90/p95 and worst-case samples. Do not set production resources from
one average alone.

## Implementation order

1. Finish/drain the current run.
2. Preserve state and deploy/restart the new ZSlurm manager and chiefs.
3. Verify priority, leases, directional transfer pools, schemas, and accounting.
4. Commit/test the current plugin and directional pipeline changes separately.
5. Implement the single-rule `start_sample` routing.
6. Add instrumentation before changing SSD formulas further.
7. Pilot local-SSD staging on one high-I/O rule.
8. Pilot the per-read-group alignment fusion with monotonic high-to-low leases.
9. Re-measure on a small BAM/CRAM/FASTQ test matrix.
10. Dry-run and execute a small cohort.
11. Compare integrity, runtime, network I/O, SSD peak, CPU, and memory.
12. Roll out to a complete cohort only after acceptance criteria pass.

## Validation matrix

Test at least:

- input already active;
- input in archive;
- input in dCache;
- BAM with one read group;
- BAM with multiple read groups;
- CRAM with required reference;
- paired FASTQ;
- a small sample and worst-case large sample;
- SSD available and SSD temporarily full;
- download pool full while upload has capacity;
- upload pool full while download has capacity;
- two pipelines with different priorities;
- lease shrink;
- lease reacquire immediate, delayed, and timed out;
- manager/client rolling compatibility.

## Acceptance criteria

- No start/transfer storage deadlock.
- Already-active inputs retain an active reservation.
- Output/checksum equivalence with the current pipeline.
- No high-resource phase starts without its lease.
- Released CPU/memory is reused without starving reacquisition.
- No SSD overfill in the validation matrix.
- Measured shared-storage traffic decreases for selected fused/localized paths.
- Download saturation does not block upload-only work and vice versa.
- The higher-priority pipeline wins eligible waiting work.
- Old-manager rolling fallback remains conservative until restart.
- Restart/accounting procedure is documented and reproducible.
- Rule time, CPU, memory, and SSD requests reflect measured distributions.

## Operational note for the 2026-07-31 manager allocation

The manager/local engine is inside Slurm job `24920582` on `fcn76`, ending at
2026-07-31 13:32:46 local time. ZSlurm reports an incorrect generic 120-hour
local-engine timeleft and can therefore admit work that will outlive the real
allocation.

At the measurement point three `split_alignments_by_readgroup` jobs were
running. Their historical distribution indicates they are likely to finish
inside the remaining allocation.

Recommended handover:

1. mark `fcn76` phase-out in the TUI (`o`, then `fcn76`);
2. let the three current split jobs finish;
3. confirm no new long job was backfilled;
4. switch/restart on the new allocation before 13:32;
5. verify manager, queue, budgets, and outputs after handover.

## Implementation update: measured SSD sort pilot (2026-08-09)

The first isolated implementation is the existing `sort_bam_alignment` SSD
candidate. It was selected because it can be tested independently without
changing sample routing, source files, or completed production outputs.

Implemented:

- `scripts/run_samtools_sort_ssd.py` now creates the samtools-sort temporary
  directory below the exact `/scratch-node/<user>.<SLURM_JOB_ID>` allocation.
- The runner fails when `ssd_use="required"` has no writable assigned node SSD;
  it no longer borrows the newest unrelated `/scratch-node/<user>.*` directory
  or silently treats GPFS `/scratch-local` as SSD.
- `scripts/io_profile.py` records JSON measurements for the process tree and
  local work directory: timestamps, input/output sizes, requested SSD/threads/
  memory, RSS/PSS, process I/O counters, local peak bytes, filesystem free
  space, return code, signal state, and cleanup result.
- `sort_bam_alignment` writes these measurements to
  `logs/Aligner/<sample>.<readgroup>.sort_bam_alignment.io.json` without adding
  a new pipeline completion target.
- Success and command failure both remove the per-job scratch subdirectory.

Validation performed:

- Python syntax checks passed.
- The focused profiler/runner tests passed: 4/4.
- The complete repository test suite passed: 28/28.
- The production `Snakefile` parsed successfully in a dry-run against the
  already-complete `pipeline.done` target.
- A separate native-ZSlurm dry-run resolved exactly two executable pilot jobs;
  only `sort_on_node_ssd` requested `ssd_use=required, ssd_gb=1`.
- The isolated ZSlurm execution completed as jobs `671622`
  (`make_unsorted_bam`) and `671623` (`sort_on_node_ssd`) on `fcn77`.
- The profiler confirmed the real assigned path
  `/scratch-node/hulsmanm.25374087/...`, return code zero, and successful
  removal of the job scratch directory.
- `samtools quickcheck` and the generated index succeeded. The two deliberately
  unsorted fixture reads were emitted in coordinate order (positions 101 and
  501).

The tiny correctness pilot fits entirely in memory, so its measured sort-spill
peak is zero; it validates placement, measurement, output integrity, and
cleanup but is not evidence for changing the production `ssd_gb` formula. The
synthetic failure test did exercise non-zero local usage and verified cleanup.
A representative real read-group BAM is still required before retuning the
current factor.

The running manager/chiefs were deliberately not restarted: unrelated jobs
are still active and the manager remains the old in-memory deployment.
Production activation of lease-dependent fusion therefore remains pending
until the controlled manager/chief handover described above.

## Implementation update: active-storage fusions (2026-09-15)

Implemented in the isolated `codex/active-storage-fusion-20260915` branch:

- `markdup` now consumes all validated readgroup BAMs. For multiple readgroups,
  `samtools merge` writes only to assigned node SSD; `samtools markdup` consumes
  that local BAM and atomically publishes only the final BAM, index and stats.
- The job starts at the former merge CPU reservation and monotonically shrinks
  to the former markdup reservation. A single-readgroup sample starts directly
  at the markdup reservation. Memory never grows between phases.
- `cram_encrypt_fused` writes CRAM/CRAI to assigned SSD, encrypts there and
  uploads both products directly from SSD. Only the checksum and `.copied`
  receipt cross back to active storage.
- Over 19,000 historical jobs, upload occupied 3.8% of combined CRAM,
  encryption and upload time. The entire job therefore holds a weighted
  `dcache_upload_slots=0.05` reservation: the default pool of four admits 80
  fused jobs, while existing pure transfers retain their full-slot limit.
- The CRAM job monotonically shrinks from conversion to the combined
  encryption/upload tail. It never waits inside a worker for a phase-time
  global slot acquisition.
- Initial SSD reservations include simultaneous local input/intermediate/final
  bytes plus fixed headroom. Runner metrics record per-phase scratch high-water
  use for production calibration; these first estimates must not be reduced
  from tiny-fixture measurements.

The normal durable filenames are unchanged: `bams/{sample}.markdup.bam`, its
index/statistic, `cram/{sample}.mapped_hg38.cram.ADLER32`, and
`cram/{sample}.mapped_hg38.cram.copied`. The old GPFS-only merged BAM,
plaintext CRAM, encrypted CRAM and CRAI are no longer DAG products. Removing
the former transfer payload also removes 15% from the initial active-storage
reservation; the remaining add/remove accounting balances at markdup and
`finished_sample`.

Focused validation covers one/multiple readgroups, no-dedup, duplicate,
supplementary, secondary and unmapped records; old/new samtools record and
normalized markdup-stat equivalence; real CRAM 3.1 plus Crypt4GH round-trip;
monotonic lease calls, failure publication boundaries and scratch cleanup.
Existing direct/legacy upload tests continue to cover the unchanged transfer
job. Full-suite and measured canary results are recorded with the
branch before deployment.

## Implementation update: routed start and alignment fusion (2026-08-09)

Implemented, but deliberately feature-gated and not submitted to production:

- `start_sample` now owns the active lifecycle claim and exactly one active,
  archive, or dCache materialization route. The old `archive_to_active` and
  `dcache_to_active` jobs are removed.
- All routes publish one atomic `source/<sample>.route_ready` contract only
  after validation. External destinations use manifests, partial directories,
  atomic promotion, and recoverable quarantine of incomplete prior attempts.
- The dCache route alone requests `dcache_download_slots=1`; active/archive
  routes do not consume a transfer slot. Legacy route markers remain readable
  for one-time adoption.
- Active input bytes are included in lifecycle accounting for every route,
  including already-active samples.
- `fuse_alignment_phases=true` selects one `align_reads_fused` producer for
  DRAGMAP, merge, conditional dechimer, and both checks. The aligned BAM and
  merge/dechimer intermediates remain on the assigned `/scratch-node` and only
  terminal outputs are atomically copied back.
- The fused job starts at 22.75 cores/40 GB, preflights the lease before
  DRAGMAP, then issues the absolute/idempotent target of 6 cores and 8 GB (WGS)
  or 6.5 GB (exome). Low-phase tool threads are reduced to one per component.
- Coordinate sort remains a separate downstream job; the first fusion is
  monotonic high-to-low and never needs to reacquire.
- Phase JSON records lease responses, process-tree memory/I/O, local peak,
  return codes, requested resources, and scratch cleanup. The first SSD factor
  is explicitly provisional until a representative sample is measured.

Validation performed:

- Python syntax validation passed for the routed and fused runners.
- A temporary routed-materialization test passed through manifest validation,
  atomic promotion, and completion markers.
- A fake-tool end-to-end fused test passed through alignment, lease `status`,
  absolute lease `set`, merge/check, conditional dechimer/check, atomic output
  publication, phase metrics, and scratch cleanup.
- A targeted production-DAG dry-run for one existing sample/read group passed
  in both modes. Legacy mode selected 13 jobs including separate `align_reads`
  and merge; fused mode selected 12 jobs with only `align_reads_fused`.
- That dry-run exposed and fixed a pre-existing `getsize()` failure for a
  checkpoint-generated CRAM. The resource function now uses a conservative
  full-sample upper bound until the split file exists.

Rollout gate:

1. Keep `fuse_alignment_phases=false` while health schema 1 chiefs are active.
2. After the controlled manager/chief handover, verify lease status and a tiny
   real ZSlurm child before enabling `alignment_lease_mode=required`.
3. Run a single small real read group, compare BAM/check/stat outputs with the
   legacy chain, and inspect the phase metrics before a small-cohort pilot.

## Implementation update: extended opt-in fusions (2026-08-10)

This section supersedes the 2026-08-09 alignment-tail details above. The
legacy rules remain available and every new producer is disabled by default.
No production job was submitted during this implementation.

Implemented:

- `fuse_alignment_phases=true` now keeps coordinate sort in the same
  per-read-group job. The final sorted BAM/index plus all legacy stats/check
  contracts are published atomically. The monotonic low lease is 6 cores and
  15 GB, which covers the measured sort-memory tail without reacquisition.
- `fuse_kmer_sex=true` keeps the KMC database on assigned SSD and shrinks from
  2 cores/36 GB to 0.5 core/3 GB for sex validation.
- `fuse_deepvariant_phasing=true` keeps raw regional DeepVariant VCF/gVCF
  output on assigned SSD and continues through the existing Whatshap, merge,
  statistics, and exome-extraction output contracts. It shrinks from
  8 cores/10 GB to 1 core/9 GB after DeepVariant.
- `fuse_external_adapter=true` extracts each BAM/CRAM read group once and
  performs adapter processing against the SSD-local FASTQs. It also publishes
  the legacy temporary raw FASTQs because the later BAM-tag merge consumes
  them; this prevents the duplicate extraction exposed by the combined DAG.
- `fuse_chrm_extract_align=true` combines both chrM/NUMT extraction branches
  and four BWA alignments. Four paired FASTQ intermediates and alignment
  scratch stay on assigned SSD; the eight legacy BAM/index paths are unchanged.
- Each runner fails closed when `ssd_use=required` has no assigned writable
  scratch directory, removes job-local scratch, records phase I/O metrics, and
  publishes final files through atomic rename.

Validation performed:

- All 35 repository tests pass, including fake-tool end-to-end tests for all
  five fused runners.
- Python syntax checks pass for all five runners and `git diff --check` is
  clean.
- A forced combined Snakemake dry-run for one real three-read-group sample
  completed through regional DeepVariant and the final annotated chrM gVCF.
  Its 17-job DAG contains three `external_adapter_fused`, one
  `kmer_sex_fused`, three `align_reads_fused`, one
  `deepvariant_phasing_fused`, and one `chrm_extract_align_fused` job. It
  contains none of their replaced legacy extraction/adapter/KMC/alignment,
  sort, phasing, or chrM extraction/alignment rules.
- The combined dry-run found and fixed two composition defects before any
  execution: an overlapping alignment provenance-output ambiguity and a
  duplicate BAM/CRAM-to-FASTQ extraction.

Rollout gate and next action:

1. Keep every fusion flag false on the currently running schema-1 manager.
2. Deploy/restart the schema-2 manager and chiefs, then verify one native
   `zslurm_lease status`/absolute `set` smoke job and accounting.
3. Enable one fusion at a time for one small sample, compare output checksums,
   stats, peak SSD, RSS, wall time, and lease history with the legacy rules.
4. Only after those comparisons, enable the combined five-flag pilot on a
   tiny mixed BAM/CRAM cohort.
5. Treat the remaining Mutect/merge/BP chrM tail and small QC/report chains as
   a separate next tranche; their broader output/restart surface is not yet
   fused.

## Implementation update: chrM tail and BAM-QC fusion (2026-08-10)

Implemented behind disabled-by-default feature flags:

- 'fuse_chrm_mutect_tail=true' keeps the four standard Mutect calls, shifted
  liftover, Mutect-stat/VCF merges, filters, four BP-resolution calls,
  normalization, and FILTER annotation in one assigned-SSD job. Only the two
  final annotated gVCFs and indexes are published.
- 'fuse_bam_qc=true' stages the markdup BAM/index once and runs VerifyBamID2,
  HsMetrics, sequencing-artifact/OxoG metrics, both samtools-stat scans, both
  sampled bamstats scans, and mosdepth as six independent parallel task
  groups. It atomically publishes all 18 existing QC output contracts.
- The chrM tail keeps a fixed resource envelope. BAM-QC uses a 12-core/12-GB
  initial lease; each of its six parallel task groups atomically releases its
  relative share at completion, including on task failure. Stable per-task
  release ids make a repeated response harmless. `bam_qc_lease_mode=required`
  therefore needs the upgraded manager/chief; use `disabled` only for local
  tests.

Validation performed:

- All 38 repository tests pass. The BAM-QC tests observe exactly six unique
  relative release calls totalling 12 cores and 10,496 MB, and verify that a
  failing task releases its share before returning failure.
- All 40 tracked ZSlurm tests pass, including atomic concurrent relative
  releases, release-id idempotence, the observed-memory floor, and a real
  `zslurm_lease release` CLI/Unix-socket round trip.
- An isolated one-sample Snakemake dry-run selects exactly one
  `bam_qc_fused` job, no legacy QC producer, and resolves the runner with
  12 cores/12 GB, required SSD, and `--lease-mode required`.
- A combined one-sample dry-run selects all seven fused producers and no
  overlapping legacy producers. It exposed and fixed one remaining
  'align_reads' versus 'align_reads_fused' provenance-output ambiguity before
  execution.
- ZSlurm job '927177' ran the fused chrM tail on the existing small chrM/NUMT
  BAMs in 247 seconds. All 16,291 chrM and 16,285 NUMT output records are
  identical to the existing legacy outputs; assigned scratch was removed.
- ZSlurm job '927179' ran all fused BAM-QC consumers on an isolated 144 MB
  chr21 slice. The runner published all 18 outputs with 'success=true', removed
  assigned scratch, and peaked at 6.28 GB RSS. Both bamstats outputs were
  byte-identical to independently recomputed outputs; both samtools-stat
  outputs were identical after excluding only the expected command-line header
  containing the SSD-local BAM path.
- The pilot used a regional dbSNP slice because the production dbSNP contains
  1.1 billion records; the first attempted small pilot was cancelled when its
  full dbSNP scan made clear that it was not a small test. Production retains
  the full configured dbSNP input.
