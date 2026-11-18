# Dependency inventory and Spider/ZSlurm migration plan

Date: 2026-09-07; implementation update: 2026-09-17. Steps 1 and 2 and the code
portion of step 3 are published on the default pipeline and ZSlurm branches.
They are not automatically deployed into running controllers/managers. A real
Spider preprocessing canary passed, but the workflow is still **not
end-to-end tested on Spider**. Running pipelines and schedulers are not
reconfigured by the repository merge.

## Conclusions

1. Production DeepVariant defaults to **native execution**, not Apptainer.
   `Snakefile` defaults to `caller=Deepvariant`; `deepvariant_phasing_fused`
   uses `get_deepvariant_native_runner`. The Apptainer comparison target has
   its own output namespace and does not feed production phasing. A profile
   with `use-apptainer: true` does not change that selection. GLnexus joint
   calling is a separate container dependency.
2. **Cloning the pipeline and copying the resource folder is not sufficient
   today.** Software installation/builds, relocatable path configuration,
   scheduler + executor installation, input manifests and authorized storage
   access are also required. Some important paths are outside the resource
   folder. The remote default branch now contains the cleanup/portability
   implementation, but still requires the separately pinned components and
   site installation described below.
3. Keep one workflow implementation, with site configuration for Snellius and
   Spider. Do not fork the biological algorithms, increase reservations, or
   force all jobs through Apptainer as part of this migration.

## What belongs in an installation

| Component | Current dependency/evidence | Copying resources sufficient? |
|---|---|---|
| Workflow + restart handling | This repository, including `restart_state.py`, `scripts/prepare_restart.py`, runners and delivery templates | Clone the tested revision, including the cleanup/restart commits |
| Controller/worker Snakemake | Local fork reports `9.13.8.dev16`; existing installations at `83c79c1d` passed the Snellius canary, while the component lock pins `50bd8694` with duplicate-rule validation restored | No; retain the fork's checkpoint-temp, Conda-`run:`, metadata-cleanup and Apptainer-symlink fixes. Stock Snakemake was a compatibility test, not a replacement |
| ZSlurm manager/chiefs | Separate `holstegelab/zslurm`, release commit `3b88987`; own `env.yaml`/installation | No; install into the environment used by manager and pilots; repository merge does not update a running manager process |
| Native executor | Separate package `snakemake-executor-plugin-zslurm`; release fix commit `9ed793e` transports declared/resolved environment values to shell-free workers | No; install the pinned commit into the Snakemake environment on controller and workers |
| Controller Python modules | `Snakefile`, `common.py`, `read_stats.py`: pandas, numpy, PyYAML, h5py, Snakemake and its plugins | No; a complete pinned controller environment is not supplied by this pipeline |
| Rule environments | `envs/*.yaml` and **post-deploy scripts**; environments currently live outside resources under the user's `.snakemake` prefix | No; recreate at the final installation prefix, then smoke-test |
| Custom native tools | `scripts/Makefile`: `bam_merge`, `bam_dechimer`, `fix_bam_rg_pairs`, `fastcheck` and `fastcheck_hts` extensions | No; binaries/extensions are Git-ignored. The preprocessing post-deploy hook now builds them automatically on target against that environment's Python/htslib ABI |
| Native DeepVariant 1.9.0 | `software/deepvariant-1.9.0-native` under resources; `scripts/prepare_deepvariant_native.py` | Not by moving the installed runtime: launchers embed absolute paths. Recreate at final prefix from pinned image/archive and verify |
| Containers | GLnexus v1.4.1 for joint calling; optional DeepVariant comparison image 1.9.0 | Not normally: images/cache are in a separate `.apptainer` prefix. Provision Apptainer and cache verified images if those endpoints are selected |
| References/databases | See resource groups below, including all indexes and region files | Usually data can be copied, but verify completeness, symlink targets, hashes and licenses |
| Input data | Sample `.tsv`, `.source`, optional `.target`/`.exclude`, FASTQ/BAM/CRAM or authorized dCache/S3 objects | No; not part of the reference resource bundle |
| Existing-run state | Frozen restart manifest, completion receipts, retained gVCFs/QC/sex files, statistics archives and workdir | No; needed only to continue an existing run. A fresh run should not inherit old `.adat` caches containing absolute paths |
| Transfers/authentication | `dcache_cp` + `rclone`, ADA/staging helpers; AWS CLI for S3; optional iRODS tools/credentials | No; install clients and provision project/user credentials separately |
| Encryption | Crypt4GH in preprocessing env; sender private key and intended recipient public keys (`Encrypt.smk`) | No automatic sharing: provision approved keys separately |
| Host utilities/build tools | Bash/process substitution, coreutils, pigz/gzip/bzip2, tar, git, make, GCC/G++, zlib/htslib development libraries; wget/unzip/skopeo for relevant deployments | Verify on target; not all are covered by each environment YAML |

The revision checks above identify the inspected installations, not a claim
that these repositories are clean, pushed, or fully reproducible. Release
preparation must record dirty diffs and package lockfiles as well as HEADs.

### Resource bundle groups

- Alignment/reference: `hg38_phix/` FASTA, FAI, dictionaries, male/female
  DRAGMAP indexes and masks; `cram_refs/` for decoding source CRAMs.
- Calling/QC: `intervals/`, WGS/WES regional bins and padded intervals,
  capture-kit BED/interval lists, precomputed capture-kit JSON, known sites,
  dbSNP, VerifyBamID SVD data, windows, mappability tracks and ploidy priors.
- Adapters/sex/classification: adapter list, four k-mer reference sets,
  Kraken/Bracken database with matching database files.
- chrM: original and shifted mitochondrial references/indexes, NUMT regions,
  chain/reference preparation products.
- Optional endpoints: GLnexus custom preset, annotation databases
  (gnomAD/ClinVar/REVEL and licensed ANNOVAR), patched bcftools/Java/GATK and
  other software selected by HaplotypeCaller/CNV/annotation workflows.

This is an inventory of code dependencies, **not a full recursive checksum or
completeness audit of the resource tree**. The inspected male/female reference
symlinks are relative and remain valid when their parent tree is copied; this
does not certify every symlink in every database/software subtree.

Do **not** distribute the resource folder blindly: `.agh/` and `.c4gh/` are
credential/key locations referenced by the workflow. Separate public/reference
data, software and secrets. Never put private keys, macaroons, AWS credentials
or cached signed URLs in the repository or a general-purpose resource archive.
No credential contents were read for this inventory.

## Concrete portability gaps in pipeline code

| Location | Original coupling | Release status / remaining change |
|---|---|---|
| `constants.py: RESOURCES` | Hardcoded `/gpfs/work3/0/qtholstg/hg38_res_v2/`; most paths were derived at import time | **Implemented:** schema-1 `site_config.py` is selected before importing `common`/`constants`; its absolute filename is declared to workers |
| `Snakefile.paths.yaml` | Historical/stale path file; not the general configuration source for the main pipeline. `GLnexus_to_iRODS.smk` reads some metadata from it | Replace/supersede with an actually loaded site configuration and consistent provenance; changing this YAML alone is insufficient |
| `Aligner.smk` / external adapter defaults | Corrected hg19/b37 chrY reference was outside the bundle at `.../marc/genome/hg19_b37chrY.fa` | **Implemented:** configured path plus checksum-pinned FASTA/FAI supplement; header-driven selection is retained |
| `constants.py`, `precompute_capture_kits.py`, `decrypt_c4gh.py`, `GLnexus.smk` | Additional fixed resource/tool/preset paths; old user-specific token defaults | **Implemented for selected paths:** resource, preset, CRAM references, transfer executables/remotes and key filenames resolve through the shared site configuration |
| `scripts/prepare_deepvariant_native.py` | Loader, Python, library and executable paths baked into launchers | Rebuild runtime at destination prefix; record source image digest and verify `--version` plus tiny WGS/WES jobs |
| `profiles/zslurm/config.yaml` | User-specific Conda/Apptainer cache paths and `/scratch-node` bind | **Implemented:** `profiles/spider` uses project prefixes and no false bind. ZSlurm exports the exact child scratch path. Logical `partition=compute` stays unchanged |
| `pipeline_runtime.assigned_scratch` | Accepted only Snellius `/scratch-node` layouts | **Implemented:** consumes `ZSLURM_SCRATCH_DIR`, retains the legacy path and explicit `SLURM_TMPDIR`, and never guesses from generic cross-site `TMPDIR` |
| `common.node_ssd_base`, `DBImport.smk`, `Deepvariant_apptainer.smk`, level-2 tar script | Multiple scratch-selection implementations; some scanned other `/scratch-node/$USER.*` directories | **Implemented:** all route through the shared resolver; foreign-job scanning is removed |
| child scratch cleanup | Multiple jobs share one pilot allocation | **Implemented:** the chief owns and cleans only the unique `.zslurm/jobs/job-<zslurm-id>-<suffix>` attempt directory; rule-level cleanup stays below that private root |
| `constants.py` / `tmpdir` | Shared temporary directory is relative to run workdir; `/scratch-local` secondary fallback | Explicit durable shared-temp root on target project filesystem; node-temp paths only for data whose consumers are inside the same job |
| Archive input route | `/opt/dacommands/bin/daget`/`dals` and Snellius archive paths | Do not enable unchanged on Spider. Select dCache/S3 input route or implement/test a genuine site archive backend |
| Encryption / transfers | Private sender key fixed under resources; tokens/default remotes spread across modules | **Path separation implemented:** explicit protected filenames and tools/remotes. Endpoint capability checks remain part of the later preflight step |

### Environments are more than their YAML filenames

Do not replace the custom environments with generic packages during the port:

- `align_fused.post-deploy.sh` rebuilds pinned DRAGMAP with fixes, including
  the release-assertion/allocator and tokenizer issues. It expects system
  GCC/G++ and network/build tools.
- `kmc.post-deploy.sh` builds a pinned KMC revision with thread/I/O/ISA fixes.
  Merely installing `kmc` does not reproduce it.
- `vcf_handling.yaml` declares GATK 4.4, but its post-deploy script installs
  and selects **4.5.0.0**. Record the executed binary, not just the YAML value.
- Several Python/tool dependencies are not fully pinned. Generate lockfiles
  and capture the exact environment + post-deploy source hashes after tests.
- Build `fastcheck` for each Python ABI that uses it; verify loader/fallback
  behavior in preprocessing/alignment and the distinct QC Python/PyPy env.
  Do not copy `.so` files between Python 3.9/3.11/3.12 installations.
  The isolated baseline adapter test demonstrated that an on-demand extension
  build can print into the FASTQ stream; prebuilding is a correctness
  prerequisite, not just a startup-time optimization.

## Spider-specific findings

Spider documents shared **CephFS** home/project storage and SSD-local scratch.
Project `Data`, `Share` and `Software` directories are available on login and
worker nodes. Use project space for references, installed software and durable
workflow state. CephFS metadata-heavy workloads (including Conda trees) deserve
measurement; the documentation discourages recursive `du` for usage checks.
[SURF: storage on Spider](https://doc.spider.surfsara.nl/en/latest/Pages/storage_on_spider.html)

Inside a Spider Slurm job, use `$TMPDIR`: documentation describes `/tmp` bound
to `/scratch/slurm.<JOBID>` and cleaned after the allocation. Do not infer that
an arbitrary `/tmp` on a login node is usable SSD. Published partition names
include `normal`/`short`; scratch and memory grants depend on allocation.
The documentation's newer hardware table and older per-core examples differ,
so verify the actual grants with `scontrol`, environment/cgroup and site policy
before sizing pilots. Do not advertise the full physical SSD to every partial
pilot just because `df` shows it.
[SURF: compute on Spider](https://doc.spider.surfsara.nl/en/latest/Pages/compute_on_spider.html)

### ZSlurm implementation on main, preprocessing-certified only

The inspected `cluster_manager/config/sites/spider.yaml` and
`zslurm._cluster_defaults` already provide `normal`, SSD feature `ssd`, a
dynamic 2--30-core/8-GB-per-core pilot profile, and demand-backed compute
autogrow capped at 15 running plus queued pilots. Use this as a starting
template, **not as measured current capacity**. Temporarily disable autogrow or
cap it at one pilot during the site-owned production canary.

Implemented in the ZSlurm release:

1. Scratch discovery is site-gated: Spider accepts its pilot `$TMPDIR`, while
   Snellius cannot accidentally advertise ordinary `/tmp` as SSD.
2. Scratch feature detection uses the configured `ssd`/`scratch-node` feature;
   all 29 currently visible Spider `normal` nodes advertise `ssd`.
3. Every child attempt receives a unique mode-0700 directory via
   `ZSLURM_SCRATCH_DIR`; `TMPDIR`, `TMP`, `TEMP`, and `TEMPDIR` point there and
   controller paths are removed.
4. Spider partial pilots advertise 73 GiB per allocated core, additionally
   bounded by physical total/free space. A 30-core pilot therefore exposes at
   most 2,190 GiB, not the whole physical device.
5. The nonexistent `staging` partition and staging autogrow are disabled;
   archive sources remain unsupported until a real Spider backend exists.
6. GPFS RDMA monitoring is disabled in the Spider site configuration.

Real `short` pilots verified manager ↔ worker RPC, lease resizing, private
scratch creation/removal and clean engine shutdown. Still required before
production are full cgroup/failure/restart coverage, Apptainer visibility,
dCache publication and alignment/calling jobs using the installed Spider
resource tree.

The implemented contract has the chief validate scratch, export the allocation
root/capacity and give each attempt a unique
`job-<zslurm_id>-<suffix>` subdirectory through `ZSLURM_SCRATCH_DIR`.
The pipeline consumes this contract. Snellius and Spider only differ in how
the chief discovers the allocation root; pipeline phases need no site branch.

`ssd_use=required` must retain its meaning. Do not fix detection by declaring
every job SSD-optional or by accepting a writable generic `/tmp`. Optional
rules may use configured shared temporary storage. Publish outputs atomically
to shared storage before success; scratch must never hold the only copy of
products needed by a later rule or a restarted workflow.

## Staged implementation plan and acceptance criteria

### 1. Freeze a distributable baseline

Implementation status: **merged to the default branch**. See
`RELEASE_BASELINE.md`, `deployment/component-lock.yaml` and
`deployment/environment-inputs.sha256`. The executor environment fix is a
separate tested commit. Solver-complete Conda locks remain desirable because
some upstream YAML dependencies are not version-pinned.

Finish/review the isolated fused-only cleanup; retain restart contracts.
Record tested commits + any patches for pipeline, Snakemake, ZSlurm and plugin.
Provide controller and rule environment locks, build instructions and a
runtime/tool-version manifest. Keep this distinct from optimization changes.

Acceptance: a fresh clone builds all required custom tools and passes tests
without borrowing Git-ignored binaries from another checkout.

### 2. Separate resources, software and secrets

Implementation status: **implemented for the pipeline and release layout**.
`site_config.py` loads one validated YAML before path derivation and passes its
absolute filename to workers. Snellius and Spider examples are under
`config/sites/`; `DEPLOYMENT.md` documents the public resource archives and
separate private key provisioning. Destination installation and checksum
verification are still site-operator actions.

Generate an endpoint-specific resource manifest: logical name, relative path,
size, checksum, indexes, symlink target, upstream/version/license. Include the
external hg19/b37 chrY reference. Export references without credentials and
verify the copy at destination. Supply secure credential/key provisioning
instructions separately; configure intended recipient public keys explicitly.

Add one site configuration, loaded before Python imports and identically
passed into workers. Proposed fields: resource root, native runtime prefix,
shared temp root, environment/container prefixes, transfer tools/remotes,
approved credential file paths and encryption key paths. Resolve derived paths
once; retire misleading duplicate defaults.

Acceptance: no mandatory selected-stage path resolves to another user's home,
Snellius GPFS, old Conda prefix or a resource symlink outside the supplied bundle.

### 3. Complete the scratch contract in ZSlurm and pipeline

Implementation status: **merged to the default pipeline and ZSlurm branches;
unit-tested and exercised by real Spider manager/worker, lease, cleanup and
paired-FASTQ preprocessing canaries; not deployed for a production cohort**.

Implement the discovery/environment/capacity fixes above in a scheduler
checkout, with simulated Snellius/Spider/partial allocations and unavailable
scratch tests. Update all pipeline scratch users, including markdup cleanup,
and add a Spider Snakemake profile with correct optional container binds.
Do not tune CPU/memory reservations or change tool command lines in this step.

Acceptance: two child jobs in one pilot and two partial pilots on one host
cannot claim each other's files or SSD quota; end/failure of one child cannot
delete another's scratch. Required jobs wait when capacity is unavailable;
optional jobs have a correct shared fallback. Leases still release/reacquire
correctly without changing the algorithm's tool parallelism.

### 4. Build software and add a preflight

At final Spider paths, recreate environments including post-deploy actions;
the preprocessing hook builds the native tools/extensions automatically.
Prepare native DeepVariant and cache GLnexus only if needed. Test imports and
actual binary versions from worker envs.

Add a read-only preflight before expensive DAG construction: selected endpoint
dependencies, executable/ABI/ISA checks, references/indexes, intended target
permissions, storage budget configuration and restart manifest consistency.
Check source/output authorization without logging credentials; actual transfer
tests need explicitly designated test objects/destinations, not production
files. No full recursive `stat`/`du` sweep on every start.

Acceptance: preflight reports actionable missing dependencies before hundreds
of jobs are submitted; normal run initialization never builds software in an
active shared checkout.

### 5. Single-pilot Spider smoke tests, then controlled scale-up

Use a separately named manager instance, override the supplied autogrow default
to off or a one-pilot cap, start one deliberately sized pilot and use a small
quota-limited run directory. No Snellius instance name or run state should be
reused accidentally.

- Paired FASTQ/BAM/CRAM, one/multiple readgroups, ERF correction, dedup/rescue,
  WGS/WES, male/female (including skipped Y), chrM on/off.
- Full default native DeepVariant tail, QC/retained stats, Crypt4GH and a
  checksum-verified test upload; GLnexus and HaplotypeCaller only when selected.
- Worker-level rule reparsing with inherited site config and frozen restart
  manifest; completed samples excluded from processing but included in cohort
  statistics. Restart after stage interruption without reopening 1,000 samples.
- Compare alignment flags/tags/checksums (supplementary, secondary, unmapped and
  duplicate records), decoded FASTQs, normalized variants, QC/stat schema and
  record counts against the accepted Snellius baseline. Compressed byte hashes
  alone are not the criterion for BAM/gzip/VCF equivalence.
- Measure scratch high-water usage, shared filesystem traffic, publication
  overlap, active-storage reservation release, average CPU/RSS/PSS and pilot
  entitlement. Revisit budgets only with these measurements, separately from
  the port. Keep final gVCFs/statistics durable until cohort completion.

Acceptance: successful end-to-end canary + restart, clean scratch, no missing
outputs or silent sample skips, bounded reservations and correct RPC/lease
behavior. Only then restore the reviewed compute-autogrow cap (15 in the
supplied Spider ZSlurm configuration).

## Scope and open verification

Spider login/scheduler queries and real `short` allocations were used for the
scratch, RPC, lease, cleanup and paired-FASTQ preprocessing canaries described
above. They do not replace a full alignment/calling, restart, transfer and
container canary. Existing Snellius test results are recorded separately in
`FUSED_CLEANUP_VALIDATION.md`. This plan intentionally does not migrate patient
data, share credentials or alter a running pipeline merely by merging code.
The merge+markdup and CRAM+encryption+upload fusions are implemented on the
pipeline default branch; they still require measured production-like canaries
before deployment at a new site.
