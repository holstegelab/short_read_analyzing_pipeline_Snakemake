# Deployment guide

This guide deploys the same pipeline implementation on another site. It does
not fork biological logic between Snellius and Spider. Work from the exact
component versions in [`deployment/component-lock.yaml`](deployment/component-lock.yaml).

## 1. Obtain the pinned source revisions

Clone the pipeline, Snakemake fork, ZSlurm and native executor into stable,
shared locations visible to the controller and workers. Check out the pinned
commits or the pipeline release tag; do not deploy from a moving default branch.

Do not silently substitute stock Snakemake for the pinned fork. Besides the
historical duplicate-rule bypass, that fork contains independent fixes used by
this workflow, including checkpoint-safe temp cleanup, Conda support for
`run:` rules, incomplete-metadata cleanup, and the symlinked-working-directory
Apptainer bind. This list is not a replacement for the fork history. The
workflow no longer needs duplicate rule names, so the pinned release revision
restores only that validation on top of the pinned fork commit; it does not
rebase the deployment onto stock Snakemake.

The pipeline remains compatible with the immediately preceding fork revision
`83c79c1d3d8dfd68bc42f4ac9ffa3dbb97d8c8ed`: the Snellius release canary
built its DAG and completed real worker jobs with that revision. Updating to
the revision in the component lock is nevertheless recommended because it
restores duplicate-rule validation and prevents the ambiguity from returning.

Install the Snakemake fork and executor plugin into the same controller
environment. Install ZSlurm into the environment used by the manager and pilot
chiefs. Confirm imports, rather than assuming the checkout is the imported code:

```bash
python -c 'import snakemake; print(snakemake.__file__)'
python -c 'import snakemake_executor_plugin_zslurm as p; print(p.__file__)'
python -c 'import zslurm_shared; print(zslurm_shared.__file__)'
python -m snakemake --version
```

The executor revision is required. It carries workflow `envvars:` and
storage-provider values in the XML-RPC environment used by shell-free ZSlurm
workers. Without it, a worker can lose the selected site or restart manifest.

## 2. Install resources without credentials

Download the required resource archive plus the manifests archive. Use `core`
for the production alignment/calling/QC path; add `joint_annotation` for the
default GLnexus annotation tail and `gcnv_optional` only when those rules are
selected. Verify the adjacent SHA-256 before extraction:

```bash
sha256sum -c short_read_pipeline_resources_core_c386a7d_20260914.tar.zst.sha256
tar --zstd -xf short_read_pipeline_resources_core_c386a7d_20260914.tar.zst \
  -C /project/holstegelab/Software/short_read_pipeline
```

The archive member paths start at `hg38_res_v2/`. Verify the installed files
against the `*.paths` manifest from the manifests archive. The release hashes
are also recorded in the component lock.

Install the `hg19_b37_chry` supplement too. It places the corrected FASTA used
for exceptional CRAM dictionaries at
`hg38_res_v2/cram_refs/hg19_b37chrY.fa` together with its FAI. The original
core archive was audited before that external file was folded into the resource
tree; the supplement and both member hashes are pinned in the component lock.

Never extract `.agh` tokens or private `.c4gh` keys into a public/shared
resource installation by accident. The original Crypt4GH trust unit is in a
separate password-encrypted archive on a separate, private endpoint. Its
password is delivered out of band. A site may instead provision new approved
sender/recipient keys.

## 3. Create the site configuration

Copy [`config/sites/spider.example.yaml`](config/sites/spider.example.yaml) to
a protected deployment location and edit every path, remote and credential
filename. Snellius has a concrete reference in
[`config/sites/snellius.yaml`](config/sites/snellius.yaml). The YAML contains
paths only; it must never contain a token, private-key contents or a password.

Schema 1 separates:

- `paths`: resources, software, shared temporary storage, environment/cache
  prefixes and CRAM reference fallbacks;
- `storage`: dCache read/processed config files and the default read remote;
- `tools`: transfer and archive command locations;
- `encryption`: sender, recipients and decryption key-file locations.

Select it before invoking Snakemake:

```bash
export SHORT_READ_SITE_CONFIG=/absolute/protected/path/spider.yaml
python -c 'from site_config import configure; print(dict(configure().values))'
```

The loader rejects unknown keys, unsupported schema versions and unresolved
environment variables. It resolves the selector to an absolute path and the
root Snakefile declares it as an `envvars:` dependency. Every worker therefore
loads the same file before importing modules that derive resource paths.

CLI workflow settings still take precedence for their established keys. For
example, `--config deepvariant_native_prefix=...` overrides the corresponding
site default for that invocation.

## 4. Bootstrap the site once

This section is for the site administrator, not for every cohort run. Choose
durable Conda and Apptainer prefixes in the site configuration and check the
reproducibility inputs first:

```bash
sha256sum -c deployment/environment-inputs.sha256
```

Install a named, site-bound Snakemake profile once. The installer copies the
Spider execution settings, records this pinned Snakefile, takes the Conda and
Apptainer prefixes from the site configuration, passes that configuration to
every controller/worker parse, and selects the production-safe `mtime` rerun
trigger:

```bash
PIPELINE=/absolute/path/to/short_read_analyzing_pipeline_Snakemake
python "$PIPELINE/scripts/install_site_profile.py" \
  --site-config /absolute/protected/path/spider.yaml \
  --profile-name zslurm2
```

After that one-time step, both explicit environment preparation and normal
workflow execution use the named profile. Run environment preparation from a
configured cohort directory containing its sample TSV and sidecars; the real
DAG determines exactly which environments are required. An empty directory has
no sample jobs and therefore no environments to create. Explicit preparation
is optional because a normal workflow invocation creates missing environments:

```bash
cd /absolute/path/to/configured/cohort-run
snakemake --profile zslurm2 --cores 1 --conda-create-envs-only \
  --config END_POINT=gVCF caller=Deepvariant
```

The `preprocess` environment's post-deploy hook automatically builds
`fastcheck`, `fastcheck_hts`, `bam_merge`, `bam_dechimer` and
`fix_bam_rg_pairs` in an isolated temporary directory and atomically installs
them into this pipeline checkout. Native-source checksums participate in the
environment hash. Do not locate a hashed environment manually and do not run
`conda run --prefix ... make`.

On Spider, a non-interactive Bash descended from an SSH login can reload the
Conda function from `.bashrc`. If that function still points at an old base
Conda, it takes precedence over a newer executable at the front of `PATH`.
Snakemake requires Conda 24.7.1 or newer. This is a controller installation
problem, not part of each pipeline invocation. Activate the intended
controller and verify what it actually sees:

```bash
source /absolute/path/to/controller-environment/bin/activate
snakemake --version
conda --version
type -a conda
```

The controller may be a Python environment containing Snakemake while Conda is
provided separately; it does not have to contain `bin/conda`. If `conda
--version` is too old, update/fix the site's login Conda before continuing.
Merely checking `which conda` does not expose a shell function.

There is intentionally no second deployment Snakefile or duplicated list of
environments. To prepare another endpoint, invoke the main Snakefile with that
endpoint against a representative configured cohort. Every `*.post-deploy.sh`
belonging to an environment selected by that DAG runs when the environment is
first created.
Important post-deploy behavior includes the pinned patched DRAGMAP and KMC
builds, GATK 4.5 installation and GATK-gCNV Python setup. `GATK_CNV_ROOT` is
derived from the selected site software root.

Prepare native DeepVariant at its final prefix. Its launchers contain absolute
paths and the directory must not be moved afterwards:

```bash
python scripts/prepare_deepvariant_native.py \
  --prefix /project/holstegelab/Software/short_read_pipeline/deepvariant-1.9.0-native
```

## 5. ZSlurm and Spider boundary

The Spider profile is [`profiles/spider/config.yaml`](profiles/spider/config.yaml).
It intentionally has no `/scratch-node` bind. `partition=compute` is a logical
ZSlurm job class, not a physical Spider partition. The matching ZSlurm Spider
profile discovers the pilot allocation's `$TMPDIR`, caps the advertised SSD
entitlement, and exports one private `ZSLURM_SCRATCH_DIR` per child job.

The pinned ZSlurm release discovers Spider's allocation `$TMPDIR`,
gives every ZSlurm child a private directory through
`ZSLURM_SCRATCH_DIR`, rewrites all ordinary temp variables to that directory,
and cleans only that child directory. Partial pilots advertise 73 GiB per
allocated core, bounded again by actual filesystem free space; a 30-core pilot
therefore advertises at most 2,190 GiB. Generic `TMPDIR` discovery is enabled
only in the Spider ZSlurm site file; it is not a cross-site pipeline guess. The
nonexistent physical `staging` partition and its autogrow path are disabled.
Demand-backed compute autogrow is enabled with a cap of 15 running plus queued
pilots. GPFS RDMA telemetry is disabled for CephFS.

A real Spider `short` allocation confirmed that the pipeline resolves the
exported child path below the allocation's private `/tmp` XFS bind. A second
canary exercised manager/worker RPC, dynamic lease resize, private scratch
cleanup and a real two-job paired-FASTQ preprocessing DAG through the native
executor. Its output FASTQs passed gzip validation and the expected read/base
counts. The adapter rule used its normal five-core claim; only its walltime was
reduced to fit Spider's 30-minute `short` limit.

This is not yet a full alignment/calling production certification. The canary
account could not read the pre-existing resource/software tree under
`/project/cardseq`, which belongs to another Spider project; CardSeq is not a
pipeline component. Alignment, DeepVariant, dCache writes, failure/restart and
Apptainer therefore still require site-owned canaries before allowing the
configured autogrow limit to scale a cohort run. Disable compute autogrow or
temporarily cap it at one pilot during those canaries.

Snellius archive sources additionally require its `daget`, `dals` and
`darelease` commands. On Spider, use a tested dCache/S3 route unless an
equivalent archive backend is installed.

### Large-write quota backoff (Snellius only)

The shared gate in `scripts/fastq_quota.py`, exposed through
`scripts/storage_quota.py` for direct shell rules, checks only the destination
Snellius fileset/mount with `mmlsquota -e -Y -j FILESET FILESYSTEM`. The default
is `storage_quota_mode=auto`; home and scratch quotas are not queried. GPFS
`in_doubt` usage is included. Checks precede bulk writes in:

- `external_adapter_fused`: publication of raw and trimmed FASTQs from SSD;
- `markdup`: publication of the completed BAM, index and duplicate statistics;
- `deepvariant_phasing_fused`: shared-scratch computation and final publication;
- `start_sample_archive/dcache/s3`: new source materializations;
- native-FASTQ `adapter_removal`: before preparing FASTQs on shared storage;
- `split_alignments_by_readgroup`: actual splitting or ERF correction;
- BAM restoration: encrypted download, decryption, and CRAM-to-BAM expansion.

Validated/adopted materializations, active inputs, single-RG hardlinks, shared
tape staging batches, and alignment/upload/cleanup paths are exempt.
DeepVariant on local SSD checks project quota only when publishing its results.

At 97% or above (storage or inodes), or if a quota query fails, the job keeps
any SSD data/reservation, requests a small existing CPU/RAM lease, and retries
after a random 4–6 minutes. Once quota recovers, it waits a random 0–10 minutes
and checks again before writing. Each caller's complete CPU/RAM work profile
is reacquired through the existing lease interface, capped at its original
reservation, and quota is rechecked after reacquisition. Healthy/disabled
checks do not query the scheduler; initially healthy jobs proceed immediately.
Local execution without a lease can wait but retains its resources. No
manager/chief changes, shared locks or atomic space reservations are involved;
this is best-effort
backpressure, not a guarantee against other writers filling the quota.

The Spider profile explicitly sets `storage_quota_mode=disabled` and
`external_adapter_quota_mode=disabled`; auto mode also skips non-Snellius
project paths. Other profiles can disable all gates with
`--config storage_quota_mode=disabled`, or require recognized Snellius project
destinations with `...=required`. `external_adapter_quota_mode` overrides the
global setting for that rule alone. Already-running processes do not reload
these changes; already-expanded queued shell commands also retain their old
form until the workflow reparses. Frozen checkouts must be updated separately.
Waiting remains subject to the node allocation's lifetime.

## 6. Start only after the readiness checks pass

Once the site profile is installed, a normal cohort start is deliberately
short and does not repeat any deployment paths:

```bash
cd /absolute/path/to/configured/cohort-run
snakemake --profile zslurm2 --zslurm-instance zslurm_site_controller
```

Endpoint/caller overrides can still be appended with `--config`. The installed
profile already supplies the Snakefile, site configuration, Conda/Apptainer
prefixes, executor and normal safety settings. Snellius can use the installer
with `--profile-template "$PIPELINE/profiles/zslurm/config.yaml"`; its
established `zslurm2` profile remains valid.

For Spider, use `profiles/spider` at production scale only after the scratch
work above passes a single-pilot canary. Temporarily disable compute autogrow
or cap it at one pilot for that canary; then restore the configured cap of 15.
Before scaling, verify:

- every selected reference/index and executable resolves below the intended
  project or protected credential root;
- controller and worker import the pinned Snakemake/plugin/ZSlurm revisions;
- native tools load under the worker ABI and expected ISA;
- source and destination remotes can list/read/write the designated test area;
- no secret is printed in logs or readable by group/other;
- FASTQ, BAM and CRAM canaries cover one/multiple readgroups, supplementary,
  secondary and unmapped records, restart, output checksums and scratch cleanup.

See [`RESTARTING.md`](RESTARTING.md) before resuming an existing run. Never
reuse another site's manager instance name, run directory, frozen restart
manifest or `.snakemake` metadata.
