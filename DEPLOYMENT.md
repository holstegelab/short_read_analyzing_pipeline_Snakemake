# Deployment guide

This guide deploys the same pipeline implementation on another site. It does
not fork biological logic between Snellius and Spider. Work from the exact
component versions in [`deployment/component-lock.yaml`](deployment/component-lock.yaml).

## 1. Obtain the pinned source revisions

Clone the pipeline, Snakemake fork, ZSlurm and native executor into stable,
shared locations visible to the controller and workers. Check out the pinned
commits or the pipeline release tag; do not deploy from a moving default branch.

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

## 4. Create environments and native tools

Choose durable Conda and Apptainer prefixes. The supplied profiles contain the
known project defaults; change them when the target project path differs. Check
the reproducibility inputs first:

```bash
sha256sum -c deployment/environment-inputs.sha256
```

Use Snakemake to create the environments for the selected endpoints so every
matching `*.post-deploy.sh` runs. Important post-deploy behavior includes the
pinned patched DRAGMAP and KMC builds, GATK 4.5 installation and GATK-gCNV
Python setup. `GATK_CNV_ROOT` is derived from the selected site software root.

The Git-ignored native tools must be built in the Python/htslib environment
that runs preprocessing, not copied from another checkout or Python ABI:

```bash
conda run --prefix /path/to/resolved/preprocess-env \
  make -C scripts all-hts PYTHON=python
conda run --prefix /path/to/resolved/preprocess-env \
  python -c 'import sys; sys.path.insert(0, "scripts"); import fastcheck, fastcheck_hts'
```

`all-hts` builds `fastcheck`, `fastcheck_hts`, `bam_merge`, `bam_dechimer` and
`fix_bam_rg_pairs`. This release was tested from a fresh worktree, so those
artifacts were not borrowed from the active pipeline checkout.

Prepare native DeepVariant at its final prefix. Its launchers contain absolute
paths and the directory must not be moved afterwards:

```bash
python scripts/prepare_deepvariant_native.py \
  --prefix /project/holstegelab/Software/short_read_pipeline/deepvariant-1.9.0-native
```

## 5. ZSlurm and Spider boundary

The Spider profile is [`profiles/spider/config.yaml`](profiles/spider/config.yaml).
It intentionally has no `/scratch-node` bind. `partition=compute` is a logical
ZSlurm job class, not a physical Spider partition.

Do not run the production Spider workflow yet. Step 3 of
[`PORTABILITY_PLAN.md`](PORTABILITY_PLAN.md) remains required:

- chief discovery and propagation of the pilot allocation `$TMPDIR`;
- a bounded allocation scratch capacity, including partial pilots;
- unique child-job scratch roots and ownership-safe cleanup;
- removal of remaining `/scratch-node` assumptions in selected rules;
- replacement/disablement of the nonexistent physical `staging` partition;
- real manager/worker RPC, lease and failure/restart smoke tests.

Snellius archive sources additionally require its `daget`, `dals` and
`darelease` commands. On Spider, use a tested dCache/S3 route unless an
equivalent archive backend is installed.

## 6. Start only after the readiness checks pass

For Snellius, an explicit site-configured invocation is:

```bash
export SHORT_READ_SITE_CONFIG=/path/to/checkout/config/sites/snellius.yaml
python -m snakemake \
  --profile /path/to/checkout/profiles/zslurm \
  --snakefile /path/to/checkout/Snakefile \
  --config END_POINT=Genotype caller=Deepvariant Combine_gVCF_method=GLnexus
```

For Spider, use `profiles/spider` only after the scratch work above passes a
single-pilot canary. Before scaling, verify:

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
