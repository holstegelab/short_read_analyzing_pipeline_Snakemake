# Project MINE CRAM delivery

This package contains the tools and instructions for decrypting the
Crypt4GH-encrypted, mapped hg38 CRAM files produced by the short-read pipeline.
It can also restore paired FASTQ files from a decrypted CRAM.

## Security

This package deliberately contains **no** dCache macaroon, private Crypt4GH
key, or key passphrase. Obtain the dCache read macaroon and the private key
matching one of the encryption recipients through separate approved channels.
Do not email a private key or passphrase with this package or the data.

Store credentials with owner-only permissions:

```bash
chmod 600 recipient.key dcache.conf
```

## Contents

- `decrypt_cram.py`: atomic Crypt4GH decryption wrapper
- `bam_revert.py`: standalone CRAM/BAM/SAM to paired FASTQ converter
- `environment.yml`: reproducible Conda/Mamba environment
- `KEY_PACKAGE_INSTRUCTIONS.md`: instructions for the separately encrypted key
- `MANIFEST.sha256`: SHA-256 checksums for this package

The target dCache layout is:

```text
/processed/mine/mine_hiseq2000_sample_listing/
├── cram/
│   ├── SAMPLE.mapped_hg38.cram.c4gh
│   └── SAMPLE.mapped_hg38.cram.crai
├── reference/
│   ├── GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa
│   ├── GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa.fai
│   └── GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.dict
└── tools/
    ├── bam_revert.py
    └── projectmine_cram_decryption_package.zip
```

## 1. Install the software

With Micromamba:

```bash
micromamba create -f environment.yml
micromamba activate projectmine-cram-decryption
```

Conda or Mamba can be used instead by replacing `micromamba` in these
commands.

## 2. Obtain the recipient key

The private Crypt4GH key is delivered separately as
`projectmine_recipient_key.tar.gz.gpg`. Follow
`KEY_PACKAGE_INSTRUCTIONS.md` and obtain its random password through a
different secure channel. Never upload the key package or extracted private
key to dCache.

## 3. Download the data and reference

Use the separately supplied read-only dCache macaroon configuration. The
remote name before the colon depends on that configuration:

```bash
mkdir -p encrypted reference
rclone --config dcache.conf copy \
  REMOTE:/processed/mine/mine_hiseq2000_sample_listing/reference reference/
rclone --config dcache.conf copy \
  REMOTE:/processed/mine/mine_hiseq2000_sample_listing/cram encrypted/
```

For a single sample, use `rclone copyto` with its complete remote and local
file names. Keep the `.cram.crai` file beside the decrypted `.cram`; this index
was built for the plaintext CRAM and is not encrypted.

## 4. Decrypt a CRAM

The default output name is the input name without `.c4gh`:

```bash
python decrypt_cram.py \
  --sk /secure/path/recipient.key \
  encrypted/SAMPLE.mapped_hg38.cram.c4gh
```

Crypt4GH asks interactively for the key passphrase when the private key is
encrypted. For an unattended batch job, a protected one-line passphrase file
can be used:

```bash
chmod 600 /secure/path/key.passphrase
python decrypt_cram.py \
  --sk /secure/path/recipient.key \
  --password-file /secure/path/key.passphrase \
  encrypted/SAMPLE.mapped_hg38.cram.c4gh
```

The script writes to a private temporary file and only renames it to the final
CRAM after successful decryption. It refuses to replace an existing output
unless `--replace` is supplied.

## 5. Check that the CRAM is readable

Set the reference path once:

```bash
REF=reference/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa
CRAM=encrypted/SAMPLE.mapped_hg38.cram
samtools quickcheck -v "$CRAM"
samtools view -T "$REF" -c "$CRAM"
```

The second command reads all alignment records and prints their count. A
non-zero exit status or diagnostic on standard error means the check failed.

## 6. Restore paired reads when needed

The delivered CRAM is mapped and therefore needs the delivered reference for
decoding. `bam_revert.py` name-collates the records, restores hard-clipped
bases and original read orientation, and writes paired gzipped FASTQ:

```bash
mkdir -p fastq
python bam_revert.py \
  --input "$CRAM" \
  --reference "$REF" \
  --f1 fastq/SAMPLE_R1.fastq.gz \
  --f2 fastq/SAMPLE_R2.fastq.gz \
  --stats fastq/SAMPLE.bam_revert.stats.tsv \
  --threads 4
```

Run this on storage with enough free space for both FASTQ files. Delete
plaintext CRAM and FASTQ files securely when they are no longer required,
according to the applicable data-management policy.
