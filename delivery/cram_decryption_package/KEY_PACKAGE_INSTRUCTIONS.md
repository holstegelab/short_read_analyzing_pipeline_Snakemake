# Project MINE encrypted recipient key

The separately delivered file `projectmine_recipient_key.tar.gz.gpg` contains
the private Crypt4GH recipient key needed by `decrypt_cram.py`. It is protected
with a random password using GnuPG symmetric AES-256 encryption. Obtain that
password through a different secure channel; it must not be emailed with the
encrypted key package.

Do not upload the encrypted or extracted key package to dCache.

## 1. Decrypt the key package

On a workstation with GnuPG installed:

```bash
gpg --output projectmine_recipient_key.tar.gz \
  --decrypt projectmine_recipient_key.tar.gz.gpg
```

GnuPG will ask for the package password interactively.

For a headless system, put the separately received password in an owner-only
file and use:

```bash
chmod 600 /secure/path/key-package.password
gpg --batch --yes --pinentry-mode loopback \
  --passphrase-file /secure/path/key-package.password \
  --output projectmine_recipient_key.tar.gz \
  --decrypt projectmine_recipient_key.tar.gz.gpg
```

## 2. Extract and verify the recipient key

```bash
mkdir -m 700 projectmine_recipient_key
tar -xzf projectmine_recipient_key.tar.gz \
  -C projectmine_recipient_key
cd projectmine_recipient_key
sha256sum --check MANIFEST.sha256
chmod 600 recipient1.key
cd ..
rm projectmine_recipient_key.tar.gz
```

The unencrypted `.tar.gz` also contains the private key, so remove it
immediately after extraction. `recipient1.pub` is included only to identify
the matching public recipient.

## 3. Use the key

From the software package directory:

```bash
python decrypt_cram.py \
  --sk /secure/path/projectmine_recipient_key/recipient1.key \
  encrypted/SAMPLE.mapped_hg38.cram.c4gh
```

Keep the private key owner-readable only. Delete the extracted key and any
local password file when access is no longer required, according to the
applicable data-management policy.
