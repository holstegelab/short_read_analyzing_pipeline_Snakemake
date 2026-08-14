# Project MINE CRAM-bestanden ontsleutelen

## Waar staat alles?

Op de target-dCache:

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

Apart aangeleverd, dus niet op dCache:

```text
projectmine_recipient_key.tar.gz.gpg
projectmine_recipient_key.tar.gz.gpg.sha256
```

Het wachtwoord voor het GPG-bestand wordt via een ander veilig kanaal
doorgegeven.

## 1. Controleer en ontsleutel het key-pakket

```bash
sha256sum --check projectmine_recipient_key.tar.gz.gpg.sha256

gpg --output projectmine_recipient_key.tar.gz \
  --decrypt projectmine_recipient_key.tar.gz.gpg
```

GPG vraagt nu om het apart ontvangen wachtwoord.

Pak de key uit:

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

De verwijderde `.tar.gz` bevatte de onversleutelde private key. Bewaar de
uitgepakte `recipient1.key` alleen op beveiligde opslag en zet hem nooit op
dCache.

## 2. Pak het softwarepakket uit

Download `projectmine_cram_decryption_package.zip` van de map `tools/` op
dCache en pak het uit:

```bash
mkdir -m 700 cram_tools
unzip projectmine_cram_decryption_package.zip -d cram_tools
cd cram_tools
sha256sum --check MANIFEST.sha256
cd ..
```

Maak eventueel de meegeleverde omgeving:

```bash
micromamba create -f cram_tools/environment.yml
micromamba activate projectmine-cram-decryption
```

## 3. Download CRAM, index en referentie

Gebruik de apart verstrekte read-only dCache-configuratie. Vervang `REMOTE`
door de remote-naam uit die configuratie:

```bash
mkdir -p encrypted reference

rclone --config dcache.conf copy \
  REMOTE:/processed/mine/mine_hiseq2000_sample_listing/reference \
  reference/

rclone --config dcache.conf copyto \
  REMOTE:/processed/mine/mine_hiseq2000_sample_listing/cram/SAMPLE.mapped_hg38.cram.c4gh \
  encrypted/SAMPLE.mapped_hg38.cram.c4gh

rclone --config dcache.conf copyto \
  REMOTE:/processed/mine/mine_hiseq2000_sample_listing/cram/SAMPLE.mapped_hg38.cram.crai \
  encrypted/SAMPLE.mapped_hg38.cram.crai
```

## 4. Ontsleutel de CRAM

```bash
python cram_tools/decrypt_cram.py \
  --sk projectmine_recipient_key/recipient1.key \
  encrypted/SAMPLE.mapped_hg38.cram.c4gh
```

Dit maakt:

```text
encrypted/SAMPLE.mapped_hg38.cram
```

De meegeleverde `.cram.crai` is de index voor deze ontsleutelde CRAM.

## 5. Controleer of de CRAM leesbaar is

```bash
REF=reference/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa
CRAM=encrypted/SAMPLE.mapped_hg38.cram

samtools quickcheck -v "$CRAM"
samtools view -T "$REF" -c "$CRAM"
```

Het tweede commando leest de hele CRAM en print het aantal alignments.

## 6. Zet de CRAM eventueel terug naar paired FASTQ

```bash
mkdir -p fastq

python cram_tools/bam_revert.py \
  --input "$CRAM" \
  --reference "$REF" \
  --f1 fastq/SAMPLE_R1.fastq.gz \
  --f2 fastq/SAMPLE_R2.fastq.gz \
  --stats fastq/SAMPLE.bam_revert.stats.tsv \
  --threads 4
```

De FASTQ-bestanden kunnen aanzienlijk groter zijn dan de CRAM. Verwijder
private keys, plaintext CRAMs en FASTQ-bestanden wanneer ze niet meer nodig
zijn, volgens het geldende databeleid.
