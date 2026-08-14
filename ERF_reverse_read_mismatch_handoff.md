# ERF Reverse-Read BAM/FASTQ Mismatch Handoff

## Scope

This note summarizes the current investigation into BAM versus FASTQ mismatches seen for ERF/external samples during the merge/dechimer stage.

The failing path is in `merge_bam_alignment_dechimer` in [Aligner.smk](Aligner.smk), where:

1. adapter-removed FASTQs are aligned with `dragen-os`
2. the aligned BAM is streamed into `scripts/bam_merge`
3. `bam_merge` compares BAM records against the adapter-removed FASTQs and restores clipped sequence/quality tags
4. `bam_stats_compare_hts.py` later compares aggregate BAM/CRAM content against FASTQ checksums

The aggregate checksum/count mismatch is downstream fallout from `bam_merge` aborting early on the first fragment-level BAM/FASTQ mismatch.

## What Was Already Fixed Earlier

The earlier ERF-specific QNAME issue in `scripts/fix_bam_rg_pairs.c` was fixed so that legacy names like `/1` and `/3` are canonicalized to the same base QNAME before `samtools fastq` pairing.

That fixed the all-singletons issue in `external_alignments_to_fastq`.

This current problem appears different.

## Current Observations

The mismatch is seen mostly in ERF/external samples and appears infrequent per sample, but common enough across these samples to matter.

The current fragment-level error messages look like this.

### Example 1

```text
Quality mismatch BAM-FASTQ in fragment HWI-IT879:71:6:2308:1126:118894#0:10 (rev)
  BAM qual (96): #DDHDFHFGGIDIE>GHHEHI;81CGFG@@D8)8;DD=?<FHFCHH4@F@GGGIA=A>6?3??BCC@@;@;;>;;CA@>>A>>CC35@C#######
  FASTQ qual raw (101): #######C@53CC>>A>>@AC;;>;;@;@@CCB??3?6>A=AIGGG@F@4HHCFHF<?=DD;8)8D@@GFGC18;IHEHHG>EIDIGGFHFDHDDDA=1#?
  FASTQ qual oriented (101): ?#1=ADDDHDFHFGGIDIE>GHHEHI;81CGFG@@D8)8;DD=?<FHFCHH4@F@GGGIA=A>6?3??BCC@@;@;;>;;CA@>>A>>CC35@C#######
  FASTQ qual compared (96): DDDHDFHFGGIDIE>GHHEHI;81CGFG@@D8)8;DD=?<FHFCHH4@F@GGGIA=A>6?3??BCC@@;@;;>;;CA@>>A>>CC35@C#######
BAM/CRAM (htslib) vs FASTQ comparison:
  read1_nrow: bam=122776 expect=2379564 match=False
  read2_nrow: bam=122776 expect=2379564 match=False
  read1_nbases: bam=12400376 expect=240335964
  read2_nbases: bam=12400376 expect=240335964
  read1_seq_checksum: bam=0x16b22e9fc6689012 expect=0x4501eda155408549 match=False
  read1_qual_checksum: bam=0x36c19c3f3330c518 expect=0xb8cf4dd855a2e8be match=False
  read2_seq_checksum: bam=0x3748791e6eab6df1 expect=0xade00a91cc61c5ae match=False
  read2_qual_checksum: bam=0x14f79b1febf6cd72 expect=0x7662357d0b029524 match=False
```

### Example 2

```text
Quality mismatch BAM-FASTQ in fragment HWI-ST164:320:1:1101:1074:70095#0:10 (rev)
  BAM qual (99): CF#2=BFHHHHJJFIJFIGIJJIHIJIIJIJBHIJJJIFGBGIIEHIJICHHIIJDGIECH?BFFFFFDCECCCCDDDDDCDDDDDD<ABDDCD#DDC8
  FASTQ qual raw (101): 8CDDDDCDDBA<DDDDDDCDDDDDCCCCECDFFFFFB?HCEIGDJIIHHCIJIHEIIGBGFIJJJIHBJIJIIJIHIJJIGIFJIFJJHHHHFB=2#FCCB
  FASTQ qual oriented (101): BCCF#2=BFHHHHJJFIJFIGIJJIHIJIIJIJBHIJJJIFGBGIIEHIJICHHIIJDGIECH?BFFFFFDCECCCCDDDDDCDDDDDD<ABDDCDDDDC8
  FASTQ qual compared (99):     CF#2=BFHHHHJJFIJFIGIJJIHIJIIJIJBHIJJJIFGBGIIEHIJICHHIIJDGIECH?BFFFFFDCECCCCDDDDDCDDDDDD<ABDDCDDDDC8
BAM/CRAM (htslib) vs FASTQ comparison:
  read1_nrow: bam=52 expect=1953784 match=False
  read2_nrow: bam=52 expect=1953784 match=False
  read1_nbases: bam=5252 expect=197332184
  read2_nbases: bam=5252 expect=197332184
  read1_seq_checksum: bam=0x94a7b86abba9c9c6 expect=0x7d54341a9fdd5c3 match=False
  read1_qual_checksum: bam=0xa3acb26ccd6a52fe expect=0x8efb79f606b8626f match=False
  read2_seq_checksum: bam=0xe4adbc51442a996f expect=0xd518952f04de6ebd match=False
  read2_qual_checksum: bam=0x9fcac1a01d8ed32d expect=0
```

### Example 3

```text
Quality mismatch BAM-FASTQ in fragment HWI-ST164:321:5:1208:1088:35245#0:11 (rev)
  BAM qual (99): @#4ADD?FHAH>?EHGHIIIIIIIDHDGII?FGBGGFBGGGIIIH'5<>ABCCC>ACCCC?CCCCCCACACCC:@CC:C@CCC<@<B3ABBCCCC#B##
  FASTQ qual raw (101): ##B@CCCCBBA3B<@<CCC@C:CC@:CCCACACCCCCC?CCCCA>CCCBA><5'HIIIGGGBFGGBGF?IIGDHDIIIIIIIHGHE?>HAHF?DDA4#@@@
  FASTQ qual oriented (101): @@@#4ADD?FHAH>?EHGHIIIIIIIDHDGII?FGBGGFBGGGIIIH'5<>ABCCC>ACCCC?CCCCCCACACCC:@CC:C@CCC<@<B3ABBCCCC@B##
  FASTQ qual compared (99): @#4ADD?FHAH>?EHGHIIIIIIIDHDGII?FGBGGFBGGGIIIH'5<>ABCCC>ACCCC?CCCCCCACACCC:@CC:C@CCC<@<B3ABBCCCC@B##
BAM/CRAM (htslib) vs FASTQ comparison:
  read1_nrow: bam=157875 expect=2396661 match=False
  read2_nrow: bam=157875 expect=2396661 match=False
  read1_nbases: bam=15945375 expect=242062761
  read2_nbases: bam=15945375 expect=242062761
  read1_seq_checksum: bam=0x917e9cbae5b763f4 expect=0xc6aa221e972b4077 match=False
  read1_qual_checksum: bam=0x411b3c7e6b6b4202 expect=0xfdb01005eb9ba759 match=False
  read2_seq_checksum: bam=0x1e5f00698e2fd7ad expect=0x4dbc11286f6e4353 match=False
  read2_qual_checksum: bam=0xbc0522046200a4cb expect=0x69471a5cd04fe1da match=False
```

## Interpreting These Examples

The important signal is not the later checksum block.

That block is expected to look catastrophically wrong once `bam_merge` aborts on the first bad fragment, because the aggregate BAM/CRAM comparison then sees only a truncated prefix of the stream.

The useful signal is the fragment-level mismatch.

In the examples above:

1. the read names are matching
2. the mismatch is on reverse reads
3. the BAM quality and the FASTQ-derived quality are mostly the same
4. the disagreement is concentrated near one end of the reverse-oriented window
5. the raw/oriented FASTQ lengths are often 101 while the BAM quality length is shorter, for example 96 or 99

That pattern does **not** look like random wrong-pair matching.

It looks more like a wrong assumption about which substring of the original FASTQ corresponds to the BAM read for some reverse/clipped ERF-derived reads.

## Current Working Hypothesis

The most likely bug is in `scripts/bam_merge.c` and/or `scripts/bam_merge.py` during restoration of clipped sequence/quality for reverse reads.

The current code assumes that for reverse reads, after orienting the FASTQ into BAM orientation, the BAM read should match the **suffix** of that oriented FASTQ, and the missing clipped sequence belongs on the opposite side.

If that assumption is wrong for a subset of ERF reads, then:

1. the fragment lookup still succeeds
2. the read largely matches
3. the mismatch appears only near the clipping boundary
4. aggregate checksums become wrong only because the merge aborts early

## Could This Be Caused By Reverse-Complement Handling Before Adapter Removal?

Possibly, but it does not look like the first explanation to try.

Why it is possible:

1. these are reverse-read mismatches
2. the BAM/FASTQ comparison is sensitive to sequence orientation and clipping restoration
3. an upstream transformation that changed orientation assumptions could produce consistent reverse-window errors later

Why it seems less likely as the primary root cause:

1. the failures are infrequent within a sample rather than global
2. the read names still match and the qualities mostly line up
3. the mismatches cluster near the clipped boundary rather than the entire read being reversed or unrelated
4. if reads had been globally reverse-complemented incorrectly before adapter removal, a much broader failure pattern would be expected

So:

- `adapter_removal` or pre-adapter FASTQ orientation should still be checked
- but the strongest current hypothesis remains a **reverse-read restoration-window bug** in `bam_merge`

## Other Plausible ERF-Specific Explanations

The fact that this seems to happen mainly in ERF/external samples suggests at least one of these may matter:

1. the external BAM/CRAM content already has hard clipping / soft clipping / OQ-style history that does not match the assumptions in `bam_merge`
2. recalibrated or externally processed reads may have quality strings that correspond to a different clipping state than the current sequence window
3. `fix_bam_rg_pairs.c` may have repaired flags and QNAMEs correctly, but the resulting BAM still carries read records whose clipping/orientation pattern differs from what `bam_merge` expects
4. some ERF reads may have undergone a transform before the pipeline sees them, for example recalibration or reversion details, that changes qualities only for the aligned/clipped portion

## Diagnostics Added So Far

The merge executable `scripts/bam_merge` was rebuilt multiple times with more detailed mismatch diagnostics.

The current diagnostics in `scripts/bam_merge.c` print:

1. fragment name
2. mismatch mode, for example `rev`
3. BAM flag
4. CIGAR
5. BAM sequence length
6. FASTQ sequence length
7. computed clip length
8. BAM sequence
9. FASTQ raw sequence
10. FASTQ oriented sequence
11. FASTQ compared sequence window
12. BAM quality
13. FASTQ raw quality
14. FASTQ oriented quality
15. FASTQ compared quality window
16. the actual offset where the BAM sequence is found in the oriented FASTQ sequence
17. the offset that the code assumed

This last offset comparison is the key next discriminator.

If the actual offset differs from the assumed offset, then the current reverse-window logic is wrong.

## What Information Is Needed Next

To take the next step efficiently, the following would be most useful.

### Highest Priority

1. A failing diagnostic block from the **current** rebuilt `scripts/bam_merge`, including:
   - `flag`
   - `cigar`
   - `bam_len`
   - `fastq_len`
   - `computed_clip`
   - `BAM seq`
   - `FASTQ seq raw`
   - `FASTQ seq oriented`
   - `FASTQ seq compared`
   - `FASTQ oriented match offset`
   - `FASTQ assumed match offset`
   - all four quality lines

2. For one failing fragment, the original source alignment record before any ERF repair if available:
   - the original BAM/CRAM SAM line for that read
   - preferably both mates

3. For the same fragment, the corresponding FASTQ entries from:
   - the extracted FASTQ from `external_alignments_to_fastq`
   - the adapter-removed FASTQ used for alignment

### Very Helpful

4. Whether the failing sample was provided as:
   - `bam`
   - `cram`
   - `recalibrated_bam`
   - `recalibrated_cram`
   - `extracted_bam`
   - `extracted_cram`

5. Whether the sample is known to be recalibrated or otherwise externally transformed before entering the pipeline.

6. The sample config row or metadata for one failing ERF sample, including any special flags.

7. If available, confirmation whether the same fragment mismatch exists in:
   - the repaired readgroup BAM after `fix_bam_rg_pairs`
   - the aligned BAM from `align_reads`

### Nice To Have

8. The AdapterRemoval settings/log for a failing sample/readgroup.

9. A count of how often this happens per sample or per million reads.

10. Whether the same phenomenon ever appears in non-ERF samples.

## Most Likely Next Implementation Step

If the newly added offset diagnostics show that the BAM sequence occurs at a different position in the oriented FASTQ than the code assumes, then the right fix is probably:

1. stop assuming that reverse reads always map to the final `bam_len` bases of the oriented FASTQ
2. locate the BAM sequence within the oriented FASTQ
3. derive restored prefix/suffix tags from the **actual** match interval instead of an assumed interval
4. keep the current hard failure if multiple equally plausible placements exist

That would be a root-cause fix rather than a relaxation.

## Status

Current status: unresolved, but narrowed.

What seems unlikely now:

1. pure checksum bug
2. pure QNAME pairing bug
3. random FASTQ/BAM desynchronization

What seems most likely now:

1. reverse-read clipping restoration bug in `bam_merge`
2. potentially triggered only by ERF/external sample characteristics such as preprocessed or recalibrated alignments
