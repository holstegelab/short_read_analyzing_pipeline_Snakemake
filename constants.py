import os
import getpass

#parameters
pj = os.path.join
RESOURCES = '/gpfs/work3/0/qtholstg/hg38_res_v2/'

#region files
INTERVALS_DIR = pj(RESOURCES,'intervals')
MERGED_CAPTURE_KIT_BED = pj(INTERVALS_DIR, 'MERGED_INTERSECT', 'merged_capture_kits_cds.bed')
MERGED_CAPTURE_KIT_IVL = pj(INTERVALS_DIR, 'MERGED_INTERSECT', 'merged_capture_kits_cds.interval_list')
MERGED_CAPTURE_KIT_IVL_CNV = pj(INTERVALS_DIR, 'preprocessed_intervals_for_GATK_CNV', 'merged_capture_kits_cds.interval_list')
INTERSECT_CAPTURE_KIT_BED = pj(INTERVALS_DIR, 'MERGED_INTERSECT', 'intersect_non_focused_capture_kits.bed')

INTERSECT_CAPTURE_KIT_AUTO_BED = pj(INTERVALS_DIR, 'MERGED_INTERSECT_SEX_SPECIFIC', 'intersect_non_focused_capture_kits.auto.bed')
INTERSECT_CAPTURE_KIT_X_BED = pj(INTERVALS_DIR, 'MERGED_INTERSECT_SEX_SPECIFIC', 'intersect_non_focused_capture_kits.chrX.bed')
INTERSECT_CAPTURE_KIT_Y_BED = pj(INTERVALS_DIR, 'MERGED_INTERSECT_SEX_SPECIFIC', 'intersect_non_focused_capture_kits.chrY.bed')

INTERSECT_CAPTURE_KIT_IVL = pj(INTERVALS_DIR, 'MERGED_INTERSECT', 'intersect_non_focused_capture_kits.interval_list')
HARD_MAPPABILITY_TRACK = pj(RESOURCES, 'k24.umap.bed.gz')
TARGETS_BED = pj(INTERVALS_DIR, 'MERGED_INTERSECT', 'gencode_43_cds.bed')
TARGETS_IVL = pj(INTERVALS_DIR,'MERGED_INTERSECT', 'gencode_43_cds.interval_list')
PL_PR_TABLE = pj(RESOURCES, 'ploydi_priors_table_hg38.tsv')
MAIN_CHRS_BED = pj(RESOURCES, 'only_main_chr.bed')
#resource folder with cram reference fasta files
CRAMREFS = pj(RESOURCES,'cram_refs')
GENOME_FILE = pj(INTERVALS_DIR, 'hg38.genome')
PRECOMPUTEED_BED = pj(INTERVALS_DIR, 'precomputed_kits.json')
# dir with this file
SNAKEMAKE_DIR_PATH = os.path.dirname('')

#conda env's paths
CONDA_VERIFYBAMID = 'envs/verifybamid.yaml'
CONDA_MAIN = 'envs/preprocess.yaml'
CONDA_VCF = 'envs/vcf_handling.yaml'
CONDA_PYPY = 'envs/pypy.yaml'
CONDA_KMC = 'envs/kmc.yaml'
CONDA_KRAKEN = 'envs/kraken.yaml'
CONDA_MOSDEPTH = 'envs/mosdepth.yaml'
CONDA_PCA = 'envs/PCA.yaml'
CONDA_GATK_CNV = 'envs/gatk_gcnv.yaml'
CONDA_ANNOVAR = 'envs/annovar.yaml'
CONDA_DRAGMAP = 'envs/dragenos.yaml'
CONDA_ALIGN_FUSED = 'envs/align_fused.yaml'
CONDA_QC_FUSED = 'envs/qc_fused.yaml'
CONDA_CK_FINDER = 'envs/capture_kit_finder.yaml'

DEFAULT_JAVA_OPTIONS = ' -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4 '


#OUTPUT FOLDERS
SOURCEDIR= 'source'
SAMPLEINFODIR= 'sampleinfo'
FETCHDIR= 'fetch'
LOG= 'logs'
BENCH= 'benchmark'
BAM= 'bams'
GVCF= 'gvcf_conv'
GVCF_TAR= 'gvcf_tar'
VCF= 'vcfs'
VCF_Final= 'Final_VCF'
STAT= 'stats'
KRAKEN= 'kraken'
READGROUPS= 'readgroups'
FQ= 'fq'
FQ_BADMAP='fq_badmap'
uBAM= 'uBAM'
uCRAM= 'uCRAM'
KMER= 'kmer'
CRAM= 'cram'
TARGET= 'cnvkit/target'
CNVKIT= 'cnvkit'
DELLY= 'SV_delly'
MULTICOHORT= 'Multicohort'
DEEPVARIANT= 'deepvariant'
chrM= 'chrM_analysis'
GATK_gCNV = 'GATK_gCNV'


#programs
SOFTWARE = pj(RESOURCES, 'software')
DEEPVARIANT_NATIVE_RUNTIME = pj(SOFTWARE, 'deepvariant-1.9.0-native')
gatk= 'gatk'
samtools= 'samtools'
bcftools= 'bcftools'
dragmap= 'dragen-os'
verifybamid2= 'verifybamid2'
java_cnv = pj(SOFTWARE, 'java/jdk-17.0.7/bin/java')
gatk_cnv = pj(SOFTWARE, 'gatk_4.4/build/bundle-files-collected/gatk-package-4.4.0.0-27-gabe8148-SNAPSHOT-local.jar')
annovar = pj(RESOURCES, "annovar/annovar/table_annovar.pl")
annovar_db = pj(RESOURCES, "annovar/annovar/humandb/")
ada = pj(SOFTWARE, 'SpiderScripts/ada/ada')
bcftools_patched = pj(SOFTWARE, 'bcftools-1.8/bcftools')
#custom scripts (encapsulate in srcdir())
BAMMERGE= 'scripts/bam_merge'
BAMCHECK='scripts/bam_check_fastq.py'
BAMSTATS= 'scripts/bam_stats.py'
DECHIMER= 'scripts/bam_dechimer'
DECHIMER_THRESHOLD= 0.005
MERGEPHASE = 'scripts/merge_phasing.py'
MERGEPHASEDIRECT = 'scripts/merge_phasing_direct.py'
CHECKEMPTY = '/gpfs/work3/0/qtholstg/hg38_res_v2/scripts/check_empty.py'
SLOPSCRIPT = 'scripts/slop_start_stop.py'
CAPTURE_KIT_CHECKER = 'scripts/capture_kit_cheker.py'
BED_PRECOMP = 'scripts/precompute_capture_kits.py'
ADA = 'ada'

#path to kmer files
KMER_CHRY= pj(RESOURCES,'kmer_sex/k32.chrY.diff')
KMER_CHRX= pj(RESOURCES,'kmer_sex/k32.chrX.diff')
KMER_CHRM= pj(RESOURCES,'kmer_sex/k32.chrM.diff')
KMER_AUTO= pj(RESOURCES,'kmer_sex/k32.auto.diff')

#path to ref and add ref files
REF =  pj(RESOURCES, 'hg38_phix/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa')
AUTO_ONLY_BED = pj(RESOURCES, 'hg38_phix/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.onlyauto.bed')
REF_DIR = pj(RESOURCES, 'hg38_phix')

SHIFTED_MT= pj(RESOURCES,'MT_ref_shifted')
SHIFTED_MT_fa= pj(RESOURCES,'MT_ref_shifted/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta')
SHIFTED_MT_dict= pj(RESOURCES,'MT_ref_shifted/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict')
SHIFTED_MT_fai= pj(RESOURCES,'MT_ref_shifted/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai')
MT_CHAIN= pj(RESOURCES,'MT_ref_shifted/ShiftBack.chain')
ORIG_MT_fa= pj(RESOURCES,'MT_ref/reference.fasta')
ORIG_MT_dict= pj(RESOURCES,'MT_ref/reference.dict')
ORIG_MT_fai= pj(RESOURCES,'MT_ref/reference.fasta.fai')
ORIG_MT= pj(RESOURCES,'MT_ref')
NUMTs= pj(RESOURCES, 'databases/NUMT_list_hg38.bed')
mask_bed= 'Ref_PhiX_Male/hg38_alt_mask.male.bed'
str_ref= 'hg38_phix/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.str.zip'

REF_FEMALE = REF
REF_FEMALE_DIR = pj(RESOURCES, 'hg38_phix/female/')
REF_FEMALE_STR = pj(RESOURCES, str_ref)
REF_FEMALE_DICT = pj(RESOURCES, 'hg38_phix/female/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.dict')
REF_FEMALE_BED = pj(RESOURCES, 'hg38_phix/female/hg38-ht_mask_bed-v3-female.bed')
REF_FEMALE_HASH = pj(RESOURCES, 'hg38_phix/female/hash_table.cfg')
REF_FEMALE_FAI = pj(RESOURCES, 'hg38_phix/female/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa.fai')

REF_MALE = REF
REF_MALE_DIR = pj(RESOURCES, 'hg38_phix/male/')
REF_MALE_BED = pj(RESOURCES, 'hg38_phix/male/hg38-ht_mask_bed-v3.bed')
REF_MALE_DICT = pj(RESOURCES, 'hg38_phix/male/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.dict')
REF_MALE_HASH = pj(RESOURCES, 'hg38_phix/male/hash_table.cfg')
REF_MALE_STR = pj(RESOURCES, str_ref)
REF_MALE_FAI = pj(RESOURCES, 'hg38_phix/male/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa.fai')

PAR_and_CENTROMERIC = pj(RESOURCES, 'PAR_and_centromeric_regions_hg38.bed')
PAR = pj(RESOURCES, 'hg38_phix/GRCh38_PAR.bed')
#verifybamid files
VERIFYBAMID_EXOME = pj(RESOURCES,'verifybamid/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat')
VERIFYBAMID_WGS = pj(RESOURCES, 'verifybamid/wgs/1000g.phase3.100k.b38.vcf.gz.dat')

# VQSR DBses
DBSNP = pj(RESOURCES,'databases/HG38_dbSNP_v155_updatesd.vcf.gz')
HAPMAP = pj(RESOURCES,'databases/hapmap_3.3.hg38.vcf.gz')
OMNI = pj(RESOURCES,'databases/1000G_omni2.5.hg38.vcf.gz')
KILO_G = pj(RESOURCES,'databases/1000G_phase1.snps.high_confidence.hg38.vcf.gz')
MILLS = pj(RESOURCES,'databases/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz')
DBSNP_INDEL = pj(RESOURCES,'databases/Homo_sapiens_assembly38.known_indels.vcf.gz')

# ANNOTATIONS
REVEL = pj(RESOURCES, 'REVEL/revel_for_bcftools.tab.gz')
REVEL_header = pj(RESOURCES, 'REVEL/revel.hdr')
multiallelic_hdr = pj(SOFTWARE, 'bcftools-1.8/multi_allele.hdr')
CLINVAR = pj(RESOURCES, 'databases/clinvar_20240708_renamed_chrs.vcf.gz')
GNOMAD_4 = pj(RESOURCES, 'databases/gnomad/gnomad4.genomes.full_genome.vcf.gz')
GNOMAD_2 = pj(RESOURCES, 'databases/gnomad_v2/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz')
#path to file with adapters
ADAPTERS = pj(RESOURCES, 'databases/Adapters_illumina.txt')

#windows
WINDOWS = pj(INTERVALS_DIR, 'windows/all.selected.sorted.3.bed')
WINDOWS_ANNOTATED = pj(INTERVALS_DIR, 'windows/all.selected.sorted.bed')

AGH_DCACHE_CONFIG =  pj(RESOURCES, ".agh/agh_processed.conf")
# kraken db
KRAKEN_DB = pj(RESOURCES, 'kraken/pluspf_20230605')

#tmp folders
TMPDIR = 'tmp' #do not use scratch, amount of storage is limited
TMPDIR_ALT = '/scratch-local'
tmpdir = pj(TMPDIR,getpass.getuser())
tmpdir_alternative = pj(TMPDIR_ALT,getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

current_dir = os.getcwd()


# --- expected job runtimes (seconds) -------------------------------------
# Per-rule (WGS, exome) wall-clock estimates, consumed via common.get_time()
# as the Snakemake `time` resource. The zslurm scheduler uses these to pack
# jobs onto nodes and to decide what still fits before a node's walltime ends,
# so they are UPPER estimates: the WGS column is the p99 of observed runtimes,
# rounded up to 5 minutes.
#
# WGS column: measured over 2 production runs on Snellius genoa,
#   ~/data/exome_runs/report-2026-07-14_22-48.tsv (ADSP-FUS1, 1899 samples)
#   ~/data/exome_runs/report-2026-07-06_11-33.tsv (AMP-AD,    1030 samples)
#   successful jobs only (retcode==0); n ranges from 114 to 40966 per rule.
# EXOME column: no exome run was available to measure. For rules whose runtime
#   correlates with input size (r>=0.80 in those reports) it is the fitted
#   line evaluated at 1/6 of the WGS input volume; all other rules keep the
#   WGS value. Treat the exome column as a placeholder to be replaced by
#   measurements from the first real exome run.
#
# Rules absent from this table fall back to the profile default (time=3600).
RUNTIME = {
    'merge_rgs':                               ( 23400,  23400),
    'align_reads':                             ( 16200,   3600),
    'merge_bam_alignment_dechimer':            ( 13800,  13800),
    # Sum of the two phase budgets until fused measurements are available.
    'align_reads_fused':                       ( 30000,  17400),
    'markdup':                                 ( 13200,   2700),
    'adapter_removal':                         ( 12900,   2700),
    'external_alignments_to_fastq':            ( 12300,  12300),
    'external_adapter_fused':                  ( 25200,  15000),
    'sort_bam_alignment':                      ( 10500,   2100),
    'extract_and_tar_deepvariant_level2_wgs':  (  9300,   9300),
    'artifacts_and_oxog_metrics':              (  9300,   9300),
    'bam_qc_fused':                            ( 10800,  10800),
    'split_alignments_by_readgroup':           (  8700,   8700),
    'mCRAM':                                   (  8700,   8700),
    'archive_get':                             (  8400,   8400),
    'dcache_get':                              ( 90000,  90000),
    'dcache_to_active':                        ( 43200,  43200),
    'samtools_stats':                          (  7200,   7200),
    'kmer_reads':                              (  6600,   2100),
    'hs_stats':                                (  6300,   2100),
    'get_validated_sex':                       (  5700,   5700),
    'kmer_sex_fused':                          ( 12300,   7800),
    'deepvariant':                             (  4500,   4500),
    'DVWhatshapPhasingMerge':                  (  3300,   3300),
    'deepvariant_phasing_fused':               (  7800,   7800),
    'bamstats_all_and_exome':                  (  2700,   2700),
    'coverage':                                (  2100,   2100),
    'mutect_bp_resolution_both':               (  1800,   1800),
    'verifybamid':                             (  1500,   1500),
    'extract_and_tar_deepvariant_level2_wes':  (  1200,   1200),
    'kraken':                                  (  1200,   1200),
    'extract_NUMTs_reads':                     (  1200,   1200),
    'chrm_extract_align_fused':                (  2400,   2400),
    'chrm_mutect_tail_fused':                  (  3000,   3000),
    'Encrypt_crams':                           (  1200,   1200),
    'archive_to_active':                       (   900,    900),
    'merge_rgs_badmap':                        (   900,    900),
    'mutect_calls_both':                       (   900,    900),
    'align_chrM_and_NUMTs':                    (   600,    600),
    'copy_to_dcache':                          (   600,    600),
    'tar_badmap_fastqs':                       (   600,    600),
    'extract_chrM_reads':                      (   600,    600),
    'copy_deepvariant_wgs_region_to_dcache':   (   300,    300),
    'copy_deepvariant_wes_region_to_dcache':   (   300,    300),
    'copy_badmap_to_dcache':                   (   300,    300),
    'merge_and_filter_both':                   (   300,    300),
    'chrM_and_numt_read_stats':                (   300,    300),
    'whatsHap_phase_stats':                    (   300,    300),
    'tar_stats_per_sample':                    (   300,    300),
    'bracken':                                 (   300,    300),
    'kraken_summary':                          (   300,    300),
    'finished_sample':                         (   300,    300),
    'start_sample':                            (   300,    300),
    'get_readgroups':                          (   300,    300),
}
