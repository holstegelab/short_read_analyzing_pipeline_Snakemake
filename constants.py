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
CONDA_GATK_CNV = pj(RESOURCES, 'software', 'gatk_4.4', 'build', 'gatkcondaenv.yml')
CONDA_ANNOVAR = 'envs/annovar.yaml'
CONDA_DRAGMAP = 'envs/dragenos.yaml'
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
VCF= 'vcfs'
VCF_Final= 'Final_VCF'
STAT= 'stats'
KRAKEN= 'kraken'
READGROUPS= 'readgroups'
FQ= 'fq'
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
BAMMERGE= 'scripts/bam_merge.py'
BAMCHECK='scripts/bam_check_fastq.py'
BAMSTATS= 'scripts/bam_stats.py'
DECHIMER= 'scripts/bam_dechimer.py'
DECHIMER_THRESHOLD= 0.005
MERGEPHASE = 'scripts/merge_phasing.py'
MERGEPHASEDIRECT = 'scripts/merge_phasing_direct.py'
CHECKEMPTY = '/gpfs/work3/0/qtholstg/hg38_res_v2/scripts/check_empty.py'
SLOPSCRIPT = 'scripts/slop_start_stop.py'
CAPTURE_KIT_CHECKER = 'scripts/capture_kit_cheker.py'
BED_PRECOMP = 'scripts/precompute_capture_kits.py'

#path to kmer files
KMER_CHRY= pj(RESOURCES,'kmer_sex/k32.chrY.diff')
KMER_CHRX= pj(RESOURCES,'kmer_sex/k32.chrX.diff')
KMER_CHRM= pj(RESOURCES,'kmer_sex/k32.chrM.diff')
KMER_AUTO= pj(RESOURCES,'kmer_sex/k32.auto.diff')

#path to ref and add ref files
REF =  pj(RESOURCES, 'hg38_phix/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa')
REF_DIR = pj(RESOURCES, 'hg38_phix')

SHIFTED_MT= pj(RESOURCES,'MT_ref_shifted')
SHIFTED_MT_fa= pj(RESOURCES,'MT_ref_shifted/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta')
SHIFTED_MT_dict= pj(RESOURCES,'MT_ref_shifted/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict')
SHIFTED_MT_fai= pj(RESOURCES,'MT_ref_shifted/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai')
MT_CHAIN= pj(RESOURCES,'MT_ref_shifted/ShiftBack.chain')
ORIG_MT_fa= pj(RESOURCES,'MT_ref/chrM_hg38.fasta')
ORIG_MT_dict= pj(RESOURCES,'MT_ref/chrM_hg38.dict')
ORIG_MT_fai= pj(RESOURCES,'MT_ref/chrM_hg38.fasta.fai')
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


# kraken db
KRAKEN_DB = pj(RESOURCES, 'kraken/pluspf_20230605')

#tmp folders
TMPDIR = 'tmp' #do not use scratch, amount of storage is limited
TMPDIR_ALT = '/scratch-local'
tmpdir = pj(TMPDIR,getpass.getuser())
tmpdir_alternative = pj(TMPDIR_ALT,getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

current_dir = os.getcwd()
