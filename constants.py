import os
import getpass

#parameters
pj = os.path.join
RESOURCES = '/gpfs/work3/0/qtholstg/hg38_res/'

#region files
INTERVALS_DIR = pj(RESOURCES,'intervals_v2')
MERGED_CAPTURE_KIT_BED = pj(INTERVALS_DIR, 'merged_capture_kits_cds.bed')
MERGED_CAPTURE_KIT_IVL = pj(INTERVALS_DIR, 'merged_capture_kits_cds.interval_list')
TARGETS_BED = pj(INTERVALS_DIR,'gencode_43_cds.bed')
TARGETS_IVL = pj(INTERVALS_DIR,'gencode_43_cds.interval_list')
PL_PR_TABLE = pj(RESOURCES, 'ploydi_priors_table_hg38.tsv')

#resource folder with cram reference fasta files
CRAMREFS = pj(RESOURCES,'cram_refs')

# dir with this file
SNAKEMAKE_DIR_PATH = os.path.dirname('')

#conda env's paths
CONDA_VERIFYBAMID = 'envs/verifybamid.yaml'
CONDA_MAIN = 'envs/preprocess.yaml'
CONDA_MAIN_RUN = pj(RESOURCES, 'envs/preprocess.yaml')
CONDA_VCF = 'envs/vcf_handling.yaml'
CONDA_PYPY = 'envs/pypy.yaml'
CONDA_PYPY_RUN = pj(RESOURCES, 'envs/pypy.yaml')
CONDA_KMC = 'envs/kmc.yaml'
CONDA_KMC_RUN = pj(RESOURCES, 'envs/kmc.yaml')
CONDA_MOSDEPTH = 'envs/mosdepth.yaml'
CONDA_PCA = 'envs/PCA.yaml'
CONDA_GATK_CNV = pj(RESOURCES, 'software', 'gatk_4.4', 'build', 'gatkcondaenv.yml')

DEFAULT_JAVA_OPTIONS = ' -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4 '


#OUTPUT FOLDERS
SOURCEDIR= 'source'
SAMPLEINFODIR= 'sampleinfo'
FETCHDIR= 'fetch'
LOG= 'logs'
BENCH= 'benchmark'
BAM= 'bams'
GVCF= 'gvcf'
VCF= 'vcfs'
VCF_Final= 'Final_VCF'
STAT= 'stats'
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

#custom scripts (encapsulate in srcdir())
BAMMERGE= 'scripts/bam_merge.py'
BAMSTATS= 'scripts/bam_stats.py'
DECHIMER= 'scripts/bam_dechimer.py'
DECHIMER_THRESHOLD= 0.005
MERGEPHASE = 'scripts/merge_phasing.py'


#path to kmer files
KMER_CHRY= pj(RESOURCES,'kmer_sex/k32.chrY.diff')
KMER_CHRX= pj(RESOURCES,'kmer_sex/k32.chrX.diff')
KMER_CHRM= pj(RESOURCES,'kmer_sex/k32.chrM.diff')
KMER_AUTO= pj(RESOURCES,'kmer_sex/k32.auto.diff')

#path to ref and add ref files
REF =  pj(RESOURCES, 'Ref/GRCh38_full_analysis_set_plus_decoy_hla.fa')
REF_DIR = pj(RESOURCES, 'Ref')

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
str_ref= 'Ref_PhiX_Male/GRCh38.str.zip'

REF_FEMALE = pj(RESOURCES, 'Ref_PhiX_Female/GRCh38_full_analysis_set_plus_decoy_hla.fa')
REF_FEMALE_DIR = pj(RESOURCES, 'Ref_PhiX_Female/')
REF_FEMALE_BED = pj(RESOURCES, 'Ref_PhiX_Female/hg38_alt_mask.female.bed')
REF_FEMALE_STR = pj(RESOURCES, 'Ref_PhiX_Female/GRCh38_full_analysis_set_plus_decoy_hla.str.zip')
REF_FEMALE_DICT = pj(RESOURCES, 'Ref_PhiX_Female/GRCh38_full_analysis_set_plus_decoy_hla.dict')
MASK_FEMALE_BED = pj(RESOURCES, 'Ref_PhiX_Female/hg38_alt_mask.female.bed')
REF_FEMALE_HASH = pj(RESOURCES, 'Ref_PhiX_Female/hash_table.cfg')
REF_FEMALE_FAI = pj(RESOURCES, 'Ref_PhiX_Female/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai')

REF_MALE = pj(RESOURCES, 'Ref_PhiX_Male/GRCh38_full_analysis_set_plus_decoy_hla.fa')
REF_MALE_DIR = pj(RESOURCES, 'Ref_PhiX_Male/')
REF_MALE_BED = pj(RESOURCES, 'Ref_PhiX_Male/hg38_alt_mask.male.bed')
REF_MALE_DICT = pj(RESOURCES, 'Ref_PhiX_Male/GRCh38_full_analysis_set_plus_decoy_hla.dict')
MASK_MALE_BED = pj(RESOURCES, 'Ref_PhiX_Male/hg38_alt_mask.male.bed')
REF_MALE_HASH = pj(RESOURCES, 'Ref_PhiX_Male/hash_table.cfg')
REF_MALE_STR = pj(RESOURCES, 'Ref_PhiX_Male/GRCh38_full_analysis_set_plus_decoy_hla.str.zip')
REF_MALE_FAI = pj(RESOURCES, 'Ref_PhiX_Male/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai')

#verifybamid files
VERIFYBAMID_EXOME = pj(RESOURCES,'verifybamid_hg38_res/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat')
VERIFYBAMID_WGS = pj(RESOURCES, 'verifybamid_hg38_res/wgs/1000g.phase3.100k.b38.vcf.gz.dat')

# VQSR DBses
DBSNP = pj(RESOURCES,'databases/HG38_dbSNP_v155_updatesd.vcf.gz')
HAPMAP = pj(RESOURCES,'databases/hapmap_3.3.hg38.vcf.gz')
OMNI = pj(RESOURCES,'databases/1000G_omni2.5.hg38.vcf.gz')
KILO_G = pj(RESOURCES,'databases/1000G_phase1.snps.high_confidence.hg38.vcf.gz')
MILLS = pj(RESOURCES,'databases/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz')
DBSNP_INDEL = pj(RESOURCES,'databases/Homo_sapiens_assembly38.known_indels.vcf.gz')

#path to file with adapters
ADAPTERS = pj(RESOURCES, 'databases/Adapters_illumina.txt')

#windows
WINDOWS = pj(INTERVALS_DIR, 'windows/all.selected.sorted.3.bed')
WINDOWS_ANNOTATED = pj(INTERVALS_DIR, 'windows/all.selected.sorted.bed')

#tmp folders
TMPDIR = 'tmp' #do not use scratch, amount of storage is limited
TMPDIR_ALT = '/scratch-local'
tmpdir = pj(TMPDIR,getpass.getuser())
tmpdir_alternative = pj(TMPDIR_ALT,getpass.getuser())

os.makedirs(tmpdir,mode=0o700,exist_ok=True)

current_dir = os.getcwd()