import os
import getpass

#parameters
pj = os.path.join
RESOURCES = '/gpfs/work3/0/qtholstg/hg38_res_v2/'

#region files
INTERVALS_DIR = pj(RESOURCES,'intervals')
MERGED_CAPTURE_KIT_BED = pj(INTERVALS_DIR, 'merged_capture_kits_cds.bed')
MERGED_CAPTURE_KIT_IVL = pj(INTERVALS_DIR, 'merged_capture_kits_cds.interval_list')
TARGETS_BED = pj(INTERVALS_DIR,'gencode_43_cds.bed')
TARGETS_IVL = pj(INTERVALS_DIR,'gencode_43_cds.interval_list')

#resource folder with cram reference fasta files
CRAMREFS = pj(RESOURCES,'cram_refs')

#conda env's paths
CONDA_VERIFYBAMID = 'envs/verifybamid.yaml'
CONDA_MAIN = 'envs/preprocess.yaml'
CONDA_VCF = 'envs/vcf_handling.yaml'
CONDA_PYPY = 'envs/pypy.yaml'
CONDA_KMC = 'envs/kmc.yaml'
CONDA_KRAKEN = 'envs/kraken.yaml'
CONDA_MOSDEPTH = 'envs/mosdepth.yaml'
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

#programs
SOFTWARE = pj(RESOURCES, 'software')
gatk= 'gatk'
samtools= 'samtools'
bcftools= 'bcftools'
dragmap= 'dragen-os'
verifybamid2= 'verifybamid2'

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
