
# Main aim of the pipeline
This pipeline is made for analyzing short-read sequncing data. 
Analyzing include several steps:
* Removing adapters with **AdapterRemoval** 

  > AdapterRemoval doesn't delete sequences from reads, so we always will have a back-up copy 
  
* Aligning reads with **Dragmap**  
* Merge different bam-files after Dragmap for a single sample *only in case if one sample has more than 1 pair of fastq files*
* Mark duplicates with **Samtools**

  > at this step we create additional **CRAM** files for storage. This files will be saved on the tape and used as back-up copy
  
* Check contamination fraction with custom script

  > Custom scritp **bam_stats.py** calculate number of *secondary aligments*
  > If there are too much (more than 0.5%) another script **bam_clean.py** will remove contamination
  
* Run **GATK HaplotypeCaller** for producing gVCF
* Combine gVCFs in unified Database with **GATK GenomicsDBImport** for the next *cohort calling*
  * Possibility of combining different databases is on the roadmap
  * Possibility of adding samples to an existing database is under testing
  * Alternative version with **Combinegvcf** instead of **GenomicsDBImport** is under development
* Genotype merged database with **GATK Genotype**
* normaliztion of SNPs with **bcftools norm**
* VQSR with **GATK VQSR** and producing final vcf *Merged_after_VQSR_norm.vcf*
* get several statistics and producing combined statistic files in dir *stat/* and in tsv-files with extension *.bam_quality.tab*

## Additional features

 * SV detection with **delly**
 * CNV detection with **cnvkit**
 * CNV detection with **GATK** (under testing)

# HOW TO USE
1. clone this repo on server
2. *If you want use Zslurm*
    1. install Zslurm according to the manual page
    2. Open Zslurm
    3. run pipeline with **snakezcluster** command
  
      > example code
      > NB! edit command
      > snakezcluster --snakefile Snakefile
    
    4. If you want to run just several steps (for example only Aligment step) - 
    choose suitible **smk** file as **--snakefile**

## Additional notes

First step (**Aligner.smk**) require high amount of RAM and runs unstable on **PARTIAL** thin nodes. To avoid bugs please use **exclusive thin nodes** or **fat (part or exclusive)** node
You can run this step on *fat* node with option **--snakefile Aligner.smk** and after it's finish run all other steps on *thin* nodes with option **--snakefile Snakefile**

