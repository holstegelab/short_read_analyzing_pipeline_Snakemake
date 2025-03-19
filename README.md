
# Main aim of the pipeline

The main aim of the pipeline is to provide a simple and easy-to-use tool for analyzing WGS and WES data. 
The pipeline is designed to be used on the server with **SLURM** workload manager and **ZSLURM** add-on to it.
The pipeline is written in **Snakemake** and uses **Conda** for managing dependencies for portable and reproducible analysis.



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
  
* Run **GATK HaplotypeCaller** or **Deepvariant** (default) for producing gVCF
    > Deepvariant is faster and more accuarte than GATK according to in-house tests
    > Test results will be published soon
* Combine gVCFs in unified Database with **GLnexus** for the next *cohort calling*
  * Alternative version with **GenommicDBImport** or **Combinegvcf** for GATK-based pipeline followed by **Genotypegvcf**
* normaliztion of SNPs with pathced version of **bcftools norm** 
* get several statistics and producing combined statistic files in dir *stat/* and in tsv-files with extension *.bam_quality.tab*

## Additional features

 * SV detection with **delly**
 * CNV detection with **cnvkit**
 * CNV detection with **GATK** 

# HOW TO USE
1. clone this repo on server
2. *If you want use Zslurm*
    1. install Zslurm according to the manual page
    2. Open Zslurm
    3. run pipeline with **snakemake** command

      > snakemake --profile ~/.config/snakemake/zslurm/ --snakefile ~/short_read_analyzing_pipeline_Snakemake/Snakefile --use-conda --rerun-incomplete --retries 0 
      
    > **NOTE ABOUT PROFILE**
   > copy zslurm.yaml to ~/.config/snakemake/zslurm/config.yaml and change conda prefix to your conda prefix
    
4. If you want to run just several steps (for example only Aligment step) - 
choose suitible **smk** file as **--snakefile**

## Additional notes


