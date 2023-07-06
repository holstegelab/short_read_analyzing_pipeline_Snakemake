configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")
gatk = config['gatk']
from common import *



rule BedToIntervalList:
    input:
        bed_file = "{path}.bed"
    output:
        interval_list = "{path}.interval_list"    
    params:
        java_options=config['DEFAULT_JAVA_OPTIONS'],
        seq_dict = os.path.join(config['RES'], config['ref_male_dict'])
    conda: "envs/vcf_handling.yaml"   
    resources: 
        n = 1,
        mem_mb = 1000       
    shell:
        """{gatk} BedToIntervalList  --java-options "-Xmx{resources.mem_mb}m {params.java_options}"  -I {input.bed_file} -O {output.interval_list} -SD {params.seq_dict} --UNIQUE"""


rule BedSplit:
    input:
        bed_file = "{folder}/{capture_kit}_chrom_{chrom}.bed"
    output:
        bed_file_folder = directory("{folder}/{capture_kit}_chrom_{chrom}.split")
    params: 
        nsplit = lambda wildcards: chromosome_splits[wildcards.chrom]
    shell: """
        mkdir -p {output.bed_file_folder}
        split {input.bed_file} -a 3 -d -n l/{params.nsplit} {output.bed_file_folder}/{capture_kit}_chrom_{chrom}.split.
        """

rule CreateBins:
    input: fai=ancient(os.path.join(config['RES'], config['ref_male_fai'])),
           mask=ancient(os.path.join(config['RES'], config['mask_male_bed'])),
           targets=ancient(os.path.join(config['RES'], config['kit_folder'], config['TARGETS'] + '.bed'))
    output: directory(os.path.join(config['RES'], config['kit_folder'], 'bins'))
    params:
        nsplit = 1000
    conda: "envs/preprocess.yaml"
    shell:"""
        mkdir -p {output}
        cd {output}
        bedtools slop -i {input.targets}  -g {input.fai} -b 100 | bedtools merge > gencode_43_cds.100padding.bed
        awk '{print $1 "\t" 0 "\t" $2}' {input.fai} | grep -v chrom > all_chroms.bed
        bedtools multiinter -i all_chroms.bed gencode_100padding.bed | cut -f1-3 > all_chroms_gencode.bed
        bedtools subtract -a all_chroms_gencode.bed -b {input.mask} > all_chroms.masked.bed
        bedtools makewindows -b all_chroms.masked.bed -w 10000 > all_chroms.masked_windows.bed
        bedtools nuc -fi ../../Ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -bed all_chroms.masked_windows.bed > all_chroms.masked_windows_annot.bed
        awk '$10 < $12' all_chroms.masked_windows_annot.bed | cut -f1-3  > all_chroms.masked_windows_annot_selected.bed

        touch all_chroms.split
        rm all_chroms.spli*

        split all_chroms.masked_windows_annot_selected.bed -a 4 -d -n l/{params.nsplit} all_chroms.split.

        ls all_chroms.split.* | xargs -I{} sh -c 'bedtools merge -i {} > {}.bed'
        rm all_chroms.split.????
        """


rule select_bed_chrom:
    input:
        interval_list = "{folder}/{capture_kit}.bed"
    output:
        interval_list_chrom = "{folder}/{capture_kit}/{capture_kit}_chrom_{chrom}.bed" 
    resources:
        n=1,
        mem_mb=1000        
    shell:"""
            grep -P '^{wildcards.chrom}\t' {input} > {output}
          """          
