from common import *


rule BedToIntervalList:
    input:
        bed_file="{path}.bed"
    output:
        interval_list="{path}.interval_list"
    params:
        seq_dict=REF_MALE_DICT
    conda: CONDA_VCF
    resources:
        n="1.0",
        mem_mb=1000
    shell:
        """{gatk} BedToIntervalList  --java-options "-Xmx{resources.mem_mb}m {DEFAULT_JAVA_OPTIONS}"  -I {input.bed_file} -O {output.interval_list} -SD {params.seq_dict} --UNIQUE"""


rule CreateBinsFullGenome:
    """Create bins for the full genome. Separate bins for X and Y.
    1. Create full genome bed file
    2. Cut regions according to merged capture kit
    3. Subtract masked regions of dragmap
    4. Create 10kb bins within bed records
    5. Annotate bins with number of N nucleotides, and remove bins that have only N nucleotides
    6. Split bins into 10000 bins for autosomes, 500 bins for X, 200 bins for Y (level 4)
    7. Create level 3 bins by merging level 4 bins (bin 000 = 0000-0009, etc.)
    8. Create level 2 bins by merging level 3 bins (bin 00 = 000-0099, etc.)
    9. Create level 1 bins by merging level 2 bins (bin 0 = 0000-0999, etc.)
    10. Create full genome bin by merging level 1 bins (bin 0 = 0000-9999, etc.)

    The 1070 level 3 bins are used for genotyping over all samples.
    The 12 level 1 bins (10 auto + x + y) are used for variant calling per sample.

    """
    input: fai=ancient(REF_MALE_FAI),
        mask=ancient(REF_MALE_BED),
        merged_kit=ancient(MERGED_CAPTURE_KIT_BED)
    output: directory(pj(INTERVALS_DIR,'wgs_bins'))
    params:
        nsplit=1000
    conda: CONDA_MAIN
    resources:
        n="1.0",
        mem_mb=250
    shell: """
        mkdir -p {output}
        cd {output}
        bedtools sort -i {input.merged_kit} > {input.merged_kit}.sorted.bed
        awk '{{print $1 "\t" 0 "\t" $2}}' {input.fai} | awk '$1 != "chrom"' > genome.bed
        bedtools sort -i genome.bed > genome.sorted.bed
        bedtools multiinter -i genome.sorted.bed {input.merged_kit}.sorted.bed | cut -f1-3 > all_chroms_gencode.bed
        bedtools subtract -a all_chroms_gencode.bed -b {input.mask} > genome.masked.bed
        bedtools makewindows -b genome.masked.bed -w 10000 > genome.masked_windows.bed
        bedtools nuc -fi {REF_MALE} -bed genome.masked_windows.bed > genome.masked_windows_annot.bed
        awk '$10 < $12' genome.masked_windows_annot.bed | cut -f1-3  > genome.bed

        bedtools sort -i genome.bed -g {input.fai} > genome.sorted.bed
        touch genome.split
        rm genome.spli*
        cat genome.sorted.bed | grep -v chrX | grep -v chrY > genome.auto.bed
        cat genome.sorted.bed | grep chrX > genome.X.bed
        cat genome.sorted.bed | grep chrY > genome.Y.bed

        split genome.auto.bed -a 4 -d -n l/10000 genome.autosplit4.
        split genome.X.bed -a 3 -d -n l/500 genome.Xsplit3.
        split genome.Y.bed -a 3 -d -n l/200 genome.Ysplit3.

        ls genome.*split?.* | xargs -I{{}} sh -c "bedtools merge -i {{}} | awk '\$2 != \$3' > {{}}.bed"
        rm genome.*split4.????
        rm genome.*split3.???

        #level 3
        for j in $(seq 0 999); do
            counter=$(printf "%03d" $j)
            cat genome.autosplit4.$counter?.bed | bedtools merge > genome.autosplit3.$counter.bed
        done

        for j in $(seq 0 49); do
            counter=$(printf "%02d" $j)
            cat genome.Xsplit3.$counter?.bed | bedtools merge > genome.Xsplit2.$counter.bed
        done
        for j in $(seq 0 19); do
            counter=$(printf "%02d" $j)
            cat genome.Ysplit3.$counter?.bed | bedtools merge > genome.Ysplit2.$counter.bed
        done


        #level 2
        for j in $(seq 0 99); do
            counter=$(printf "%02d" $j)
            cat genome.autosplit3.$counter?.bed | bedtools merge > genome.autosplit2.$counter.bed
        done
        for i in $(seq 0 4); do
            cat genome.Xsplit2.$i?.bed | bedtools merge > genome.Xsplit1.$i.bed
        done
        for i in $(seq 0 1); do
            cat genome.Ysplit2.$i?.bed | bedtools merge > genome.Ysplit1.$i.bed
        done

        #level 1
        for i in $(seq 0 9); do
            cat genome.autosplit2.$i?.bed | bedtools merge > genome.autosplit1.$i.bed
        done

        cat genome.sorted.bed | bedtools merge > genome.fullsplit0.bed
        cat genome.autosplit1.?.bed | bedtools merge > genome.autosplit0.bed
        cat genome.Xsplit1.?.bed | bedtools merge > genome.Xsplit0.bed
        cat genome.Ysplit1.?.bed | bedtools merge > genome.Ysplit0.bed
                """


rule CreateBinsExome:
    """Create bins for the exome. Separate bins for X and Y.
    1. Determine merged kit padding regions (300bp)
    2. Intersect for each level/split in CreateBinsFullGenome the padded merged kit file

    """
    input: merged_kit=ancient(MERGED_CAPTURE_KIT_BED),
        wgs_folder=pj(INTERVALS_DIR,'wgs_bins'),
        fai=ancient(REF_MALE_FAI),
    output: directory(pj(INTERVALS_DIR,'wes_bins'))
    params:
        nsplit=1000
    conda: CONDA_MAIN
    resources:
        n="1.0",
        mem_mb=250
    shell: """
        mkdir -p {output}
        cd {output}
        bedtools slop -i {input.merged_kit} -g {input.fai} -b 300 | bedtools merge | bedtools sort > {input.merged_kit}.padded.bed

        #level 3
        for j in $(seq 0 999); do
            counter=$(printf "%03d" $j)
            bedtools intersect -a {input.wgs_folder}/genome.autosplit3.$counter.bed -b {input.merged_kit}.padded.bed > merged.autosplit3.$counter.bed
        done

        for j in $(seq 0 49); do
            counter=$(printf "%02d" $j)
            bedtools intersect -a {input.wgs_folder}/genome.Xsplit2.$counter.bed -b {input.merged_kit}.padded.bed > merged.Xsplit2.$counter.bed
        done
        for j in $(seq 0 19); do
            counter=$(printf "%02d" $j)
            bedtools intersect -a {input.wgs_folder}/genome.Ysplit2.$counter.bed -b {input.merged_kit}.padded.bed > merged.Ysplit2.$counter.bed
        done



        #level 2
        for j in $(seq 0 99); do
            counter=$(printf "%02d" $j)
            bedtools intersect -a {input.wgs_folder}/genome.autosplit2.$counter.bed -b {input.merged_kit}.padded.bed > merged.autosplit2.$counter.bed
        done
        for i in $(seq 0 4); do
            bedtools intersect -a {input.wgs_folder}/genome.Xsplit1.$i.bed -b {input.merged_kit}.padded.bed > merged.Xsplit1.$i.bed
        done
        for i in $(seq 0 1); do
            bedtools intersect -a {input.wgs_folder}/genome.Ysplit1.$i.bed -b {input.merged_kit}.padded.bed > merged.Ysplit1.$i.bed
        done

        #level 1
        for i in $(seq 0 9); do
            bedtools intersect -a {input.wgs_folder}/genome.autosplit1.$i.bed -b {input.merged_kit}.padded.bed > merged.autosplit1.$i.bed
        done

        bedtools intersect -a {input.wgs_folder}/genome.autosplit0.bed -b {input.merged_kit}.padded.bed > merged.autosplit0.bed
        bedtools intersect -a {input.wgs_folder}/genome.fullsplit0.bed -b {input.merged_kit}.padded.bed > merged.fullsplit0.bed
        bedtools intersect -a {input.wgs_folder}/genome.Xsplit0.bed -b {input.merged_kit}.padded.bed > merged.Xsplit0.bed
        bedtools intersect -a {input.wgs_folder}/genome.Ysplit0.bed -b {input.merged_kit}.padded.bed > merged.Ysplit0.bed
                """


rule select_bed_chrom:
    input:
        interval_list="{folder}/{capture_kit}.bed"
    output:
        interval_list_chrom="{folder}/{capture_kit}/{capture_kit}_chrom_{chrom}.bed"
    resources:
        n="1.0",
        mem_mb=1000
    shell: """
            grep -P '^{wildcards.chrom}\t' {input} > {output}
          """  
