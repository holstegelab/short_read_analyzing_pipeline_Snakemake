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
    output: directory(pj(INTERVALS_DIR,'wgs_bins_v3'))
    params:
        nsplit=1000,
        nproc=8,
        slop_script=srcdir(SLOPSCRIPT)
    conda: CONDA_MAIN
    resources:
        n="1.0",
        mem_mb=250
    shell: r"""
        mkdir -p {output}
        cd {output}
        bedtools sort -i {input.merged_kit} > {input.merged_kit}.sorted.bed
        awk '{{print $1 "\t" 0 "\t" $2}}' {input.fai} | awk '$1 != "chrom"' > genome.orig.bed
        bedtools sort -i genome.orig.bed > genome.sorted.bed
        bedtools multiinter -i genome.sorted.bed {input.merged_kit}.sorted.bed | cut -f1-3 > all_chroms_gencode.bed
        bedtools subtract -a all_chroms_gencode.bed -b {input.mask} > genome.masked.bed
        bedtools makewindows -b genome.masked.bed -w 10000 > genome.masked_windows.bed
        bedtools nuc -fi {REF_MALE} -bed genome.masked_windows.bed > genome.masked_windows_annot.bed
        awk '$10 < $12' genome.masked_windows_annot.bed | cut -f1-3  > genome.bed

        bedtools sort -i genome.bed -g {input.fai} > genome.sorted.bed
        touch genome.split
        rm genome.spli*
        cat genome.sorted.bed | grep -v -P "chrX|chrY" > genome.auto.bed
        cat genome.sorted.bed | grep -v -P "chrX|chrY|HLA|NC|chrUn|alt|random|KMT2C|EBV|MAP2K3|KCNJ18" > genome.auto_classic.bed
        cat genome.sorted.bed | grep chrX > genome.X.bed
        cat genome.sorted.bed | grep chrY > genome.Y.bed
        cat genome.sorted.bed | grep  -P "HLA|NC|chrUn|alt|random|KMT2C|EBV|MAP2K3|KCNJ18" > genome.O.bed

        #need the full contig for all the remaining contigs, as this is required for genomicsdbimport to merge intervals.
        cat genome.orig.bed | grep -P "HLA|NC|chrUn|alt|random|KMT2C|EBV|MAP2K3|KCNJ18" > genome.other.bed

        echo "Main files prepared, now splitting"

        split genome.auto.bed -a 4 -d -n l/10000 genome.autosplit4.
        split genome.X.bed -a 3 -d -n l/500 genome.Xsplit3.
        split genome.Y.bed -a 3 -d -n l/200 genome.Ysplit3.
        split genome.other.bed -a 3 -d -n l/200 genome.Osplit3.

        echo "Merging bed records"
        ls genome.*split?.* | xargs -I{{}} -P {params.nproc} sh -c "bedtools merge -i {{}} | awk '\$2 != \$3' > {{}}.bed"
        echo "Merging bed records, filtering on classic regions"
        ls genome.autosplit4.???? | xargs -I{{}} -P {params.nproc} sh -c "bedtools merge -i {{}} | awk '\$2 != \$3' | grep -v -P 'chrX|chrY|HLA|NC|chrUn|alt|random|KMT2C|EBV|MAP2K3|KCNJ18' > {{}}.classic.bed || true" || true

        echo "Splitting done, now dropping empty/temp files"
        
        #drop files that now do not contain any regions
        find . -name "genome.autosplit4.????.classic.bed" -type f -empty -print -delete
       
        echo "Empty classic files dropped"

        rm genome.*split4.????
        rm genome.*split3.???


        
        echo "Merging split files level 4"
        #level 3
        for j in $(seq 0 999); do
            counter=$(printf "%03d" $j)
            cat genome.autosplit4.$counter?.bed | bedtools merge > genome.autosplit3.$counter.bed
            (cat genome.autosplit4.$counter?.classic.bed 2>/dev/null || true) | bedtools merge > genome.autosplit3.$counter.classic.bed
        done
        find . -name "genome.autosplit3.*.classic.bed" -type f -empty -print -delete

        for j in $(seq 0 49); do
            counter=$(printf "%02d" $j)
            cat genome.Xsplit3.$counter?.bed | bedtools merge > genome.Xsplit2.$counter.bed
        done
        for j in $(seq 0 19); do
            counter=$(printf "%02d" $j)
            cat genome.Ysplit3.$counter?.bed | bedtools merge > genome.Ysplit2.$counter.bed
        done
        for j in $(seq 0 19); do
            counter=$(printf "%02d" $j)
            cat genome.Osplit3.$counter?.bed | bedtools merge > genome.Osplit2.$counter.bed
        done

        echo "Merging split files level 3"

        #level 2
        for j in $(seq 0 99); do
            counter=$(printf "%02d" $j)
            cat genome.autosplit3.$counter?.bed | bedtools merge > genome.autosplit2.$counter.bed
            (cat genome.autosplit3.$counter?.classic.bed 2>/dev/null || true) | bedtools merge > genome.autosplit2.$counter.classic.bed
        done
        find . -name "genome.autosplit2.*.classic.bed" -type f -empty -print -delete

        for i in $(seq 0 4); do
            cat genome.Xsplit2.$i?.bed | bedtools merge > genome.Xsplit1.$i.bed
        done
        for i in $(seq 0 1); do
            cat genome.Ysplit2.$i?.bed | bedtools merge > genome.Ysplit1.$i.bed
        done

        for i in $(seq 0 1); do
            cat genome.Osplit2.$i?.bed | bedtools merge > genome.Osplit1.$i.bed
        done


        echo "Merging split files level 2"
        #level 1
        for i in $(seq 0 9); do
            cat genome.autosplit2.$i?.bed | bedtools merge > genome.autosplit1.$i.bed
            cat genome.autosplit2.$i?.classic.bed | bedtools merge > genome.autosplit1.$i.classic.bed
        done

        echo "Adding slop to level 1"

        ls genome.autosplit1.?.bed | xargs -I{{}} -P {params.nproc} sh -c "python {params.slop_script} {{}} {input.fai}"
        ls genome.autosplit1.?.classic.bed | xargs -I{{}} -P {params.nproc} sh -c "python {params.slop_script} {{}} {input.fai}"
        echo "Adding slop to level 2"
        ls genome.autosplit2.??.bed | xargs -I{{}} -P {params.nproc} sh -c "python {params.slop_script} {{}} {input.fai}"
        ls genome.autosplit2.??.classic.bed | xargs -P {params.nproc} -I{{}} sh -c "python {params.slop_script} {{}} {input.fai}"
        ls genome.Xsplit1.?.bed | xargs -P {params.nproc} -I{{}} sh -c "python {params.slop_script} {{}} {input.fai}"
        ls genome.Ysplit1.?.bed | xargs -P {params.nproc} -I{{}} sh -c "python {params.slop_script} {{}} {input.fai}"
        echo "Adding slop to level 3"
        ls genome.autosplit3.???.bed | xargs -I{{}} -P {params.nproc} sh -c "python {params.slop_script} {{}} {input.fai}"
        ls genome.autosplit3.???.classic.bed | xargs -P {params.nproc} -I{{}} sh -c "python {params.slop_script} {{}} {input.fai}"
        ls genome.Xsplit2.??.bed | xargs -P {params.nproc} -I{{}} sh -c "python {params.slop_script} {{}} {input.fai}"
        ls genome.Ysplit2.??.bed | xargs -P {params.nproc} -I{{}} sh -c "python {params.slop_script} {{}} {input.fai}"

        echo "Merging split files level 1"
        cat genome.sorted.bed | bedtools merge > genome.fullsplit0.bed
        cat genome.sorted.bed | grep -v -P "HLA|NC|chrUn|alt|random|KMT2C|EBV|MAP2K3|KCNJ18" > genome.fullsplit0.classic.bed
        cat genome.autosplit1.?.bed | bedtools merge > genome.autosplit0.bed
        cat genome.autosplit1.?.classic.bed | bedtools merge > genome.autosplit0.classic.bed
        cat genome.Xsplit1.?.bed | bedtools merge > genome.Xsplit0.bed
        cat genome.Ysplit1.?.bed | bedtools merge > genome.Ysplit0.bed
                """


rule CreateBinsExome:
    """Create bins for the exome. Separate bins for X and Y.
    1. Determine merged kit padding regions (300bp)
    2. Intersect for each level/split in CreateBinsFullGenome the padded merged kit file

    """
    input: merged_kit=ancient(MERGED_CAPTURE_KIT_BED),
        wgs_folder=pj(INTERVALS_DIR,'wgs_bins_v3'),
        fai=ancient(REF_MALE_FAI),
    output: directory(pj(INTERVALS_DIR,'wes_bins_v3'))
    params:
        nsplit=1000
    conda: CONDA_MAIN
    resources:
        n="1.0",
        mem_mb=250
    shell: r"""
        mkdir -p {output}
        cd {output}
        bedtools slop -i {input.merged_kit} -g {input.fai} -b 300 | bedtools merge | bedtools sort > {input.merged_kit}.padded.bed

        #level 3
        for j in $(seq 0 999); do
            counter=$(printf "%03d" $j)
            bedtools intersect -a {input.wgs_folder}/genome.autosplit3.$counter.bed -b {input.merged_kit}.padded.bed > merged.autosplit3.$counter.bed
        done
        
        for j in $(seq 0 999); do
            counter=$(printf "%03d" $j)
            if [ -f {input.wgs_folder}/genome.autosplit3.$counter.classic.bed ]; then
                bedtools intersect -a {input.wgs_folder}/genome.autosplit3.$counter.classic.bed -b {input.merged_kit}.padded.bed > merged.autosplit3.$counter.classic.bed
            fi
        done

        for j in $(seq 0 49); do
            counter=$(printf "%02d" $j)
            bedtools intersect -a {input.wgs_folder}/genome.Xsplit2.$counter.bed -b {input.merged_kit}.padded.bed > merged.Xsplit2.$counter.bed
            bedtools intersect -a {input.wgs_folder}/genome.Xsplit2.$counter.padded.bed -b {input.merged_kit}.padded.bed > merged.Xsplit2.$counter.padded.bed
        done
        for j in $(seq 0 19); do
            counter=$(printf "%02d" $j)
            bedtools intersect -a {input.wgs_folder}/genome.Ysplit2.$counter.bed -b {input.merged_kit}.padded.bed > merged.Ysplit2.$counter.bed
            bedtools intersect -a {input.wgs_folder}/genome.Ysplit2.$counter.padded.bed -b {input.merged_kit}.padded.bed > merged.Ysplit2.$counter.padded.bed
        done

        for j in $(seq 0 19); do
            counter=$(printf "%02d" $j)
            bedtools intersect -a {input.wgs_folder}/genome.Osplit2.$counter.bed -b {input.merged_kit}.padded.bed > merged.Osplit2.$counter.bed
        done


        #level 2
        for j in $(seq 0 99); do
            counter=$(printf "%02d" $j)
            bedtools intersect -a {input.wgs_folder}/genome.autosplit2.$counter.bed -b {input.merged_kit}.padded.bed > merged.autosplit2.$counter.bed
            bedtools intersect -a {input.wgs_folder}/genome.autosplit2.$counter.padded.bed -b {input.merged_kit}.padded.bed > merged.autosplit2.$counter.padded.bed
        done
        for j in $(seq 0 99); do
            counter=$(printf "%02d" $j)
            if [ -f {input.wgs_folder}/genome.autosplit2.$counter.classic.bed ]; then
                bedtools intersect -a {input.wgs_folder}/genome.autosplit2.$counter.classic.bed -b {input.merged_kit}.padded.bed > merged.autosplit2.$counter.classic.bed
                bedtools intersect -a {input.wgs_folder}/genome.autosplit2.$counter.classic.padded.bed -b {input.merged_kit}.padded.bed > merged.autosplit2.$counter.classic.padded.bed
            fi
        done

        for i in $(seq 0 4); do
            bedtools intersect -a {input.wgs_folder}/genome.Xsplit1.$i.bed -b {input.merged_kit}.padded.bed > merged.Xsplit1.$i.bed
            bedtools intersect -a {input.wgs_folder}/genome.Xsplit1.$i.padded.bed -b {input.merged_kit}.padded.bed > merged.Xsplit1.$i.padded.bed
        done
        for i in $(seq 0 1); do
            bedtools intersect -a {input.wgs_folder}/genome.Ysplit1.$i.bed -b {input.merged_kit}.padded.bed > merged.Ysplit1.$i.bed
            bedtools intersect -a {input.wgs_folder}/genome.Ysplit1.$i.padded.bed -b {input.merged_kit}.padded.bed > merged.Ysplit1.$i.padded.bed
        done
        for i in $(seq 0 1); do
            bedtools intersect -a {input.wgs_folder}/genome.Osplit1.$i.bed -b {input.merged_kit}.padded.bed > merged.Osplit1.$i.bed
        done

        #level 1
        for i in $(seq 0 9); do
            bedtools intersect -a {input.wgs_folder}/genome.autosplit1.$i.bed -b {input.merged_kit}.padded.bed > merged.autosplit1.$i.bed
            bedtools intersect -a {input.wgs_folder}/genome.autosplit1.$i.padded.bed -b {input.merged_kit}.padded.bed > merged.autosplit1.$i.padded.bed
        done

        for i in $(seq 0 9); do
            bedtools intersect -a {input.wgs_folder}/genome.autosplit1.$i.classic.bed -b {input.merged_kit}.padded.bed > merged.autosplit1.$i.classic.bed
            bedtools intersect -a {input.wgs_folder}/genome.autosplit1.$i.classic.padded.bed -b {input.merged_kit}.padded.bed > merged.autosplit1.$i.classic.padded.bed
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


rule precompute_bed_for_capture_kit_checker:
    input: kit_dir = INTERVALS_DIR,
            genome = GENOME_FILE
    output: precomputed = PRECOMPUTEED_BED
    params: CAPTURE_KIT_CHECKER=srcdir(CAPTURE_KIT_CHECKER)
    conda: CONDA_CK_FINDER
    resources:
        n="16",
        mem_mb=8000
    shell:
        """
        python {params.CAPTURE_KIT_CHECKER}  --precompute_mode \
            --kit_dir {input.kit_dir} \
            --genome {input.genome} \
            --precompute_output {output.precomputed} \
            --threads {resources.n} \
        """

