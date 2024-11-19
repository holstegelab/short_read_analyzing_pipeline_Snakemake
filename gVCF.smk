from common import *
onsuccess: shell("rm -fr logs/gVCF/*")

wildcard_constraints:
    sample=r"[\w\d_\-@]+",

module Tools:
    snakefile:
        "Tools.smk"
    config:
        config

use rule * from Tools


rule gVCF_all:
    input:
        expand("{gvcf}/{sample}.done", sample=sample_names, gvcf=GVCF),
    default_target: True


def get_gvcf_files(wildcards):  # {{{
    sample = wildcards["sample"]
    if 'wgs' in SAMPLEINFO[sample]['sample_type']:
        return [pj(GVCF, "exome_extract", region, f"{sample}.{region}.wg.vcf.gz") for region in level1_regions]
    else:
        return [pj(GVCF, "reblock", region, f"{sample}.{region}.wg.vcf.gz") for region in level0_regions]
# }}}


rule gvcf_sample_done:
    input:
        get_gvcf_files,
    output:
        touch(pj(GVCF, "{sample}.done")),
    resources:
        n="1.0",
        mem_mb=50,


rule ComposeSTRTableFile:
    input:
        ref_fasta="{path}.fa",
    output:
        str_table="{path}.str.zip",
    params:
        java_options=DEFAULT_JAVA_OPTIONS,
    resources:
        n="1.0",
        mem_mb=3000,
    conda:
        CONDA_VCF
    shell:
        """gatk ComposeSTRTableFile --java-options "-Xmx{resources.mem_mb}m {params.java_options}" -R {input.ref_fasta} -O {output.str_table}"""


rule CalibrateDragstrModel:
    """CalibrateDragstrModel. Estimates the parameters of the dragstr model from a set of aligned reads.
    Provides better STR variant calling.
    """
    input:
        bam=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        validated_sex=pj(KMER,"{sample}.result.yaml"),
    output:
        dragstr_model=pj(BAM, "{sample}-dragstr.txt"),
    priority: 26
    params:
        ref=get_ref_by_validated_sex,
        str_ref=get_strref_by_validated_sex,
        java_options=DEFAULT_JAVA_OPTIONS,
    resources:
        n="1.0",
        mem_mb=lambda wildcards, attempt: attempt * 1750,
    log:
        pj(LOG, "gVCF","{sample}_calibratedragstr.log"),
    conda:
        CONDA_VCF
    shell:
        """{gatk} CalibrateDragstrModel --java-options \
                    "-Xmx{resources.mem_mb}m {params.java_options}"  -R {params.ref} -I {input.bam} \
                    -O {output} -str {params.str_ref} 2>{log}"""


def read_contam_w(wildcards):  # {{{
    """Read contamination from verifybamid output file."""
    filename = pj(STAT, "contam", wildcards["sample"] + ".verifybamid.pca2.selfSM")
    with open(filename, "r", encoding="utf-8") as f:
        c = csv.reader(f, delimiter="\t")
        c = iter(c)
        headers = c.__next__()
        data = c.__next__()
        freemix = float(data[6])

    if freemix < 0.01:  # only remove contamination if it is more than 1%
        freemix = 0.0
    return freemix
# }}}


def get_mem_mb_HaplotypeCaller(wildcards, attempt):  # {{{
    """Get memory for HaplotypeCaller."""
    res = 2800 if "wgs" in SAMPLEINFO[wildcards["sample"]]["sample_type"] else 1800

    # HaplotypeCaller has exponential memory scaling on some regions. Average is very  low, but on some regions it can scale to 75 GB...
    return res * (3 ** (attempt - 1))  # aggressively reserve more memory
# }}}


def region_to_interval_file(wildcards):  # {{{
    """Converts a region to a interval file location (see common.py and Tools.smk)"""
    sample = wildcards["sample"]
    region = wildcards["region"]
    return region_to_file(
        region,
        wgs="wgs" in SAMPLEINFO[sample]["sample_type"],
        padding=True,
        extension="interval_list",
    )
# }}}

rule HaplotypeCaller:
    """HaplotypeCaller. Call SNPs and indels for each sample."""
    input:
        # check what bam file we need to use (with or without additional cleanup)
        bams=pj(BAM, "{sample}.markdup.bam"),
        bai=pj(BAM, "{sample}.markdup.bam.bai"),
        model=rules.CalibrateDragstrModel.output.dragstr_model,
        contam=pj(STAT, 'contam/{sample}.verifybamid.pca2.selfSM'),
        interval=region_to_interval_file,
        validated_sex=pj(KMER,"{sample}.result.yaml"),
    output:
        orig_gvcf=temp(pj(GVCF, "raw/{region}/{sample}.{region}.g.vcf.gz")),
        orig_gvcf_tbi=temp(pj(GVCF, "raw/{region}/{sample}.{region}.g.vcf.gz.tbi")),
        genotyped_vcf=temp(pj(GVCF, "raw/{region}/{sample}.{region}.vcf.gz")),
        genotyped_vcf_tbi=temp(pj(GVCF, "raw/{region}/{sample}.{region}.vcf.gz.tbi")),
        tmp_vcf=temp(pj(GVCF, "raw/{region}/{sample}.{region}.tmp.vcf.gz")),
        vcf=temp(pj(GVCF, "raw/{region}/{sample}.{region}.w.vcf.gz")),
        vcf_tbi=temp(pj(GVCF, "raw/{region}/{sample}.{region}.w.vcf.gz.tbi")),
        wstats=pj(STAT, "whatshap_phasing/{sample}.{region}.stats"),
        mwstats=pj(STAT, "whatshap_phasing/{sample}.{region}.merge_stats"),
        tmp_gvcf=temp(pj(GVCF, "raw/{region}/{sample}.{region}.wg.vcf")),        
        gvcf=temp(pj(GVCF, "raw/{region}/{sample}.{region}.wg.vcf.gz")),
        gvcf_tbi=temp(pj(GVCF, "raw/{region}/{sample}.{region}.wg.vcf.gz.tbi")),
    conda:
        CONDA_VCF
    resources:
        n="1.1",  #average 1.3 cores
        mem_mb=get_mem_mb_HaplotypeCaller,
        tmpdir=tmpdir_alternative,
    params:
        dbsnp=DBSNP,
        padding=0,  #padding is included in the interval file. Keep at 0 to prevent overlapping intervals
        contam_frac=read_contam_w,  # get contamination fraction per sample
        ploidy=lambda wildcards: 1 if wildcards["region"].endswith("H") else 2,
        java_options=DEFAULT_JAVA_OPTIONS,
        ref=get_ref_by_validated_sex,
        dragen_mode=lambda wildcards: "--dragen-mode true " if not "H" in wildcards["region"] else "",
        merge_script=srcdir(MERGEPHASE),
        skipsex = lambda wildcards, input: int(get_validated_sex_file(input) == 'female' and wildcards['region'].startswith('Y')),
    priority: 28
    # Overview of available annotations (GATK version 4.4)
    # AS_StandardAnnotation: AS_BaseQualityRankSumTest, AS_FisherStrand, AS_InbreedingCoeff, AS_MappingQualityRankSumTest, AS_QualByDepth, AS_RMSMappingQuality, AS_ReadPosRankSumTest, AS_StrandOddsRatio (all allele specific annotations)
    # StandardHCAnnotation: DepthPerSampleHC
    # StandardANnotation: Coverage, ChromosomeCounts, BaseQualityRankSumTest, MappingQualityRankSumTest, FisherStrand, QualByDepth, InbreedingCoeff, DepthPerAlleleBySample, ReadPosRankSumTest, ExcessHet
    #                      StrandOddsRatio, RMSMappingQuality,
    # not included in standard annotation: FragmentDepthPerAlleleBySample,StrandBiasBySample, AssemblyComplexity,  FragmentLength, OrientationBiasReadCounts, TandemRepeat, UniqueAltReadCount

    # findings:
    # FragmentDepthPerAlleleBySample, OrientationBiasReadCounts, TandemRepeat --> do not output anything in haplotypecaller, omitted

    # regarding linked-de-bruijn-graph, see comments in https://gatk.broadinstitute.org/hc/en-us/community/posts/360077647812-Why-do-a-clear-expected-variant-not-show-up-in-the-Mutect2-vcf-file
    # also: https://gatk.broadinstitute.org/hc/en-us/articles/360043491652-When-HaplotypeCaller-and-Mutect2-do-not-call-an-expected-variant
    # still, seems quite experimental, so not using it for now

    # related: --recover-all-dangling-branches
    # to consider: --population-callset in combination with gnomad,  --use-posteriors-to-calculate-qual

    # max-mnp-distance: not compatible with multiple sample gVCF calling, omitted.

    # PHASING
    # GATK: physical phasing depends on de bruijn graph, which is not always optimal. In practice this means
    # that the phasing is only short distance (~50bp). Whathap uses the full dna fragment, usually 300bp, to phase.
    # Here we use whatshap to phase the gatk called variants, and then merge the whatshap phased variants with the gatk phased variants.
    # This way we get the best of both worlds: the gatk phasing is used for short distance phasing also of variants with low read
    # support (which might become full variants during multi-sample genotyping), and the whatshap phasing is used for long distance phasing.

    shell:
        """ 
        {gatk} --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" HaplotypeCaller \
            -R {params.ref} -L {input.interval} -ip {params.padding} -D {params.dbsnp} -ERC GVCF --contamination {params.contam_frac} \
            --ploidy {params.ploidy} -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
            --annotate-with-num-discovered-alleles --adaptive-pruning \
            -A StrandBiasBySample -A AssemblyComplexity -A FragmentLength \
            -I {input.bams} -O {output.orig_gvcf}  --native-pair-hmm-threads 2  --create-output-variant-index true \
            --seconds-between-progress-updates 120 \
            {params.dragen_mode} --dragstr-params-path {input.model}

        {gatk} --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" GenotypeGVCFs \
                -R {params.ref} -V {output.orig_gvcf} -O {output.genotyped_vcf} \
        --seconds-between-progress-updates 120 


        mkdir -p `dirname {output.wstats}`
        if [ {params.ploidy} -eq 2 ]
        then
            if [ {params.skipsex} -eq 0 ]
            then 
                whatshap unphase  {output.genotyped_vcf} >  {output.tmp_vcf}
                whatshap phase  --ignore-read-groups --reference {params.ref} {output.tmp_vcf} {input.bams} -o {output.vcf} |tail -n 20
                bcftools index --tbi {output.vcf}
                whatshap stats {output.genotyped_vcf} > {output.wstats} 
                whatshap stats {output.vcf} >> {output.wstats} 
                python {params.merge_script} {output.orig_gvcf} {output.vcf} {output.tmp_gvcf} {output.mwstats} -q
                bcftools annotate -x 'FORMAT/PGT,FORMAT/PS,FORMAT/PID' {output.tmp_gvcf} -o {output.gvcf}
            else
                touch {output.tmp_vcf}
                touch {output.vcf}
                touch {output.vcf}.tbi
                whatshap stats {output.genotyped_vcf} > {output.wstats} || true
                touch {output.mwstats}
                touch {output.tmp_gvcf}
                bcftools view {output.orig_gvcf} -o {output.gvcf}
            fi
        else
            touch {output.tmp_vcf}
            touch {output.vcf}
            touch {output.vcf}.tbi
            whatshap stats {output.genotyped_vcf} > {output.wstats} || true
            touch {output.mwstats}
            touch {output.tmp_gvcf}
            bcftools view {output.orig_gvcf} -o {output.gvcf}
        fi  
        bcftools index --tbi {output.gvcf}      
        """


rule reblock_gvcf:
    input:
        gvcf=rules.HaplotypeCaller.output.gvcf,
        idx=rules.HaplotypeCaller.output.gvcf_tbi

    output:
        gvcf_reblock=ensure( pj(GVCF, "reblock/{region}/{sample}.{region}.wg.vcf.gz"), non_empty=True),
        tbi=ensure( pj(GVCF, "reblock/{region}/{sample}.{region}.wg.vcf.gz.tbi"), non_empty=True),
    conda:
        CONDA_VCF
    priority: 29
    params:
        dbsnp=DBSNP,
        java_options=DEFAULT_JAVA_OPTIONS
    resources:
        n="1.0",
        mem_mb=lambda wildcards, attempt: attempt * 1250,
    shell:
        """
    {gatk} --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" ReblockGVCF  \
        --keep-all-alts --create-output-variant-index true -D {params.dbsnp} -R {REF_MALE} --do-qual-score-approximation \
         -V {input.gvcf} -O {output.gvcf_reblock}  --seconds-between-progress-updates 120 \
        -GQB 3 -GQB 5 -GQB 8 -GQB 10 -GQB 15 -GQB 20 -GQB 30 -GQB 50 -GQB 70 -GQB 100 \
        -G StandardAnnotation -G AS_StandardAnnotation 
    """


rule extract_exomes_gvcf:
    input:
        gvcf = rules.HaplotypeCaller.output.gvcf,
        tbi = rules.HaplotypeCaller.output.gvcf_tbi
        # gvcf = pj(GVCF, "reblock/{region}/{sample}.{region}.wg.vcf.gz"),
        # tbi = pj(GVCF, "reblock/{region}/{sample}.{region}.wg.vcf.gz.tbi"),
    output:
        gvcf_exome = ensure( pj(GVCF, "exome_extract/{region}/{sample}.{region}.wg.vcf.gz"), non_empty=True),
        tbi = ensure( pj(GVCF, "exome_extract/{region}/{sample}.{region}.wg.vcf.gz.tbi"), non_empty=True),
    conda: CONDA_VCF
    params: java_options=DEFAULT_JAVA_OPTIONS,
            interval = lambda wildcards: region_to_file(region = wildcards.region, extension="interval_list"),
            padding = 500,
    resources: n= "1.0",
               mem_mb= 1500,
    run:
        if "wgs" in SAMPLEINFO[wildcards.sample]["sample_type"] :
            shell(
                """
                    gatk --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" SelectVariants \
                    -V {input.gvcf} -O {output.gvcf_exome} \
                    -L {params.interval} -ip {params.padding} --seconds-between-progress-updates 120    
                """),
        else:
            shell(
                """
                    cp {input.gvcf} {output.gvcf_exome}
                    cp {input.tbi} {output.tbi}
                """)


