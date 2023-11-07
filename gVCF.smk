from common import *

wildcard_constraints:
    sample="[\w\d_\-@]+",
    # readgroup="[\w\d_\-@]+"




module Aligner:
    snakefile: 'Aligner.smk'
    config: config
use rule * from Aligner

module Stat:
    snakefile: 'Stat.smk'
    config: config
use rule verifybamid from Stat

module Tools:
    snakefile: 'Tools.smk'
    config: config
use rule * from Tools



rule gVCF_all:
    input:
        expand("{gvcf}/{sample}.done",sample=sample_names, gvcf = GVCF),
        rules.Aligner_all.input
    default_target: True

def get_gvcf_files(wildcards):#{{{
    sample = wildcards['sample']
    regions = level1_regions if 'wgs' in SAMPLEINFO[sample]['sample_type'] else level0_regions    
    return [pj(GVCF, region, f'{sample}.{region}.wg.vcf.gz') for region in regions]
#}}}

rule gvcf_sample_done:
    input:
        get_gvcf_files
    output:
        temp(touch(pj(GVCF, "{sample}.done")))
    resources:
        n="1.0",
        mem_mb=50



rule ComposeSTRTableFile:
    input:
        ref_fasta = "{path}.fa"
    output:
        str_table = "{path}.str.zip"
    params:
        java_options=DEFAULT_JAVA_OPTIONS
    resources: 
        n = "1.0",
        mem_mb = 3000
    conda: CONDA_VCF   
    shell:
        """gatk ComposeSTRTableFile --java-options "-Xmx{resources.mem_mb}m {params.java_options}" -R {input.ref_fasta} -O {output.str_table}"""


rule CalibrateDragstrModel:
    """CalibrateDragstrModel. Estimates the parameters of the dragstr model from a set of aligned reads.
    Provides better STR variant calling.
    """
    input:
        bam = rules.markdup.output.mdbams,
        bai= rules.markdup.output.mdbams_bai,
        validated_sex = rules.get_validated_sex.output.yaml       
    output:
        dragstr_model = pj(BAM, "{sample}-dragstr.txt")
    priority: 26
    params:
        ref=get_ref_by_validated_sex,
        str_ref = get_strref_by_validated_sex,
        java_options=DEFAULT_JAVA_OPTIONS
    resources: 
        n = "1.0",
        mem_mb = lambda wildcards, attempt: attempt * 1750
    log: pj(LOG, "{sample}_calibratedragstr.log")
    benchmark: pj(BENCH, "{sample}_calibrate_dragstr.txt")
    conda: CONDA_VCF
    shell:
        """{gatk} CalibrateDragstrModel --java-options \
                    "-Xmx{resources.mem_mb}m {params.java_options}"  -R {params.ref} -I {input.bam} \
                    -O {output} -str {params.str_ref} 2>{log}"""

def read_contam_w(wildcards):#{{{
    """Read contamination from verifybamid output file."""
    filename = pj(STAT, 'contam', wildcards['sample'] +  '.verifybamid.pca2.selfSM')
    with open(filename,'r', encoding='utf-8') as f:
        c = csv.reader(f, delimiter='\t')
        c = iter(c)
        headers = c.__next__()
        data = c.__next__()
        freemix = float(data[6])

    if freemix < 0.01: #only remove contamination if it is more than 1%
        freemix = 0.0        
    return freemix
#}}}

def get_mem_mb_HaplotypeCaller(wildcards, attempt):#{{{
    """Get memory for HaplotypeCaller."""
    res = 2400 if 'wgs' in SAMPLEINFO[wildcards['sample']]['sample_type'] else 1400

    #HaplotypeCaller has exponential memory scaling on some regions. Average is very  low, but on some regions it can scale to 75 GB...
    return (res * (3 ** (attempt - 1)))  #aggressively reserve more memory
#}}}

def region_to_interval_file(wildcards):#{{{
    """Converts a region to a interval file location (see common.py and Tools.smk)"""
    sample = wildcards['sample']
    region = wildcards['region']
    return region_to_file(region, wgs='wgs' in SAMPLEINFO[sample]['sample_type'], extension='interval_list')
#}}}

rule HaplotypeCaller:
    """HaplotypeCaller. Call SNPs and indels for each sample."""
    input:
        # check what bam file we need to use (with or without additional cleanup)
        bams = rules.markdup.output.mdbams,
        bai = rules.markdup.output.mdbams_bai,
        model = rules.CalibrateDragstrModel.output.dragstr_model,
        contam= rules.verifybamid.output.VBID_stat,
        interval= region_to_interval_file,
        validated_sex = rules.get_validated_sex.output.yaml
    output:
        orig_gvcf= ensure(temp(pj(GVCF, "{region}/{sample}.{region}.g.vcf.gz")), non_empty=True),
        orig_gvcf_tbi = ensure(temp(pj(GVCF, "{region}/{sample}.{region}.g.vcf.gz.tbi")), non_empty=True),
        genotyped_vcf = temp(pj(GVCF, "{region}/{sample}.{region}.vcf.gz")),
        genotyped_vcf_tbi = temp(pj(GVCF, "{region}/{sample}.{region}.vcf.gz.tbi")),
        tmp_vcf = temp(pj(GVCF, "{region}/{sample}.{region}.tmp.vcf.gz")),
        vcf = temp(pj(GVCF, "{region}/{sample}.{region}.w.vcf.gz")),
        vcf_tbi = temp(pj(GVCF, "{region}/{sample}.{region}.w.vcf.gz.tbi")),
        wstats = pj(STAT, "whatshap_phasing/{sample}.{region}.stats"),
        mwstats = pj(STAT, "whatshap_phasing/{sample}.{region}.merge_stats"),
        tmp_gvcf= temp(pj(GVCF, "{region}/{sample}.{region}.wg.vcf")), 
        gvcf= pj(GVCF, "{region}/{sample}.{region}.wg.vcf.gz"), 
        gvcf_tbi = pj(GVCF, "{region}/{sample}.{region}.wg.vcf.gz.tbi"), 
    log:
        HaplotypeCaller=pj(LOG, "{sample}_{region}_haplotypecaller.log")
    benchmark:
        pj(BENCH, "{sample}_{region}_haplotypecaller.txt")
    conda: CONDA_VCF
    resources: 
               n="0.9", #average 1.3 cores
               mem_mb = get_mem_mb_HaplotypeCaller,
               tmpdir = tmpdir_alternative
    params:
        dbsnp = DBSNP,
        padding=0,  #padding is included in the interval file. Keep at 0 to prevent overlapping intervals
        contam_frac = read_contam_w, # get contamination fraction per sample
        ploidy = lambda wildcards: 1 if wildcards['region'].endswith('H') else 2,
        java_options=DEFAULT_JAVA_OPTIONS,
        ref=get_ref_by_validated_sex,
        dragen_mode = lambda wildcards: '--dragen-mode true' if not 'H' in wildcards['region'] else '',
        merge_script=srcdir(MERGEPHASE)
    priority: 28
    #Overview of available annotations (GATK version 4.4)
    #AS_StandardAnnotation: AS_BaseQualityRankSumTest, AS_FisherStrand, AS_InbreedingCoeff, AS_MappingQualityRankSumTest, AS_QualByDepth, AS_RMSMappingQuality, AS_ReadPosRankSumTest, AS_StrandOddsRatio (all allele specific annotations)
    #StandardHCAnnotation: DepthPerSampleHC
    #StandardANnotation: Coverage, ChromosomeCounts, BaseQualityRankSumTest, MappingQualityRankSumTest, FisherStrand, QualByDepth, InbreedingCoeff, DepthPerAlleleBySample, ReadPosRankSumTest, ExcessHet
    #                      StrandOddsRatio, RMSMappingQuality, 
    #not included in standard annotation: FragmentDepthPerAlleleBySample,StrandBiasBySample, AssemblyComplexity,  FragmentLength, OrientationBiasReadCounts, TandemRepeat, UniqueAltReadCount
    
    #findings:
    #FragmentDepthPerAlleleBySample, OrientationBiasReadCounts, TandemRepeat --> do not output anything in haplotypecaller, omitted

    #regarding linked-de-bruijn-graph, see comments in https://gatk.broadinstitute.org/hc/en-us/community/posts/360077647812-Why-do-a-clear-expected-variant-not-show-up-in-the-Mutect2-vcf-file
    #also: https://gatk.broadinstitute.org/hc/en-us/articles/360043491652-When-HaplotypeCaller-and-Mutect2-do-not-call-an-expected-variant
    #still, seems quite experimental, so not using it for now
    
    #related: --recover-all-dangling-branches
    #to consider: --population-callset in combination with gnomad,  --use-posteriors-to-calculate-qual

    #max-mnp-distance: not compatible with multiple sample gVCF calling, omitted.

    #PHASING
    #GATK: physical phasing depends on de bruijn graph, which is not always optimal. In practice this means
    #that the phasing is only short distance (~50bp). Whathap uses the full dna fragment, usually 300bp, to phase.
    #Here we use whatshap to phase the gatk called variants, and then merge the whatshap phased variants with the gatk phased variants.
    #This way we get the best of both worlds: the gatk phasing is used for short distance phasing also of variants with low read
    #support (which might become full variants during multi-sample genotyping), and the whatshap phasing is used for long distance phasing.
    
    #WORKAROUND: Due to bug in whatshap stats for empty contigs, we use || true to prevent pipeline from crashing when no heterozygous variants 
    #are found.
    shell: """ 
        {gatk} --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" HaplotypeCaller     \
            -R {params.ref} -L {input.interval} -ip {params.padding} -D {params.dbsnp} -ERC GVCF --contamination {params.contam_frac} \
            --ploidy {params.ploidy} -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
            --annotate-with-num-discovered-alleles --adaptive-pruning \
            -A StrandBiasBySample -A AssemblyComplexity -A FragmentLength \
            -I {input.bams} -O {output.orig_gvcf}  --native-pair-hmm-threads 2  --create-output-variant-index true\
            {params.dragen_mode} --dragstr-params-path {input.model} 2> {log.HaplotypeCaller}

        {gatk} --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" GenotypeGVCFs \
                -R {params.ref} -V {output.orig_gvcf} -O {output.genotyped_vcf} 

        mkdir -p `dirname {output.wstats}`
        if [ {params.ploidy} -eq 2 ]
        then 
            whatshap unphase  {output.genotyped_vcf} >  {output.tmp_vcf}
            whatshap phase  --ignore-read-groups --reference {params.ref} {output.tmp_vcf} {input.bams} -o {output.vcf}
            bcftools index --tbi {output.vcf}
            whatshap stats {output.genotyped_vcf} > {output.wstats} || true
            whatshap stats {output.vcf} >> {output.wstats} || true
            python {params.merge_script} {output.orig_gvcf} {output.vcf} {output.tmp_gvcf} {output.mwstats}
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
        bcftools index --tbi {output.gvcf}          
        """

rule reblock_gvcf:
    input:
        gvcf = rules.HaplotypeCaller.output.orig_gvcf,
        idx = rules.HaplotypeCaller.output.orig_gvcf_tbi,
        validated_sex = rules.get_validated_sex.output.yaml
    output: gvcf_reblock = ensure(pj(GVCF, "reblock/{region}/{sample}.{region}.g.vcf.gz"), non_empty=True),
            tbi = ensure(pj(GVCF, "reblock/{region}/{sample}.{region}.g.vcf.gz.tbi"), non_empty=True)
    log: Reblock=pj(LOG, "{sample}_{region}_reblock.log")
    benchmark:
        pj(BENCH, "{sample}_{region}_reblock.txt")
    conda: CONDA_VCF
    priority: 29
    params:
        dbsnp=DBSNP,
        java_options=DEFAULT_JAVA_OPTIONS,
        ref=get_ref_by_validated_sex
    resources: 
        n="2.0",
        mem_mb = lambda wildcards, attempt: attempt * 2500
    shell:
        """{gatk} --java-options "-Xmx{resources.mem_mb}M  {params.java_options}" ReblockGVCF   --keep-all-alts --create-output-variant-index true -D {params.dbsnp} -R {params.ref} -V {input.gvcf} -O {output.gvcf_reblock} -GQB 3 -GQB 5 -GQB 8 -GQB 10 -GQB 15 -GQB 20 -GQB 30 -GQB 50 -GQB 70 -GQB 100 -G StandardAnnotation -G AS_StandardAnnotation 2> {log}"""

