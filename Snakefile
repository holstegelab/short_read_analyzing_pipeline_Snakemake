import pandas as pd
import read_stats
import itertools
import os
configfile: srcdir("Snakefile.cluster.json")
configfile: srcdir("Snakefile.paths.yaml")

wildcard_constraints:
    sample="[\w\d_\-@]+",
    # readgroup="[\w\d_\-@]+"

from common import *
SAMPLE_FILES, SAMPLEFILE_TO_SAMPLES, SAMPLEINFO, SAMPLE_TO_BATCH, SAMPLEFILE_TO_BATCHES = load_samplefiles('.',config)
sample_names = SAMPLEINFO.keys()

# extract all sample names from SAMPLEINFO dict to use it rule all
sample_names = SAMPLEINFO.keys()

module Aligner:
    snakefile: 'Aligner.smk'
    config: config
module gVCF:
    snakefile: 'gVCF.smk'
    config: config
module Kraken:
    snakefile: 'Kraken.smk'
    config: config
module DBImport:
    snakefile: 'DBImport.smk'
    config: config
module Genotype:
    snakefile: 'Genotype.smk'
    config: config
module VQSR:
    snakefile: 'VQSR.smk'
    config: config
module Stat:
    snakefile: 'Stat.smk'
    config: config
module SV_delly:
    snakefile: 'SV_delly.smk'
    config: config
module CNV_with_cnvkit_Module:
    snakefile: 'CNV_with_cnvkit_Module.smk'
    config: config

module Combine_gVCF:
    snakefile: 'Combine_gVCF.smk'
    config: config
use rule * from Stat

module chrM_analysis:
    snakefile: 'chrM_analysis.smk'
    config: config
use rule * from chrM_analysis

module Encrypt:
    snakefile: 'Encrypt.smk'
    config: config
use rule * from Encrypt

module Deepvariant:
    snakefile: 'Deepvariant.smk'
    config: config

module GLnexus:
    snakefile: 'GLnexus.smk'
    config: config

module Reference_preparation:
    snakefile: "Reference_preparation.smk"
    config: config



#SV = config.get("SV", "RUN_SV")
#if SV == "RUN_SV":
#    use rule * from SV_delly
#    SV_rule = rules.SV_delly_all.input
#else:
#    SV_rule = []

#CNV = config.get("CNV", "RUN_CNV")
#if CNV == "RUN_CNV":
#    use rule * from CNV_with_cnvkit_Module
#    CNV_rule = rules.CNV_with_cnvkit_Module_all.input
#else:
#    CNV_rule = []

gVCF_combine_method = config.get("Combine_gVCF_method", "GLnexus")


gvcf_caller = config.get("caller", "both")
glnexus_filtration = config.get("glnexus_filtration", "default")


rule_all_combine = []
VQSR_rule = []
VQSR = config.get("VQSR","NO")
end_point = config.get("END_POINT", "gVCF")

if end_point == "gVCF":
    if gvcf_caller == "both":
        use rule * from gVCF
        use rule * from Deepvariant
        use rule * from Aligner
        use rule * from Kraken

        rule finished_sample:
            """Finish processing a sample. 
            
            This rule will drop the reservation of space on active storage.
            """
            input:
                pj(CRAM,"{sample}.mapped_hg38.cram.copied"),
                pj(GVCF, "{sample}.done"),
                pj(DEEPVARIANT, "{sample}.done"),
                pj(STAT, "{sample}.done"),
                pj(KRAKEN, "{sample}.bracken_report.tsv")

            output:
                pj(SOURCEDIR, "{sample}.finished")
            resources:
                active_use_remove=Aligner.calculate_active_use,
                mem_mb=50,
                n="1"
            shell: """
                touch {output}
                """
        print("You will run following steps: Aligning with dragen and gVCF calling with HaplotypeCaller and Deepvariant (both=default). "
              "To change gVCF caller selection pass '--config caller=Deepvariant' or '--config caller=HaplotypeCaller'")
    elif gvcf_caller == "HaplotypeCaller":
        use rule * from gVCF
        use rule * from Aligner
        use rule * from Kraken
        
        rule finished_sample:
            """Finish processing a sample. 
            
            This rule will drop the reservation of space on active storage.
            """
            input:
                pj(CRAM,"{sample}.mapped_hg38.cram.copied"),
                pj(GVCF, "{sample}.done"),
                pj(STAT, "{sample}.done"),
                pj(KRAKEN, "{sample}.bracken_report.tsv")
            output:
                os.path.join(config['SOURCEDIR'], "{sample}.finished")
            resources:
                active_use_remove=Aligner.calculate_active_use,
                mem_mb=50,
                n="1"
            shell: """
                touch {output}
                """

        print("You will run following steps: Aligning with dragen and gVCF calling with HaplotypeCaller (default). "
              "To change gVCF caller to deepvariant pass '--config caller=Deepvariant'")
    elif gvcf_caller == "Deepvariant":
        use rule * from Deepvariant
        use rule * from Aligner
        use rule * from Kraken
        
        rule finished_sample:
            """Finish processing a sample. 
            
            This rule will drop the reservation of space on active storage.
            """
            input:
                pj(CRAM,"{sample}.mapped_hg38.cram.copied"),
                pj(DEEPVARIANT, "{sample}.done"),
                pj(STAT, "{sample}.done"),
                pj(KRAKEN, "{sample}.bracken_report.tsv")
            output:
                pj(SOURCEDIR, "{sample}.finished")
            resources:
                active_use_remove=Aligner.calculate_active_use,
                mem_mb=50,
                n="1"
            shell: """
                touch {output}
                """
        print("You will run following steps: Aligning with dragen and gVCF calling with Deepvariant. "
              "To change gVCF caller to HaplotypeCaller pass '--config caller=HaplotypeCaller'")
    END_RULE = [expand("{source}/{sample}.finished",sample=sample_names, source = SOURCEDIR), rules.Stat_all.input]
elif end_point == 'PrepareRef':
    use rule * from Reference_preparation
    END_RULE = rules.Reference_preparation_all.input
elif end_point == "Align" or end_point == "Aligner":
    use rule * from Aligner
    END_RULE = rules.Aligner_all.input
    print("You will run following steps: Aligning with dragen")
elif end_point == "Genotype" or end_point == "Genotyper":
    use rule * from Aligner
    if gvcf_caller == "HaplotypeCaller":
        use rule * from gVCF
        if gVCF_combine_method == "DBIMPORT" or gVCF_combine_method == "COMBINE_GVCF":
            use rule * from Genotype
            END_RULE = rules.gVCF_all.input
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with GenomicDBimport and Genotyping with GATK Genotype "
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant'"
                  "* To change combining method to GATK-s Combinbegvcf pass '--config Combine_gVCF_method=COMBINE_GVCF'"
                  "* To change jointgenotyping method to GLnexus pass --config Combine_gVCF_method=GLnexus")
        elif gVCF_combine_method == "GLnexus":
            use rule * from GLnexus
            rule_all_combine = rules.GLnexus_all.input
            END_RULE = rules.GLnexus_all.input
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with GLnexus (default)"
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant'"
                  "* To change combining method to GATK-s GenomicDBimport pass '--config Combine_gVCF_method=DBIMPORT'")
        else:
            raise ValueError(
                "invalid option provided to 'Combine_gVCF_method'; please choose either 'GLnexus'(default), 'COMBINE_GVCF' or 'DBIMPORT'."
            )
    elif gvcf_caller == "Deepvariant":
        use rule * from Deepvariant
        use rule * from GLnexus
        rule_all_combine = rules.GLnexus_all.input
        END_RULE = rules.GLnexus_all.input
        print("You will run following steps: Aligning with dragen, gVCF calling with Deepvariant, merging gVCFs with GLnexus (default)"
              "* To change gVCF caller to HaplotypeCaller pass '--config caller=HaplotypeCaller'"
              "* To change combining method to GATK-s GenomicDBimport pass '--config Combine_gVCF_method=DBIMPORT'")
    elif gvcf_caller == "BOTH":
        use rule * from Genotype
        use rule * from Deepvariant
        use rule * from GLnexus
        rule_all_combine = rules.GLnexus_all.input
        END_RULE = rules.GLnexus_all.input
        print("You will run following steps: Aligning with dragen, gVCF calling with Haplotypecaller, merging gVCFs with GLnexus (default) and separete gVCF calling with Deepvariant")
    else:
        raise ValueError(
            "invalid option provided to 'caller'; please choose either 'HaplotypeCaller'(default) or 'Deepvariant'."
        )

elif end_point == "VQSR" or end_point == "VCF" or end_point == "Combine":
    use rule * from Aligner
    if gVCF_combine_method == "DBIMPORT":
        use rule * from gVCF
        use rule * from Genotype
        END_RULE = rules.Genotype_all.input
        if VQSR == "RUN_VQSR":
            use rule * from VQSR
            VQSR_rule = rules.VQSR_all.input,
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with GenomicDBimport and Genotyping with GATK Genotype. Additionallly VQSR will be done. "
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant'"
                  "* To change combining method to GATK-s Combinbegvcf pass '--config Combine_gVCF_method=COMBINE_GVCF'"
                  "* To change jointgenotyping method to GLnexus pass --config Combine_gVCF_method=GLnexus"
                  "* To remove VQSR step pass '--config VQSR=NO")
        elif VQSR == "NO" or VQSR == "NO_VQSR" or VQSR == "NO_RUN":
            VQSR_rule = []
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with GenomicDBimport and Genotyping with GATK Genotype "
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant'"
                  "* To change combining method to GATK-s Combinbegvcf pass '--config Combine_gVCF_method=COMBINE_GVCF'"
                  "* To change jointgenotyping method to GLnexus pass --config Combine_gVCF_method=GLnexus")
        else:
            raise ValueError(
                "invalid option provided to 'VQSR'; please choose either 'RUN_VQSR' or 'NO_VQSR(default)'."
            )
        use rule * from DBImport
        rule_all_combine = rules.DBImport_all.input
    elif gVCF_combine_method == "COMBINE_GVCF":
        use rule * from gVCF
        use rule * from Genotype
        END_RULE = rules.Genotype_all.input
        if VQSR == "RUN_VQSR":
            use rule * from VQSR
            VQSR_rule = rules.VQSR_all.input,
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with Combinegvcf and Genotyping with GATK Genotype.  Additionallly VQSR will be done. "
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant'"
                  "* To change combining method to GATK-s Combinbegvcf pass '--config Combine_gVCF_method=DBIMPORT'"
                  "* To change jointgenotyping method to GLnexus pass --config Combine_gVCF_method=GLnexus"
                  "* To remove VQSR step pass '--config VQSR=NO")
        elif VQSR == "NO" or VQSR == "NO_VQSR" or VQSR == "NO_RUN":
            VQSR_rule = []
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with Combinegvcf and Genotyping with GATK Genotype "
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant'"
                  "* To change combining method to GATK-s Combinbegvcf pass '--config Combine_gVCF_method=DBIMPORT'"
                  "* To change jointgenotyping method to GLnexus pass --config Combine_gVCF_method=GLnexus")
        else:
            raise ValueError(
                "invalid option provided to 'VQSR'; please choose either 'RUN_VQSR' or 'NO_VQSR(default)'."
            )
        use rule * from Combine_gVCF
        rule_all_combine = rules.Combine_gVCF_all.input
    elif gVCF_combine_method == "GLnexus":
        use rule * from Aligner
        if gvcf_caller == "HaplotypeCaller":
            use rule * from gVCF
            use rule * from GLnexus
            END_RULE = END_RULE = rules.GLnexus_all.input
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs and Genotyping with GLnexus "
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant'"
                  "* To To change combining method to GATK-s GenomicDBimport pass '--config Combine_gVCF_method=DBIMPORT'")
        elif gvcf_caller == "Deepvariant":
            use rule * from Deepvariant
            use rule * from GLnexus
            rule_all_combine = rules.GLnexus_all.input
            END_RULE = rules.GLnexus_all.input
            print("You will run following steps: Aligning with dragen, gVCF calling with Deepvariant, merging gVCFs with GLnexus (default)"
                  "* To change gVCF caller to HaplotypeCaller pass '--config caller=HaplotypeCaller'"
                  "* To change combining method to GATK-s GenomicDBimport pass '--config Combine_gVCF_method=DBIMPORT'")
        elif gvcf_caller == "BOTH":
            use rule * from Genotype
            use rule * from Deepvariant
            use rule * from GLnexus
            rule_all_combine = rules.GLnexus_all.input
            END_RULE = rules.GLnexus_all.input
            print("You will run following steps: Aligning with dragen, gVCF calling with Haplotypecaller, merging gVCFs with GLnexus (default) and separete gVCF calling with Deepvariant")
        else:
            raise ValueError(
                "invalid option provided to 'caller'; please choose either 'HaplotypeCaller'(default) or 'Deepvariant'."
            )
    else:
        raise ValueError(
            "invalid option provided to 'Combine_gVCF_method'; please choose either 'GLnexus'(default), 'COMBINE_GVCF' or 'DBIMPORT'."
        )
else:
    raise ValueError(
        "Invalid option provided to 'END_POINT'; please choose either 'gVCF(default)', 'Align', 'Genotype' or 'Combine'."
    )



rule all:
    input:
        END_RULE,
        rules.chrM_analysis_all.input,
        #SV_rule,
        #CNV_rule,
        rules.Encrypt_all.input,
    default_target: True





