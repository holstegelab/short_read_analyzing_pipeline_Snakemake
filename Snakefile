import pandas as pd
import read_stats
import itertools
import os


ruleorder: adapter_removal > fastq_bz2togz

wildcard_constraints:
    sample=r"[\w\d_\-@]+",
# readgroup="[\w\d_\-@]+"

from common import *

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

# module Genotype:
#     snakefile: 'Genotype.smk'
#     config: config
# module VQSR:
#     snakefile: 'VQSR.smk'
#     config: config
module Stat:
    snakefile: 'Stat.smk'
    config: config

use rule * from Stat

module PCA:
    snakefile: 'PCA.smk'
    config: config

# module SV_delly:
#     snakefile: 'SV_delly.smk'
#     config: config
# module CNV_with_cnvkit_Module:
#     snakefile: 'CNV_with_cnvkit_Module.smk'
#     config: config

module Combine_gVCF:
    snakefile: 'Combine_gVCF.smk'
    config: config


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
END_RULE = []
CLEAN_RULE = []
chrM = config.get("chrM","No")
if chrM == "Yes":
    use rule * from chrM_analysis
    chrM_rule = rules.chrM_analysis_all.input
else:
    chrM_rule = []
gVCF_combine_method = config.get("Combine_gVCF_method","GLnexus")

gvcf_caller = config.get("caller","Deepvariant")
glnexus_filtration = config.get("glnexus_filtration","custom")

rule_all_combine = []
VQSR_rule = []
VQSR = config.get("VQSR","NO")
end_point = config.get("END_POINT","gVCF")

print(end_point, gvcf_caller)
if end_point == "gVCF":
    if gvcf_caller == "BOTH":
        use rule * from gVCF

        use rule * from Deepvariant

        use rule * from Aligner

        use rule * from Kraken
        
        use rule * from Stat

        rule finished_sample:
            """Finish processing a sample. 

            This rule will drop the reservation of space on active storage.
            """
            input:
                pj(CRAM,"{sample}.mapped_hg38.cram.copied"),
                pj(GVCF,"{sample}.done"),
                pj(DEEPVARIANT,"{sample}.done"),
                pj(STAT,"{sample}.done"),
                pj(KRAKEN,"{sample}.bracken_report.tsv")

            output:
                pj(SOURCEDIR,"{sample}.finished")
            resources:
                active_use_remove=Aligner.calculate_active_use,
                mem_mb=50,
                n="1"
            shell: """
                touch {output}
                """



        print("You will run following steps: Aligning with dragen and gVCF calling with HaplotypeCaller and Deepvariant (both=default). \n"
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
                pj(GVCF,"{sample}.done"),
                pj(STAT,"{sample}.done"),
                pj(KRAKEN,"{sample}.bracken_report.tsv")
            output:
                os.path.join(SOURCEDIR, "{sample}.finished")
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
        use rule * from Aligner
        
        use rule * from Deepvariant

        use rule * from Kraken

        use rule * from Stat

        rule finished_sample:
            """Finish processing a sample. 

            This rule will drop the reservation of space on active storage.
            """
            input:
                pj(CRAM,"{sample}.mapped_hg38.cram.copied"),
                pj(DEEPVARIANT,"{sample}.done"),
                pj(STAT,"{sample}.done"),
                pj(KRAKEN,"{sample}.bracken_report.tsv")
            output:
                pj(SOURCEDIR,"{sample}.finished")
            resources:
                active_use_remove=Aligner.calculate_active_use,
                mem_mb=50,
                n="1"
            shell: """
                touch {output}
            """

        print("You will run following steps: Aligning with dragen and gVCF calling with Deepvariant. \n"
              "To change gVCF caller to HaplotypeCaller pass '--config caller=HaplotypeCaller'")

    END_RULE = [expand("{source}/{sample}.finished",sample=sample_names,source=SOURCEDIR), rules.Stat_all.input]
    CLEAN_RULE = [expand(pj(STAT,"{sample}.stats.tar.gz"),sample=sample_names)]

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
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with GenomicDBimport and Genotyping with GATK Genotype \n"
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant' \n"
                  "* To change combining method to GATK-s Combinbegvcf pass '--config Combine_gVCF_method=COMBINE_GVCF' \n"
                  "* To change jointgenotyping method to GLnexus pass --config Combine_gVCF_method=GLnexus")
        elif gVCF_combine_method == "GLnexus":
            use rule * from GLnexus

            rule_all_combine = rules.GLnexus_all.input
            END_RULE = rules.GLnexus_all.input
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with GLnexus (default) \n"
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant' \n"
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
        print("You will run following steps: Aligning with dragen, gVCF calling with Deepvariant, merging gVCFs with GLnexus (default) \n"
              "* To change gVCF caller to HaplotypeCaller pass '--config caller=HaplotypeCaller' \n "
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
elif end_point == "Combine":
    use rule * from Aligner

    if gVCF_combine_method == "DBIMPORT":
        use rule * from gVCF

        use rule * from DBImport

        END_RULE = rules.DBImport_all.input
        rule_all_combine = rules.DBImport_all.input
    elif gVCF_combine_method == "COMBINE_GVCF":
        use rule * from gVCF

        END_RULE = rules.Combine_gVCF_all.input

        use rule * from Combine_gVCF

        rule_all_combine = rules.Combine_gVCF_all.input
    elif gVCF_combine_method == "GLnexus":
        use rule * from Aligner

        if gvcf_caller == "HaplotypeCaller":
            use rule * from gVCF

            use rule * from GLnexus

            END_RULE = rules.GLnexus_all.input
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs and Genotyping with GLnexus \n"
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant' \n "
                  "* To To change combining method to GATK-s GenomicDBimport pass '--config Combine_gVCF_method=DBIMPORT'")
        elif gvcf_caller == "Deepvariant":
            use rule * from Deepvariant

            use rule * from GLnexus

            rule_all_combine = rules.GLnexus_all.input
            END_RULE = rules.GLnexus_all.input
            print("You will run following steps: Aligning with dragen, gVCF calling with Deepvariant, merging gVCFs with GLnexus (default) \n"
                  "* To change gVCF caller to HaplotypeCaller pass '--config caller=HaplotypeCaller' \n"
                  "* To change combining method to GATK-s GenomicDBimport pass '--config Combine_gVCF_method=DBIMPORT' \n")
        elif gvcf_caller == "BOTH":
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
elif end_point == "VQSR" or end_point == "VCF":
    use rule * from Aligner

    if gVCF_combine_method == "DBIMPORT":
        use rule * from gVCF

        use rule * from Genotype

        END_RULE = rules.Genotype_all.input
        if VQSR == "RUN_VQSR":
            use rule * from VQSR

            VQSR_rule = rules.VQSR_all.input,
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with GenomicDBimport and Genotyping with GATK Genotype. Additionallly VQSR will be done. \n"
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant' \n"
                  "* To change combining method to GATK-s Combinbegvcf pass '--config Combine_gVCF_method=COMBINE_GVCF' \n"
                  "* To change jointgenotyping method to GLnexus pass --config Combine_gVCF_method=GLnexus \n"
                  "* To remove VQSR step pass '--config VQSR=NO")
        elif VQSR == "NO" or VQSR == "NO_VQSR" or VQSR == "NO_RUN":
            VQSR_rule = []
            print("You will run following steps: Aligning with dragen, gVCF calling with HaplotypeCaller (default), merging gVCFs with GenomicDBimport and Genotyping with GATK Genotype \n"
                  "* To change gVCF caller to deepvariant pass '--config caller=Deepvariant' \n "
                  "* To change combining method to GATK-s Combinbegvcf pass '--config Combine_gVCF_method=COMBINE_GVCF' \n"
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

            END_RULE = rules.GLnexus_all.input
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

rule pipeline:
    input:
        END_RULE,
        chrM_rule,
        #SV_rule,
        #CNV_rule,
        rules.Encrypt_all.input,
    output:
        done=touch('pipeline.done')

rule all:
    input:
        CLEAN_RULE,
        rules.pipeline.output.done
    default_target: True

sample_names = SAMPLEINFO.keys()
sample_pattern = "|".join(sample_names)
onsuccess: shell(# "rm -f zslurm-*"
                 # "rm -rf logs"
                 "rm -rf tmp")
onerror:
            shell("""
            sample_pattern="{sample_pattern}"
            grep 'Error in rule' zslurm-logs/* | awk '{{print $1 "\t" $4}}' | awk -F"[/:]" '{{print$1 "\t" $2}}' | awk '{{print$1 "\t" $3}}'>> error_rules.txt
            grep -A 2 'Error in rule' zslurm-* | grep 'input' | awk -F[,] '{{print$1}}' | grep -E -o '{sample_pattern}' >> error_samples.txt
            paste error_rules.txt error_samples.txt > error.log
            """)


