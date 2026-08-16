import pandas as pd
import read_stats
import itertools
import os



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
chrM_flag = config.get("chrM","Yes")
if chrM_flag == "Yes":
    use rule * from chrM_analysis
    chrM_rule = rules.chrM_analysis_all.input
else:
    chrM_rule = []
gVCF_combine_method = config.get("Combine_gVCF_method","GLnexus")

gvcf_caller = config.get("caller","Deepvariant")
glnexus_filtration = config.get("glnexus_filtration","custom")


def finished_sample_inputs(wildcards, require_gvcf=False, require_deepvariant=False):
    files = [
        pj(CRAM, f"{wildcards.sample}.mapped_hg38.cram.copied"),
        pj(STAT, f"{wildcards.sample}.stats.tar.gz"),
        pj(KRAKEN, f"{wildcards.sample}.bracken_report.tsv"),
    ]
    if require_gvcf:
        files.append(pj(GVCF, f"{wildcards.sample}.done"))
    if require_deepvariant:
        files.append(pj(DEEPVARIANT, f"{wildcards.sample}.done"))
    if chrM_flag == "Yes":
        files.append(pj(chrM, f"{wildcards.sample}.done"))
    if SAMPLEINFO[wildcards.sample].get('from_external'):
        files.append(pj(
            SOURCEDIR, f"{wildcards.sample}.materialized_consumed"
        ))
    return files

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
                lambda wc: finished_sample_inputs(wc, require_gvcf=True, require_deepvariant=True)

            output:
                pj(SOURCEDIR,"{sample}.finished")
            resources:
                time = get_time('finished_sample'),
                active_use_remove=active_release_finished,
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
                lambda wc: finished_sample_inputs(wc, require_gvcf=True)
            output:
                os.path.join(SOURCEDIR, "{sample}.finished")
            resources:
                time = get_time('finished_sample'),
                active_use_remove=active_release_finished,
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
                lambda wc: finished_sample_inputs(wc, require_deepvariant=True)
            output:
                pj(SOURCEDIR,"{sample}.finished")
            resources:
                time = get_time('finished_sample'),
                active_use_remove=active_release_finished,
                mem_mb=50,
                n="1"
            shell: """
                touch {output}
            """

        print("You will run following steps: Aligning with dragen and gVCF calling with Deepvariant. \n"
              "To change gVCF caller to HaplotypeCaller pass '--config caller=HaplotypeCaller'")

    END_RULE = [
        expand("{source}/{sample}.finished", sample=sample_names, source=SOURCEDIR),
        rules.Stat_all.input,
        rules.kraken_tar_all.output.done,
        rules.badmap_tar_all.output.done,
        rules.stats_to_dcache_all.output.done,
        rules.excluded_to_dcache_all.output.done,
    ]
    if gvcf_caller in ("Deepvariant", "BOTH"):
        END_RULE.append(rules.deepvariant_tar_wgs_all.output.done)
        END_RULE.append(rules.deepvariant_tar_wes_all.output.done)
    if chrM_flag == "Yes":
        END_RULE.append(rules.chrM_tar_all.output.done)
    CLEAN_RULE = [
        expand(pj(STAT,"{sample}.stats.tar.gz"), sample=sample_names)
    ]

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
        # For the gVCF endpoint END_RULE already contains the durable chrM
        # upload marker (and each finished_sample contains chrM/{sample}.done).
        # Pulling in chrM_analysis_all here as well reintroduced Aligner_all's
        # temporary mapped CRAMs and loose chrM gVCFs after they had been
        # uploaded and garbage-collected, causing a resume to rebuild samples.
        [] if end_point == "gVCF" else chrM_rule,
        #SV_rule,
        #CNV_rule,
        # rules.Encrypt_all.input,
    output:
        done=touch('pipeline.done')

rule all:
    input:
        CLEAN_RULE,
        rules.pipeline.output.done
    default_target: True

sample_names = SAMPLEINFO.keys()
sample_pattern = "|".join(sample_names)

onstart:
    # Reclaim temp files that Snakemake's temp() GC leaves behind across a RESTART.
    # A consuming job that already ran in a prior invocation is skipped on rerun, so
    # its temp outputs are never collected and pile up on active storage. Two gates,
    # each removing only files whose consumers are provably done, and NEVER touching
    # the .started/.finished/.copied markers themselves. Wrapped so it can never abort.
    import glob, os, shutil
    def _rm(path):
        try:
            if os.path.isdir(path):
                shutil.rmtree(path, ignore_errors=True)
            else:
                os.remove(path)
        except OSError:
            pass
    def _rm_all(paths):
        c = 0
        for p in paths:
            if os.path.lexists(p):
                _rm(p)
                c += 1
        return c
    try:
        n = 0
        # (A) cram uploaded (.copied): the cram + its encrypted/index siblings are done
        #     -- their only consumers (Encrypt_crams/copy_to_dcache) are past, and the
        #     cram is NOT on the deepvariant/stats path (those read the bam). Safe even
        #     while the sample is still in flight, so the big crams free early instead
        #     of waiting for the whole sample to finish.
        for cop in glob.glob(os.path.join(CRAM, "*.mapped_hg38.cram.copied")):
            base = cop[:-len(".copied")]  # -> {sample}.mapped_hg38.cram
            n += _rm_all((base, base + ".crai", base + ".c4gh"))
        # (B) Per-sample .finished means the sample products exist, but cohort
        #     aggregation may still consume loose statistics and gVCFs. Keep the
        #     staged source/FASTQ/BAM recovery path intact until the durable
        #     pipeline-wide marker proves that aggregation and uploads completed.
        finished = glob.glob(os.path.join(SOURCEDIR, "*.finished"))
        if os.path.exists("pipeline.done"):
            for fin in finished:
                s = os.path.basename(fin)[:-len(".finished")]
                se = glob.escape(s)
                paths = [os.path.join(SOURCEDIR, s + ".data"),
                         os.path.join(SOURCEDIR, s + ".dcache_data"),
                         os.path.join(FQ_BADMAP, s + ".badmap.fastqs.tar.gz"),
                         os.path.join(STAT, s + ".stats.tar.gz")]
                for pat in (os.path.join(FQ, se + ".*.fq.gz"),
                            os.path.join(FQ_BADMAP, se + ".badmap.*.fastq.gz"),
                            os.path.join(FQ_BADMAP, se + ".*.badmap_*.fastq.gz"),
                            os.path.join(STAT, "cov", se + ".*"),
                            os.path.join(BAM, se + ".*.bam"),
                            os.path.join(BAM, se + ".*.bam.bai")):
                    paths.extend(glob.glob(pat))
                n += _rm_all(paths)
            # Samplefile/cohort archives are also declared temp(), but can be
            # stranded when a previous Snakemake invocation ended after their
            # upload consumer completed. pipeline.done proves every upload and
            # aggregation consumer has finished.
            aggregate_temp = []
            for pat in (
                os.path.join(STAT, "*.stats_bundle.tar.gz"),
                os.path.join(KRAKEN, "*.kraken_reports.tar.gz"),
                os.path.join(KRAKEN, "*.kraken_read_classification.tar.gz"),
                os.path.join(chrM, "tar", "chrM_gvcfs.tar.gz"),
                os.path.join(GVCF_TAR, "deepvariant_level2_*", "*.gvcf.tar"),
            ):
                aggregate_temp.extend(glob.glob(pat))
            n += _rm_all(aggregate_temp)
        elif finished:
            print("[onstart] preserving intermediates for %d finished sample(s): "
                  "pipeline.done is absent" % len(finished))
        if n:
            print("[onstart] reclaimed %d orphaned temp file(s)/dir(s)" % n)
    except Exception:
        pass

onsuccess: shell(# "rm -f zslurm-*"
                 # "rm -rf logs"
                 "rm -rf tmp")
onerror:
            shell("""
            sample_pattern="{sample_pattern}"
            # Error reporting must never mask the workflow's original failure.
            # In particular, grep exits with 1 when a valid search has no
            # matches; with Snakemake's pipefail setting that previously made
            # this onerror hook fail as a second, misleading exception.
            rm -f error_rules.txt error_samples.txt error.log
            {{ grep -r 'Error in rule' zslurm_logs/ | awk '{{print $1 "\t" $4}}' | awk -F"[/:]" '{{print$1 "\t" $2}}' | awk '{{print$1 "\t" $3}}' || true; }} > error_rules.txt
            {{ grep -r -A 2 'Error in rule' zslurm_logs/ | grep 'input' | awk -F[,] '{{print$1}}' | grep -E -o "$sample_pattern" || true; }} > error_samples.txt
            paste error_rules.txt error_samples.txt > error.log || true
            """)
