rule all:
    input:
        "capture_kit_finder_done.txt",
        #"gatk_gcnv_done.txt",
        "kmc_done.txt",
        "align_fused_done.txt",
        "annovar_done.txt",
        "preprocess_done.txt",
        "kraken_done.txt",
        "PCA_done.txt",
        "qc_fused_done.txt",
        "vcf_handling_done.txt",
        "delly_done.txt",
        "picard_done.txt",
        "mosdepth_done.txt"

rule capture_kit_finder:
    conda: "../envs/capture_kit_finder.yaml"
    output: "capture_kit_finder_done.txt"
    shell: "echo 'Running capture_kit_finder' > {output}"

# This environment is not currently used, and its post-processing may be replaceable with `conda install bioconda::gcnvkernel`.
# rule gatk_gcnv:
#     conda: "../envs/gatk_gcnv.yaml"
#     output: "gatk_gcnv_done.txt"
#     shell: "echo 'Running gatk_gcnv' > {output}"

rule kmc:
    conda: "../envs/kmc.yaml"
    output: "kmc_done.txt"
    shell: "echo 'Running kmc' > {output}"

rule align_fused:
    conda: "../envs/align_fused.yaml"
    output: "align_fused_done.txt"
    shell: "echo 'Running align_fused' > {output}"

rule annovar:
    conda: "../envs/annovar.yaml"
    output: "annovar_done.txt"
    shell: "echo 'Running annovar' > {output}"

rule preprocess:
    conda: "../envs/preprocess.yaml"
    output: "preprocess_done.txt"
    shell: "echo 'Running preprocess' > {output}"

rule kraken:
    conda: "../envs/kraken.yaml"
    output: "kraken_done.txt"
    shell: "echo 'Running kraken' > {output}"

rule PCA:
    conda: "../envs/PCA.yaml"
    output: "PCA_done.txt"
    shell: "echo 'Running PCA' > {output}"

rule qc_fused:
    conda: "../envs/qc_fused.yaml"
    output: "qc_fused_done.txt"
    shell: "echo 'Running qc_fused' > {output}"

rule vcf_handling:
    conda: "../envs/vcf_handling.yaml"
    output: "vcf_handling_done.txt"
    shell: "echo 'Running vcf_handling' > {output}"

rule delly:
    conda: "../envs/delly.yaml"
    output: "delly_done.txt"
    shell: "echo 'Running delly' > {output}"

rule picard:
    conda: "../envs/picard.yaml"
    output: "picard_done.txt"
    shell: "echo 'Running picard' > {output}"

rule mosdepth:
    conda: "../envs/mosdepth.yaml"
    output: "mosdepth_done.txt"
    shell: "echo 'Running mosdepth' > {output}"