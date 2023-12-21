from common import *

rule getallchrom:
    input:
        os.path.join('sex_chrom.tab'),
        os.path.join('bam_rg_quality.tab'),
        os.path.join('bam_quality.tab'),
        os.path.join('oxo_quality.tab')
    output:
        'result.done'
    shell: """
        touch {output}
        """

rule merge_sexstats:
    input:
        [os.path.realpath(samplefile + '.tsv')[:-4] + '.sex_chrom.tab' for samplefile in SAMPLE_FILES]
    output:
        os.path.join('sex_chrom.tab')
    run: 
        rin = " ".join(input)
        cmd = """
        awk 'FNR==1 && NR!=1 {{ getline; }} 
             1 {{print}}' {rin}  > {output}
        """
        shell(cmd)

rule merge_bamstats:
    input:
        [os.path.realpath(samplefile + '.tsv')[:-4] + '.bam_quality.tab' for samplefile in SAMPLE_FILES]
    output:
        os.path.join('bam_quality.tab')
    run: 
        rin = " ".join(input)
        cmd = """
        awk 'FNR==1 && NR!=1 {{ getline; }} 
             1 {{print}}' {rin}  > {output}
        """
        shell(cmd)

rule merge_bamrgstats:
    input:
        [os.path.realpath(samplefile + '.tsv')[:-4] + '.bam_rg_quality.tab' for samplefile in SAMPLE_FILES]
    output:
        os.path.join('bam_rg_quality.tab')
    run: 
        rin = " ".join(input)
        cmd = """
        awk 'FNR==1 && NR!=1 {{ getline; }} 
             1 {{print}}' {rin}  > {output}
        """
        shell(cmd)

rule merge_oxostats:
    input:
        [os.path.realpath(samplefile + '.tsv')[:-4] + '.oxo_quality.tab' for samplefile in SAMPLE_FILES]
    output:
        os.path.join('oxo_quality.tab')
    run: 
        rin = " ".join(input)
        cmd = """
        awk 'FNR==1 && NR!=1 {{ getline; }} 
             1 {{print}}' {rin}  > {output}
        """
        shell(cmd)

