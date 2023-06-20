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
    conda: "envs/gatk.yaml"   
    resources: 
        n = 1,
        mem_mb = 1000       
    shell:
        """{gatk} BedToIntervalList  --java-options "-Xmx{resources.mem_mb}m {params.java_options}"  -I {input.bed_file} -O {output.interval_list} -SD {params.seq_dict} --UNIQUE"""




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