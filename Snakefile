configfile: "config.yaml"  

import io 
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#

SAMPLES = pd.read_table(config["input_ena_table"])
INPUTDIR = config["inputDIR"]
run_accession = list(SAMPLES.run_accession)
ADAPTERS = config["adapters"]
STUDY = list(set(SAMPLES['study_accession'].tolist()))
SCRATCHDIR = config["scratch"]
OUTPUTDIR = config["outputDIR"]

#----QC DATA FILE----#

assert(len(STUDY)==1), 'This ena table contains more than one study accession' 
assert(STUDY[0]==config["study_accession"]), 'The study accession provided in the config file does not match the study accession provided in the ena table.'

pathlib.Path(OUTPUTDIR).mkdir(parents=True, exist_ok=True)

#----DEFINE RULES----#

localrules: multiqc 

rule all: 
    input:
        fastqcZIP = expand("{base}/qc/fastqc/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, sample=run_accession, num = [1,2]), 
        multiQC = OUTPUTDIR + '/qc/multiqc.html', 
        trimmedData = expand("{base}/trimmed/{sample}_{num}.trimmed.fastq.gz", base = OUTPUTDIR, sample=run_accession, num = [1,2])

rule fastqc:
    input:
        INPUTDIR + "/{sample}/{sample}_{num}.fastq.gz" 
    output:
        html = OUTPUTDIR + '/qc/fastqc/{sample}_{num}_fastqc.html', 
        zip = OUTPUTDIR + '/qc/fastqc/{sample}_{num}_fastqc.zip'
    params: ""
    log: 
        OUTPUTDIR + '/logs/fastqc/{sample}_{num}.log'
    wrapper:
        "0.27.1/bio/fastqc"

rule multiqc:
    input:
        expand("{base}/qc/fastqc/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, sample = run_accession, num = [1,2])
    output:
        html = OUTPUTDIR + "/qc/multiqc.html", 
        stats = OUTPUTDIR + "/qc/multiqc_general_stats.txt"
    conda: 
        "envs/multiqc-env.yaml"
    shell: 
        """
        multiqc -n multiqc.html {input}
        mv multiqc.html {output.html}
        mv multiqc_data/multiqc_general_stats.txt {output.stats} 
        rm -rf multiqc_data
        """

rule trimmomatic: 
    input:
        r1 = INPUTDIR + "/{sample}/{sample}_1.fastq.gz", 
        r2 = INPUTDIR + "/{sample}/{sample}_2.fastq.gz" 
    output:
        r1 = OUTPUTDIR + "/trimmed/{sample}_1.trimmed.fastq.gz",
        r2 = OUTPUTDIR + "/trimmed/{sample}_2.trimmed.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired = OUTPUTDIR + "/trimmed/{sample}_1.unpaired.fastq.gz",
        r2_unpaired = OUTPUTDIR + "/trimmed/{sample}_2.unpaired.fastq.gz"
    log:
        OUTPUTDIR +  "/logs/trimmomatic/{sample}.log"
    params:
        # UPDATE TRIMMING DETAILS -- SH 
        trimmer=[ "LEADING:20 MINLEN:25"],
        # EXTRA FLAGS? 
        extra=""
    wrapper:
        "0.27.1/bio/trimmomatic/pe" 
