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
        # QC DATA
        fastqcZIP_raw = expand("{base}/qc/fastqc/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, sample=run_accession, num = [1,2]), 
        fastqcZIP_trimmed = expand("{base}/qc/fastqc/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, sample=run_accession, num = [1,2]), 
        multiQC_raw = OUTPUTDIR + '/qc/raw_multiqc.html', 
        multiQC_trimmed = OUTPUTDIR + '/qc/trimmed_multiqc.html', 
        #TRIM DATA
        trimmedData = expand("{base}/trimmed/{sample}_{num}.trimmed.fastq.gz", base = OUTPUTDIR, sample=run_accession, num = [1,2]), 
        #NORMALIZE DATA
        normalizedData = expand("{base}/normalized/{sample}_{num}.trimmed.normalized.fastq.gz", base= OUTPUTDIR, sample = run_accession, num=[1,2]),

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
        trimmer=[ "ILLUMINACLIP:ADAPTERS:2:30:7 LEADING:2 TRAILING:2 \
                SLIDINGWINDOW:4:2 MINLEN:50"],
        extra=""
    wrapper:
        "0.27.1/bio/trimmomatic/pe"
 
rule fastqc_trimmed:
    input:
        OUTPUTDIR + "/trimmed/{sample}_{num}.trimmed.fastq.gz" 
    output:
        html = OUTPUTDIR + '/qc/fastqc/{sample}_{num}.trimmed_fastqc.html', 
        zip = OUTPUTDIR + '/qc/fastqc/{sample}_{num}.trimmed_fastqc.zip'
    params: ""
    log: 
        OUTPUTDIR + '/logs/fastqc/{sample}_{num}.trimmed.log'
    wrapper:
        "0.27.1/bio/fastqc"

rule multiqc:
    input:
        raw = expand("{base}/qc/fastqc/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, sample = run_accession, num = [1,2]), 
        trimmed = expand("{base}/qc/fastqc/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, sample = run_accession, num = [1,2]) 
    output:
        html_raw = OUTPUTDIR + "/qc/raw_multiqc.html", 
        stats_raw = OUTPUTDIR + "/qc/raw_multiqc_general_stats.txt",
        html_trimmed = OUTPUTDIR + "/qc/trimmed_multiqc.html", 
        stats_trimmed = OUTPUTDIR + "/qc/trimmed_multiqc_general_stats.txt"
    conda: 
        "envs/multiqc-env.yaml"
    shell: 
        """
        multiqc -n multiqc.html {input.raw}
        mv multiqc.html {output.html_raw}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_raw} 
        rm -rf multiqc_data

        multiqc -n multiqc.html {input.trimmed}
        mv multiqc.html {output.html_trimmed}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_trimmed} 
        rm -rf multiqc_data
        """
rule normalized_samples: 
    input: 
        r1 = OUTPUTDIR + "/trimmed/{sample}_1.trimmed.fastq.gz",   
        r2 = OUTPUTDIR + "/trimmed/{sample}_2.trimmed.fastq.gz"
    output: 
        r1 = OUTPUTDIR + "/normalized/{sample}_1.trimmed.normalized.fastq.gz",
        r2 = OUTPUTDIR + "/normalized/{sample}_2.trimmed.normalized.fastq.gz", 
        inhist = OUTPUTDIR + "/normalized/{sample}.inhist",
        outhist = OUTPUTDIR + "/normalized/{sample}.outhist",
    params: 
        bb_targetDepth = 30, 
        bb_minDepth = 2,  
        bb_otherparams = "", 
        bb_threads = "", 
        r1 = OUTPUTDIR + "/normalized/{sample}_1.trimmed.normalized.fastq",
        r2 = OUTPUTDIR + "/normalized/{sample}_2.trimmed.normalized.fastq",
    conda: 
        'envs/bbmap-env.yaml'
    shell: 
        """ 
        bbnorm.sh in1={input.r1} in2={input.r2} out={params.r1} out2={params.r2} target={params.bb_targetDepth} min={params.bb_minDepth} hist={output.inhist} histout={output.outhist} {params.bb_otherparams} 
        pigz {params.r1}
        pigz {params.r2}
        """
