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
        multiQC = OUTPUTDIR + '/qc/multiqc.html', 
        fastqcHTML = expand("{base}/qc/fastqc/{sample}_{num}.html", base = OUTPUTDIR, sample=run_accession, num = [1,2])

rule fastqc:
    input:
        expand( "{base}/{sample}/{sample}_{num}.fastq.gz", base=INPUTDIR, sample=run_accession, num = [1,2]) 
    output:
        html = OUTPUTDIR + '/qc/fastqc/{sample}_{num}.html', 
        zip = OUTPUTDIR + '/qc/fastqc/{sample}_{num}.zip'
    params: ""
    log: 
        OUTPUTDIR + '/qc/fastqc/log/{sample}_{num}.log'
    wrapper:
        "0.27.1/bio/fastqc"

rule make_multiQC: 
    input: 
        html = expand("{base}/qc/fastqc/{sample}_{num}.html", sample= run_accession, num=[1,2], base=OUTPUTDIR),  
    output: 
        OUTPUTDIR+'/qc/filelist.txt'
    shell: """
        echo "{OUTPUTDIR}/qc/fastqc/"> {OUTPUTDIR}/qc/filelist.txt
        """ 
rule multiqc:
    input:
        expand("{base}/qc/fastqc/{sample}_{num}.zip", base = OUTPUTDIR, sample = run_accession, num = [1,2])
    output:
        report = OUTPUTDIR + '/qc/multiqc.html'
    params:
        qc_dir = OUTPUTDIR + '/qc/fastqc/' # Optional: extra parameters for multiqc.
    log:
	    OUTPUTDIR + '/logs/multiqc.log'
    shell: 
        """
        multiqc -f -o {OUTPUTDIR}/qc/multiqc.html {params.qc_dir} 2> {log} 1>&2  
        """
