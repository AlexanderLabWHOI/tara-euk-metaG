configfile: "config.yaml"  

import io 
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#

SAMPLES = pd.read_table(config["input_ena_table"])
INPUTDIR = config["inputDIR"]
# Make table for multiqc
run_accession = list(SAMPLES.run_accession)
multiQC_file = open('input/multi_QC_tmp.txt', 'w')
for i in run_accession:
    fastq1 = i+"_1.fastq.gz"
    fastq2 = i+"_2.fastq.gz"
    f1p = os.path.join(INPUTDIR, fastq1)
    f2p = os.path.join(INPUTDIR, fastq2)
    multiQC_file.write('\n'.join([f1p, f2p]))
    multiQC_file.write('\n')
ADAPTERS = config["adapters"]
STUDY = list(set(SAMPLES['study_accession'].tolist()))

assert(len(STUDY)==1), 'This ena table contains more than one study accession' 
assert(STUDY[0]==config["study_accession"]), 'The study accession provided in the config file does not match the study accession provided in the ena table.'

#STUDY=STUDY[0]
#RUNS = SAMPLES['run_accession'].tolist()
SCRATCHDIR = config["scratch"]
OUTPUTDIR = config["outputDIR"]
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
        #html = expand("{base}/qc/fastqc/{sample}_{num}.html", sample= run_accession, num=[1,2], base=OUTPUTDIR), 
        #zip =   expand("{base}/qc/fastqc/{sample}_{num}.zip", sample= run_accession, num=[1,2], base=OUTPUTDIR), 
        html = OUTPUTDIR + '/qc/fastqc/{sample}_{num}.html'
    params: ""
    log:
        #expand("{base}/qc/fastqc/log/{sample}_{num}.log", sample= run_accession, num=[1,2], base=OUTPUTDIR) 
        OUTPUTDIR + '/qc/fastqc/{sample}_{num}.log'
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
        OUTPUTDIR + '/qc/filelist.txt'    
    output:
        OUTPUTDIR + '/qc/multiqc.html'
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
	    OUTPUTDIR + '/qc/multiqc.log'
    wrapper:
        "0.27.1/bio/multiqc"
