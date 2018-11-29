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

localrules: make_directories, md5sum, move_fastq1, move_fastq2

rule all: 
    input: 
	    OUTPUTDIR + '/qc/multiqc.html'

rule multiqc:
    input:
	    "input/multi_QC_tmp.txt"
    output:
        OUTPUTDIR + '/qc/multiqc.html'
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
	    OUTPUTDIR + '/qc/multiqc.log'
    wrapper:
        "0.27.1/bio/multiqc"
