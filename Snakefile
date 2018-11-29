configfile: "config.yaml"  

import io 
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#

SAMPLES = pd.read_table(config["input_ena_table"])
SAMPLES[['fastq_ftp1', 'fastq_ftp2']] = SAMPLES.fastq_ftp.str.split(';', expand =True)
STUDY = list(set(SAMPLES['study_accession'].tolist()))

assert(len(STUDY)==1), 'This ena table contains more than one study accession' 
assert(STUDY[0]==config["study_accession"]), 'The study accession provided in the config file does not match the study accession provided in the ena table.'

STUDY=STUDY[0]
RUNS = SAMPLES['run_accession'].tolist()
SCRATCHDIR = config["scratch"]
OUTPUTDIR = config["outputDIR"]
DOWNLOAD_DICT = dict(zip(SAMPLES.run_accession, SAMPLES.fastq_ftp))

#----DEFINE RULES----#

rule all: 
    input: 
        read1 = expand("{outdir}/{study}/{run}/{run}_1.fastq.gz", outdir = OUTPUTDIR, study = STUDY, run = RUNS),
        read2 = expand("{outdir}/{study}/{run}/{run}_2.fastq.gz", outdir = OUTPUTDIR, study = STUDY, run = RUNS), 
        md5sum = expand("{outdir}/{study}/{run}/md5sum.tab", outdir = OUTPUTDIR, study = STUDY, run = RUNS)

localrules: make_directories, md5sum, move_fastq1, move_fastq2

rule make_directories:
    output: directory(expand("{outdir}/{study}/{run}/", outdir = OUTPUTDIR, study = STUDY, run=RUNS)) 

rule download_fastq:
    params: ftp= lambda wildcards: DOWNLOAD_DICT[wildcards.run].split(';')[0]
    output: SCRATCHDIR+'/'+STUDY+'/{run}/{run}_1.fastq.gz'   
    shell:
        """
        curl -L {params.ftp} --create-dirs --output {output} 
        """
rule download_fastq2:
    params: ftp= lambda wildcards: DOWNLOAD_DICT[wildcards.run].split(';')[1] 
    output: SCRATCHDIR+'/'+STUDY+'/{run}/{run}_2.fastq.gz'
    shell: 
        """
        curl -L {params.ftp} --create-dirs --output {output} 
        """

rule move_fastq1: 
    input:  SCRATCHDIR+'/'+STUDY+'/{run}/{run}_1.fastq.gz' 
    output: OUTPUTDIR+'/'+STUDY+'/{run}/{run}_1.fastq.gz' 
    shell: 'mv {input} {output}' 

rule move_fastq2: 
    input:  SCRATCHDIR+'/'+STUDY+'/{run}/{run}_2.fastq.gz' 
    output: OUTPUTDIR+'/'+STUDY+'/{run}/{run}_2.fastq.gz' 
    shell: 'mv {input} {output}'

rule md5sum:
    input: read1 = OUTPUTDIR+'/'+STUDY+'/{run}/{run}_1.fastq.gz',
            read2 = OUTPUTDIR+'/'+STUDY+'/{run}/{run}_2.fastq.gz'
    output: OUTPUTDIR+'/'+STUDY+'/{run}/md5sum.tab'
    shell: 
         """
         md5sum {input.read1} > {output} 
         md5sum {input.read2} >> {output}
        """
