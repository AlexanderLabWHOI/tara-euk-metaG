configfile: "config-test.yaml"  

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
SAMPLELIST = pd.read_table(config["sample_list"], index_col='Assembly_group')
ASSEMBLYGROUP = list(SAMPLELIST.index)
MEGAHIT_CPU = config["megahit_cpu"]
MEGAHIT_MIN_CONTIG = config["megahit_min_contig"]
MEGAHIT_MEM = config["megahit_mem"]


#----FUNCTIONS----#

def identify_read_groups(assembly_group_name, FORWARD=True):
    outlist=[] 
    ERR_list = SAMPLELIST.loc[assembly_group_name, 'ERR_list'].split(', ')
    if FORWARD: 
        num = 1
    else: 
        num = 2
    for E in ERR_list: 
        outlist.append(SCRATCHDIR + "/trimmed/{}_{}.trimmed.fastq.gz".format(E, num)) 
    return(outlist)

def get_sample_list(assembly_group):
    outlist = []
    for assembly_group_name in assembly_group: 
        ERR_list = SAMPLELIST.loc[assembly_group_name, 'ERR_list'].split(', ')
        for E in ERR_list:
            out = SCRATCHDIR + "/mapping/{}/{}.bam".format(assembly_group_name, E) 
            outlist.append(out)
    return(outlist)

def get_sample_list_onegroup(assembly_group_name):
    outlist = []
    ERR_list = SAMPLELIST.loc[assembly_group_name, 'ERR_list'].split(', ')
    for E in ERR_list:
        out = SCRATCHDIR + "/mapping/{}/{}.bam".format(assembly_group_name, E)
        outlist.append(out)
    return(outlist)        

#----QC DATA FILE----#

assert(len(STUDY)==1), 'This ena table contains more than one study accession' 
assert(STUDY[0]==config["study_accession"]), 'The study accession provided in the config file does not match the study accession provided in the ena table.'

pathlib.Path(OUTPUTDIR).mkdir(parents=True, exist_ok=True)

#----DEFINE RULES----#

localrules: multiqc, copy_bwa_index 

rule all: 
    input:
        # QC DATA
        fastqcZIP_raw = expand("{base}/qc/fastqc/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, sample=run_accession, num = [1,2]), 
        fastqcZIP_trimmed = expand("{base}/qc/fastqc/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, sample=run_accession, num = [1,2]), 
        multiQC_raw = OUTPUTDIR + '/qc/raw_multiqc.html', 
        multiQC_trimmed = OUTPUTDIR + '/qc/trimmed_multiqc.html', 
        #TRIM DATA
        trimmedData = expand("{base}/trimmed/{sample}_{num}.trimmed.fastq.gz", base = SCRATCHDIR, sample=run_accession, num = [1,2]), 
        #errtrimmedData = expand("{base}/errtrim/{sample}_{num}.trimmed.errtrim.fastq.gz", base = OUTPUTDIR, sample=run_accession, num = [1,2]),
        # #NORMALIZE DATA
        # normalizedData = expand("{base}/normalized/{sample}_{num}.trimmed.normalized.fastq.gz", base= OUTPUTDIR, sample = run_accession, num=[1,2]),
        #CALCULATE SOURMASH
        signature = expand("{base}/sourmash/{sample}.10k.sig", base = SCRATCHDIR, sample = run_accession),
        #ASSEMBLE
        assembly = expand("{base}/megahit/{assembly_group}/final.contigs.fa", base = OUTPUTDIR, assembly_group = ASSEMBLYGROUP),  
        assembly_copy = expand("{base}/bwa_index/{assembly_group}.fa", base = SCRATCHDIR, assembly_group = ASSEMBLYGROUP),  
        #BWA INDEX
        bwa_index = expand("{base}/bwa_index/{assembly_group}.fa.{bwa_tail}", base = SCRATCHDIR, assembly_group = ASSEMBLYGROUP, bwa_tail = ["amb", "ann", "bwt", "pac", "sa"]), 
        #BWA MAPPING:
        bwa_mem = get_sample_list(ASSEMBLYGROUP),  
        #METABAT2 
        jgi_abund = expand("{base}/metabat2/{assembly_group}/jgi_abund.txt", base = OUTPUTDIR, assembly_group = ASSEMBLYGROUP)
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
        r1 = SCRATCHDIR + "/trimmed/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{sample}_2.trimmed.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired = SCRATCHDIR + "/trimmed/{sample}_1.unpaired.fastq.gz",
        r2_unpaired = SCRATCHDIR + "/trimmed/{sample}_2.unpaired.fastq.gz"
    log:
        OUTPUTDIR +  "/logs/trimmomatic/{sample}.log"
    params:
        trimmer=["ILLUMINACLIP:{}:2:30:7".format(ADAPTERS), "LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:50"],
        extra=""
    wrapper:
        "0.27.1/bio/trimmomatic/pe"
 
rule fastqc_trimmed:
    input:
        SCRATCHDIR + "/trimmed/{sample}_{num}.trimmed.fastq.gz" 
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

rule compute_sigs:
    input:
        r1 = SCRATCHDIR + "/trimmed/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{sample}_2.trimmed.fastq.gz" 
    output: 
        SCRATCHDIR + "/sourmash/{sample}.10k.sig"
    conda: 
        "envs/sourmash.yaml"
    log:
         OUTPUTDIR +  "/logs/sourmash/sourmash_{sample}.log"
    shell: 
        """
        zcat {input.r1} {input.r2} | sourmash compute -k 21,31,51\
            --scaled 10000  --track-abundance \
            -o {output} - 2> {log}
        """


rule megahit_assembly: 
    input: r1 = lambda wildcards: identify_read_groups("{assembly_group}".format(assembly_group=wildcards.assembly_group)), 
           r2 = lambda wildcards: identify_read_groups("{assembly_group}".format(assembly_group=wildcards.assembly_group), FORWARD=False) 
    output: 
       OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa"  
    conda: 
        "envs/megahit.yaml"
    log: 
        OUTPUTDIR + "/logs/megahit/{assembly_group}.log" 
    params: 
        inputr1 = lambda wildcards, input: ','.join(input.r1),
        inputr2 = lambda wildcards, input: ','.join(input.r2),
        min_contig_len = MEGAHIT_MIN_CONTIG,
        cpu_threads = MEGAHIT_CPU, 
        memory = MEGAHIT_MEM, 
        other_options = "--continue", 
        megahit_output_name = lambda wildcards: "{}/megahit/{}".format(OUTPUTDIR, wildcards.assembly_group),
        megahit_output_prefix = lambda wildcards: "{}".format(wildcards.assembly_group),
    shell: 
        """
        megahit -1 {params.inputr1} -2 {params.inputr2} --min-contig-len {params.min_contig_len} --memory {params.memory} -t {params.cpu_threads} --out-dir {params.megahit_output_name} {params.other_options}  >> {log} 2>&1
        """

rule bwa_index:
    input:
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa"
    output:
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.amb",
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.ann",
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.bwt",
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.pac",
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.sa"
    log:
        OUTPUTDIR + "/logs/bwa_index/{assembly_group}.log"
    params:
        algorithm="bwtsw"
    conda:
        "envs/metabat-env.yaml"
 
    shell:
        """
        bwa index {input} 2> {log}
        """

rule copy_bwa_index:
    input: 
        OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa"
    output:
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa"
    shell:
        """
        cp {input} {output}
        """

rule bwa_mem:
    input:
        amb = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.amb", 
        ann = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.ann", 
        bwt = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.bwt", 
        pac = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.pac",  
        sa = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.sa", 
        reference = SCRATCHDIR + "/bwa_index/{assembly_group}.fa", 
        r1 = SCRATCHDIR + "/trimmed/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{sample}_2.trimmed.fastq.gz", 
    output:
        SCRATCHDIR + "/mapping/{assembly_group}/{sample}.bam"
    log:
        OUTPUTDIR + "/logs/bwa_mem/{assembly_group}/{sample}.log"
    params: 
        extra="", 
        #pipe_cmd = "samtools sort -o {output} -",  
        threads = 8
    conda: 
        "envs/metabat-env.yaml"
    shell:
        """ 
        bwa mem -t {params.threads} {params.extra} {input.reference} {input.r1} {input.r2} | samtools sort -o {output} - >> {log} 2>&1
        """ 
       
rule metabat_abundance:
    input:
        lambda wildcards: get_sample_list_onegroup("{assembly_group}".format(assembly_group=wildcards.assembly_group))
    output:
        OUTPUTDIR + "/metabat2/{assembly_group}/jgi_abund.txt"
    conda:
         "envs/metabat-env.yaml"
    log:
        OUTPUTDIR + "/logs/metabat2/{assembly_group}.abun.log"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input}
        """

rule metabat_binning:
    input:
        assembly = OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa",
        depth = OUTPUTDIR + "/metabat2/{assembly_group}/jgi_abund.txt"
    output:
        directory("")
    conda:
    params: 
        other = ""
    shell:
        """
        metabat2 {params.other} -i {input.assembly} -a {input.depth} -o {output}
        """
