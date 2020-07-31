configfile: "config.yaml"  
import glob
import io 
import os
import pandas as pd
import numpy as np
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#
METAG_ACCESSION = config["metaG_accession"]
METAT_ACCESSION = config["metaT_accession"]
METAG_SAMPLES = pd.read_table(config["metaG_ena_table"])
METAT_SAMPLES = pd.read_table(config["metaT_ena_table"])
INPUTDIR = config["inputDIR"]
ADAPTERS = config["adapters"]
METAG_STUDY = list(set(METAG_SAMPLES["study_accession"].tolist()))
METAT_STUDY = list(set(METAT_SAMPLES["study_accession"].tolist()))
SCRATCHDIR = config["scratch"]
OUTPUTDIR = config["outputDIR"]
METAG_SAMPLELIST = pd.read_table(config["metaG_sample_list"], index_col="Assembly_group")
METAG_ASSEMBLYGROUP = list(METAG_SAMPLELIST.index)
METAT_SAMPLELIST = pd.read_table(config["metaT_sample_list"], index_col="Assembly_group")
METAT_ASSEMBLYGROUP= list(METAT_SAMPLELIST.index)
ASSEMBLYGROUP = METAG_ASSEMBLYGROUP
MAG_SCRATCHDIR=config['scratch2']

#----COMPUTE VAR----#
MEGAHIT_CPU = config["megahit_cpu"]
MEGAHIT_MIN_CONTIG = config["megahit_min_contig"]
MEGAHIT_MEM = config["megahit_mem"]
MEGAHIT_OTHER = config["megahit_other"]

#----ACCESSION NUMS----#
metaG_run_accession = list(METAG_SAMPLES.run_accession)
metaT_run_accession = list(METAT_SAMPLES.run_accession)


#----FUNCTIONS----#

def identify_read_groups(assembly_group_name, STUDY, FORWARD=True):
    outlist=[] 
    ERR_list = METAG_SAMPLELIST.loc[assembly_group_name, 'ERR_list'].split(', ')
    if FORWARD: 
        num = 1
    else: 
        num = 2
    for E in ERR_list: 
        outlist.append(SCRATCHDIR + "/trimmed/{}/{}_{}.trimmed.fastq.gz".format(STUDY,E, num)) 
    return(outlist)

def get_sample_list(assembly_group, SAMPLELIST, STUDY):
    outlist = []
    for assembly_group_name in assembly_group: 
        ERR_list = SAMPLELIST.loc[assembly_group_name, 'ERR_list'].split(', ')
        for E in ERR_list:
            out = SCRATCHDIR + "/mapping/{}/{}/{}.bam".format(assembly_group_name,STUDY, E) 
            outlist.append(out)
    return(outlist)


def get_sample_list_onegroup(assembly_group_name, SAMPLELIST, STUDY):
    outlist = []
    ERR_list = SAMPLELIST.loc[assembly_group_name, 'ERR_list'].split(', ')
    for E in ERR_list:
        out = SCRATCHDIR + "/mapping/{}/{}/{}.bam".format(assembly_group_name, STUDY, E)
        outlist.append(out)
    return(outlist)        

#----QC DATA FILE----#

assert(len(METAT_STUDY)==1), 'This metaT ena table contains more than one study accession'
assert(len(METAG_STUDY)==1), 'This metaG ena table contains more than one study accession' 
assert(METAG_STUDY[0]==METAG_ACCESSION), 'The study accession provided in the config file does not match the study accession provided in the metaG ena table.'
assert(METAT_STUDY[0]==METAT_ACCESSION), 'The study accession provided in the config file does not match the study accession provided in the ena table.'
#assert(np.all(METAG_SAMPLELIST.index == METAT_SAMPLELIST.index)), 'The sample list provided for the MetaG does not match the MetaT'

pathlib.Path(OUTPUTDIR).mkdir(parents=True, exist_ok=True)

#----DEFINE RULES----#

localrules: multiqc, copy_bwa_index 

rule all: 
    input:
        # QC DATA
        #fastqcZIP_rawG = expand("{base}/qc/fastqc/{study}/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, study = METAG_ACCESSION, sample=metaG_run_accession, num = [1,2]),
        #fastqcZIP_rawT = expand("{base}/qc/fastqc/{study}/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, study = METAT_ACCESSION, sample=metaT_run_accession, num = [1,2]),  
        #fastqcZIP_trimmedG = expand("{base}/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, study = METAG_ACCESSION, sample=metaG_run_accession, num = [1,2]),  
        #fastqcZIP_trimmedT = expand("{base}/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, study = METAT_ACCESSION, sample=metaT_run_accession, num = [1,2]),  
        #MULTIQC
#        html_rawG = OUTPUTDIR + "/qc/rawG_multiqc.html",
#        stats_rawG = OUTPUTDIR + "/qc/rawG_multiqc_general_stats.txt",
#        html_trimmedG = OUTPUTDIR + "/qc/trimmedG_multiqc.html",
#        stats_trimmedG = OUTPUTDIR + "/qc/trimmedG_multiqc_general_stats.txt",
#        html_rawT = OUTPUTDIR + "/qc/rawT_multiqc.html",
#        stats_rawT = OUTPUTDIR + "/qc/rawT_multiqc_general_stats.txt",
#        html_trimmedT = OUTPUTDIR + "/qc/trimmedT_multiqc.html",
#        stats_trimmedT = OUTPUTDIR + "/qc/trimmedT_multiqc_general_stats.txt",
        #TRIM DATA
        trimmedDataG = expand("{base}/trimmed/{study}/{sample}_{num}.trimmed.fastq.gz", base = SCRATCHDIR, study = METAG_ACCESSION, sample=metaG_run_accession, num = [1,2]), 
        trimmedDataT = expand("{base}/trimmed/{study}/{sample}_{num}.trimmed.fastq.gz", base = SCRATCHDIR, study = METAT_ACCESSION, sample=metaT_run_accession, num = [1,2]), 
        #K-MER TRIMMED DATA
        #kmerTrimmedDataG = expand("{base}/abundtrim/{study}/{sample}.abundtrim.fastq.gz", base = SCRATCHDIR, study = METAG_ACCESSION, sample=metaG_run_accession), 

        #CALCULATE SOURMASH
        #signatureG = expand("{base}/sourmash/{study}/{sample}.10k.sig", base = SCRATCHDIR, study = METAG_ACCESSION, sample = metaG_run_accession),
        #signatureT = expand("{base}/sourmash/{study}/{sample}.10k.sig", base = SCRATCHDIR, study = METAT_ACCESSION, sample = metaT_run_accession),

        #ASSEMBLE
        #assembly = expand("{base}/megahit/{assembly_group}/final.contigs.fa", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP),  
        #assembly_copy = expand("{base}/bwa_index/{assembly_group}.fa", base = SCRATCHDIR, assembly_group = METAG_ASSEMBLYGROUP),  

        #BWA INDEX
        #bwa_index = expand("{base}/bwa_index/{assembly_group}.fa.{bwa_tail}", base = SCRATCHDIR, assembly_group = METAG_ASSEMBLYGROUP, bwa_tail = ["amb", "ann", "bwt", "pac", "sa"]), 

        #BWA MAPPING:
        #bwa_memG = get_sample_list(METAG_ASSEMBLYGROUP, METAG_SAMPLELIST, METAG_ACCESSION), 
        #bwa_memT = get_sample_list(METAT_ASSEMBLYGROUP, METAT_SAMPLELIST, METAT_ACCESSION), 

        #BINNING 

        #METABAT2 
        #jgi_abund = expand("{base}/metabat2/{assembly_group}/jgi_abund.txt", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP),
        #metabat2_bins = expand("{base}/metabat2/{assembly_group}/{assembly_group}_bin", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP),
        #CONCOCT
        concoct = expand(os.path.join(OUTPUTDIR, 'concoct', '{assembly_group}', 'clustering_gt2000.csv'), assembly_group = METAG_ASSEMBLYGROUP),
        ccfile=expand(os.path.join(OUTPUTDIR, 'concoct', '{assembly_group}','done'), assembly_group=METAG_ASSEMBLYGROUP),
        #EUKREP
        #eukrep =  expand("{base}/eukrep/{assembly_group}/euk.final.contigs.fa", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP), 
        #metabat2_bins_euk = expand("{base}/metabat2_euk/{assembly_group}/{assembly_group}_eukbin", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP),
        assignments = expand(OUTPUTDIR + "/metabat2/{assembly_group}/assignment.csv", assembly_group = METAG_ASSEMBLYGROUP), 
        #PROTEIN PREDICITION
        #PRODIGAL
        #proteins = expand("{base}/prodigal/{assembly_group}/proteins.faa", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP), 
#        MAG_MAPPING = expand("{base}/high-quality-mapping/{study}/{sample}.bam",base=MAG_SCRATCHDIR, study = METAG_ACCESSION, sample = metaG_run_accession),
#        mag_index = expand("{base}/high-quality-mags/hq-mags-concatenated.fa.{bwa_tail}", base = MAG_SCRATCHDIR, assembly_group = METAG_ASSEMBLYGROUP, bwa_tail = ["amb", "ann", "bwt", "pac", "sa"]),
#        MAG_abundance = expand("{base}/high-quality-mag-abundance/{study}/{sample}.coverm.abundance.tab", base = OUTPUTDIR, study = METAG_ACCESSION, sample = metaG_run_accession), 
#        MAG_hist = expand("{base}/high-quality-mag-abundance/{study}/{sample}.coverm.histogram.tab", base = OUTPUTDIR, study = METAG_ACCESSION, sample = metaG_run_accession) 
rule fastqc:
    input:
        INPUTDIR + "/{study}/{sample}/{sample}_{num}.fastq.gz"     
    output:
        html = OUTPUTDIR + '/qc/fastqc/{study}/{sample}_{num}_fastqc.html', 
        zip = OUTPUTDIR + '/qc/fastqc/{study}/{sample}_{num}_fastqc.zip'
    params: ""
    log: 
        OUTPUTDIR + '/logs/fastqc/{study}/{sample}_{num}.log'
    wrapper:
        "0.27.1/bio/fastqc"

rule trimmomatic: 
    input:
        r1 = INPUTDIR + "/{study}/{sample}/{sample}_1.fastq.gz", 
        r2 = INPUTDIR + "/{study}/{sample}/{sample}_2.fastq.gz" 
    output:
        r1 = SCRATCHDIR + "/trimmed/{study}/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{study}/{sample}_2.trimmed.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired = SCRATCHDIR + "/trimmed/{study}/{sample}_1.unpaired.fastq.gz",
        r2_unpaired = SCRATCHDIR + "/trimmed/{study}/{sample}_2.unpaired.fastq.gz"
    log:
        OUTPUTDIR +  "/logs/trimmomatic/{study}/{sample}.log"
    params:
        trimmer=["ILLUMINACLIP:{}:2:30:7".format(ADAPTERS), "LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:50"],
        extra=""
    wrapper:
        "0.27.1/bio/trimmomatic/pe"

rule fastqc_trimmed:
    input:
        SCRATCHDIR + "/trimmed/{study}/{sample}_{num}.trimmed.fastq.gz" 
    output:
        html = OUTPUTDIR + '/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.html', 
        zip = OUTPUTDIR + '/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip'
    params: ""
    log: 
        OUTPUTDIR + '/logs/fastqc/{study}/{sample}_{num}.trimmed.log'
    wrapper:
        "0.27.1/bio/fastqc"

rule multiqc:
    input:
        rawG = expand("{base}/qc/fastqc/{study}/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, study = METAG_ACCESSION, sample = metaG_run_accession, num = [1,2]), 
        trimmedG = expand("{base}/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, study = METAG_ACCESSION, sample = metaG_run_accession, num = [1,2]), 
        rawT = expand("{base}/qc/fastqc/{study}/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, study = METAT_ACCESSION, sample = metaT_run_accession, num = [1,2]), 
        trimmedT = expand("{base}/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, study = METAT_ACCESSION, sample = metaT_run_accession, num = [1,2]) 
    output:
        html_rawG = OUTPUTDIR + "/qc/rawG_multiqc.html", 
        stats_rawG = OUTPUTDIR + "/qc/rawG_multiqc_general_stats.txt",
        html_trimmedG = OUTPUTDIR + "/qc/trimmedG_multiqc.html", 
        stats_trimmedG = OUTPUTDIR + "/qc/trimmedG_multiqc_general_stats.txt",
        html_rawT = OUTPUTDIR + "/qc/rawT_multiqc.html", 
        stats_rawT = OUTPUTDIR + "/qc/rawT_multiqc_general_stats.txt",
        html_trimmedT = OUTPUTDIR + "/qc/trimmedT_multiqc.html", 
        stats_trimmedT = OUTPUTDIR + "/qc/trimmedT_multiqc_general_stats.txt"
    conda: 
        "envs/multiqc-env.yaml"
    shell: 
        """
        multiqc -n multiqc.html {input.rawG}
        mv multiqc.html {output.html_rawG}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_rawG} 
        rm -rf multiqc_data

        multiqc -n multiqc.html {input.trimmedG}
        mv multiqc.html {output.html_trimmedG}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_trimmedG} 
        rm -rf multiqc_data
        
        multiqc -n multiqc.html {input.trimmedT}
        mv multiqc.html {output.html_trimmedT}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_trimmedT} 
        rm -rf multiqc_data
        
        multiqc -n multiqc.html {input.trimmedT}
        mv multiqc.html {output.html_trimmedT}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_trimmedT} 
        rm -rf multiqc_data
        """ 


rule kmer_trim_reads:
    input:
        r1 = SCRATCHDIR + "/trimmed/{study}/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{study}/{sample}_2.trimmed.fastq.gz" 
    output: 
        SCRATCHDIR + "/abundtrim/{study}/{sample}.abundtrim.fastq.gz"
    conda: 
        "envs/trim_low_abund.yaml"
    log:
         OUTPUTDIR +  "/logs/abundtrim/{study}/trim_low_abund_{sample}.log"
    shell: 
        """
        interleave-reads.py {input.r2} {input.r1} | trim-low-abund.py --gzip -C 3 -Z 18 -M 30e9 -V - -o {output} 2> {log} 
        """

rule compute_sigs:
    input:
        r1 = SCRATCHDIR + "/trimmed/{study}/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{study}/{sample}_2.trimmed.fastq.gz" 
    output: 
        SCRATCHDIR + "/sourmash/{study}/{sample}.10k.sig"
    conda: 
        "envs/sourmash.yaml"
    log:
         OUTPUTDIR +  "/logs/sourmash/{study}/sourmash_{sample}.log"
    shell: 
        """
        zcat {input.r1} {input.r2} | sourmash compute -k 21,31,51\
            --scaled 10000  --track-abundance \
            -o {output} - 2> {log}
        """

rule megahit_assembly: 
    input: r1 = lambda wildcards: identify_read_groups("{assembly_group}".format(assembly_group=wildcards.assembly_group), METAG_ACCESSION), 
           r2 = lambda wildcards: identify_read_groups("{assembly_group}".format(assembly_group=wildcards.assembly_group), METAG_ACCESSION, FORWARD=False) 
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
        other_options = MEGAHIT_OTHER,
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
        r1 = SCRATCHDIR + "/trimmed/{study}/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{study}/{sample}_2.trimmed.fastq.gz", 
    output:
        SCRATCHDIR + "/mapping/{assembly_group}/{study}/{sample}.bam"
    log:
        OUTPUTDIR + "/logs/bwa_mem/{assembly_group}/{study}/{sample}.log"
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
        lambda wildcards: get_sample_list_onegroup("{assembly_group}".format(assembly_group=wildcards.assembly_group), METAG_SAMPLELIST, METAG_ACCESSION)
    output:
        OUTPUTDIR + "/metabat2/{assembly_group}/jgi_abund.txt"
    conda:
         "envs/metabat-env.yaml"
    log:
        OUTPUTDIR + "/logs/metabat2/{assembly_group}.abun.log"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input} > {log} 2>&1
        """

rule metabat_binning:
    input:
        assembly = OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa",
        depth = OUTPUTDIR + "/metabat2/{assembly_group}/jgi_abund.txt"
    output:
        OUTPUTDIR + "/metabat2/{assembly_group}/{assembly_group}_bin"
    conda:
         "envs/metabat-env.yaml"
    params: 
        other = "--saveCls",
        threads = 8
    log:
        OUTPUTDIR + "/logs/metabat2/{assembly_group}.bin.log"
    shell:
        """
        metabat2 {params.other} --numThreads {params.threads} -i {input.assembly} -a {input.depth} -o {output} > {log} 2>&1
        """

rule concoct_binning: 
    input: 
        assembly = OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa",
        depth = OUTPUTDIR + "/metabat2/{assembly_group}/jgi_abund.txt"
    output: 
        os.path.join(OUTPUTDIR, 'concoct', '{assembly_group}', 'clustering_gt2000.csv')
    params:
        outdir = os.path.join(OUTPUTDIR, 'concoct', '{assembly_group}'), length = 2000
    threads: 12
    conda: 'envs/concoct.yaml'
    shell:
        '''
        concoct --coverage_file {input.depth} \
            --composition_file {input.assembly} \
            --basename {params.outdir} \
            --length_threshold {params.length} \
            --converge_out -t 12\
        '''        

rule merge_concoct:
    input: os.path.join(OUTPUTDIR, 'concoct', '{assembly_group}', 'clustering_gt2000.csv')
    output: os.path.join(OUTPUTDIR, 'concoct', '{assembly_group}', 'clustering_merged.csv')
    conda: 'envs/concoct.yaml'
    shell: 
        '''
        merge_cutup_clustering.py {input} > {output}
        '''
rule extract_concoct_bins:
    input: assembly = OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa",
           csv = os.path.join(OUTPUTDIR, 'concoct', '{assembly_group}', 'clustering_merged.csv')
    output: done=os.path.join(OUTPUTDIR, 'concoct', '{assembly_group}','done')
    params: outdir= directory(os.path.join(OUTPUTDIR, 'concoct', '{assembly_group}','fastabins'))
    conda: 'envs/concoct.yaml'
    shell:
        '''
        mkdir -p {output}
        extract_fasta_bins.py {input.assembly} {input.csv} --output_path {params.outdir}
        touch {output.done}
        '''

rule eukcc_filter: 
    input: fastas= lambda wildcards: glob.glob(os.path.join(OUTPUTDIR, 'metabat2', '{assembly_group}','{assembly_group}*.fa').format(assembly_group=wildcards.assembly_group))
    output: OUTPUTDIR + "/metabat2/{assembly_group}/assignment.csv"
    conda: "envs/eukcc.yaml"
    params: minbpeuks=1000000, eukratio=0.2, tempdir = OUTPUTDIR + "/metabat2/{assembly_group}/tmp"
    shell:
        """
        ulimit -s 65536
        filter_euk_bins.py {input.fastas} --minbpeuks {params.minbpeuks} --eukratio {params.eukratio} --output {output} --tempdir {params.tempdir}
        """

rule eukrep:
    input: 
        assembly = OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa",
    output:
        OUTPUTDIR + "/eukrep/{assembly_group}/euk.final.contigs.fa"
    conda: 
        "envs/EukRep.yaml"
    log:
        OUTPUTDIR + "/logs/eukrep/{assembly_group}.eukrep.log"
    params: 
        prok = OUTPUTDIR + "/eukrep/{assembly_group}/prok.final.contigs.fa",
        min_contig = 1000
    shell: 
        """
        EukRep -i {input} -o {output} --prokarya {params.prok} --min {params.min_contig} > {log} 2>&1
        """

rule metabat_binning_euk:
    input:
        assembly = OUTPUTDIR + "/eukrep/{assembly_group}/euk.final.contigs.fa",
        depth = OUTPUTDIR + "/metabat2/{assembly_group}/jgi_abund.txt"
    output:
        OUTPUTDIR + "/metabat2_euk/{assembly_group}/{assembly_group}_eukbin"
    conda:
         "envs/metabat-env.yaml"
    params:
        other = "--saveCls",
        threads = 8
    log:
        OUTPUTDIR + "/logs/metabat2/{assembly_group}.eukbin.log"
    shell:
        """
        metabat2 {params.other} --numThreads {params.threads} -i {input.assembly} -a {input.depth} -o {output} > {log} 2>&1
        """

rule prodigal:
    input:
        assembly = OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa",
    output: 
        proteins = OUTPUTDIR + "/prodigal/{assembly_group}/proteins.faa",
        genes = OUTPUTDIR + "/prodigal/{assembly_group}/genes.gff"
    conda:
        "envs/prodigal.yaml"
    shell:
        """
        prodigal -i {input.assembly} -f gff -o {output.genes} -a {output.proteins} -p meta
        """

rule build_mag_index:
    input: MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa"
    output:
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.amb", 
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.ann",
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.bwt",
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.pac",
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.sa"
    log:
        OUTPUTDIR + "/logs/high-quality-mags/bwa-mags.log"
    params:
        algorithm="bwtsw"
    conda:
        "envs/metabat-env.yaml"

    shell:
        """
        bwa index {input} 2> {log}
        """

rule mag_bwa_mem:
    input:
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.amb", 
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.ann",
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.bwt",
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.pac",
        MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa.sa",
        reference = MAG_SCRATCHDIR + "/high-quality-mags/hq-mags-concatenated.fa", 
        r1 = SCRATCHDIR + "/trimmed/{study}/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{study}/{sample}_2.trimmed.fastq.gz", 
    output:
        MAG_SCRATCHDIR + "/high-quality-mapping/{study}/{sample}.bam"
    log:
        OUTPUTDIR + "/logs/high-quality-mapping/{study}/{sample}.log"
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
rule coverm_genome:
    input:
        mapping = MAG_SCRATCHDIR + "/high-quality-mapping/{study}/{sample}.bam", 
        genome_dir = OUTPUTDIR + "/high-quality-mags/"
        
    output:
        OUTPUTDIR + "/high-quality-mag-abundance/{study}/{sample}.coverm.abundance.tab" 
    conda:
        "envs/coverm.yaml"
    shell:
        """
        coverm genome --bam-files {input.mapping} --genome-fasta-directory {input.genome_dir} --genome-fasta-extension "fa"  --min-read-percent-identity 0.95     --min-read-aligned-percent 0.75  --proper-pairs-only --methods count length covered_bases covered_fraction reads_per_base mean variance trimmed_mean rpkm relative_abundance     --output-format dense     --min-covered-fraction 0     --contig-end-exclusion 75     --trim-min 0.05     --trim-max 0.95 --quiet > {output}
        """



rule coverm_genome_histogram:
    input:
        mapping = MAG_SCRATCHDIR + "/high-quality-mapping/{study}/{sample}.bam", 
        genome_dir = OUTPUTDIR + "/high-quality-mags/"
        
    output:
        OUTPUTDIR + "/high-quality-mag-abundance/{study}/{sample}.coverm.histogram.tab" 
    conda:
        "envs/coverm.yaml"
    shell:
        """
        coverm genome --bam-files {input.mapping} --genome-fasta-directory {input.genome_dir} --genome-fasta-extension "fa"  --min-read-percent-identity 0.95     --min-read-aligned-percent 0.75  --proper-pairs-only --methods coverage_histogram --output-format dense     --min-covered-fraction 0     --contig-end-exclusion 75     --trim-min 0.05     --trim-max 0.95 --quiet > {output}
        """


