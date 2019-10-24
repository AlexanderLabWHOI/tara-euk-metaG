# Snakemake workflow for analysis of Tara Euk Metagenomes

## Introduction
Molecular and genomic approaches, particularly those applied to whole, mixed communities (e.g. metagenomics, metatranscriptomics), have shed light on the ecological roles, evolutionary histories, and physiological capabilities of these organisms. We developed a scalable and reproducible pipeline to facilitate the retrieval, taxonomic assignment, and annotation of eukaryotic metagenome assembled genomes (MAGs) from mixed community metagenomes. The below pipeline uses **EukHeist** to retrieve eukaryotic (and other) metagenome assembled genomes from the *Tara Oceans* global dataset.

## Setup

Create a conda environment to run pipeline
```
conda env create --name EukHeist --file environmentv0.2.yaml
```
* conda v4.7.10
* snakemake v5.6.0


### Working directory
Provided tables & scripts include:
* ```tara-euk-metaG/input/ENA_tables```- ENA tables which are merged sample lists of the ENA download information and the Tara Ocean metadata from Pangea. [See Tara Ocean data download repo here](https://github.com/AlexanderLabWHOI/tara-download-snakemake).
* ```tara-euk-metaG/input```- List of samples to be combined for the assemblies. In the case of the Tara Oceans data this is by ocean region, province, and size fraction.
* ```tara-euk-metaG/input/Curate-Input-Data-TARA.ipynb```- jupyter notebook detailing how ENA_tables with sample information were curated (e.g., revised size fractionation annotation) and compiled into _Assembly groups_.
* There are also test scripts and files which include the code to generate the test files and the output test sample names. This is a subset of the data which can be run as a test batch. 
* ```tara-euk-metaG/input/adapters``` - illumina adapters that are used during the trimmomatic (adapter trimming) step.
* ```cluster.yaml``` and ```submit_script/``` are configuration file and shell commands to run snakemake pipeline with Slurm.
* ```envs/``` list of conda environments snakemake uses throughout the pipeline


### Download TARA data

See [repo to download raw Tara Ocean data here](https://github.com/AlexanderLabWHOI/tara-download-snakemake).

## Alternatively, set up test dataset to run EukHeist pipeline.

Test files were created using *Curate-Test-input-data-TARA.ipynb*, but subsetting 20 metagenomic and 10 metatranscriptomic samples. 
* ```tara-euk-metaG/input/SampleList_ForAssembly_meta*_python-TEST.txt```
* ```tara-euk-metaG/input/ENA_tables/PRJEB*_meta*_wenv_PE-TEST.txt```

### Output files and descriptions

***
### Troubleshooting snakemake
When throwing an error, snakemake will list log files. Each time snakemake is executed, a log file is created in ```CURRENT_SNAKEMAKE_DIR/.snakemake/log/```. These are dated and provide the printed output. Some common errors and steps to diagnose.   

*Compatibility with snakemake and conda* Note the version of snakemake and anaconda you are using. Upon conda's switch from _source activate_ to _conda activate_, snakemake was not calling on the conda environments properly. Errors associted with these were ```returned non-zero exit status 127``` and an error about *line 56* in *thread.py* like this: ```LOCATION_OF_YOUR_CONDA_ENVS/.conda/envs/snake-18S/lib/python3.6/concurrent/futures/thread.py", line 56, in run```
Update your version of snakemake. Versions listed above are compatible. This error will also be generated when there is an incomptaible conda environment called, such as an outdated [snakemake wrapper](https://snakemake-wrappers.readthedocs.io/en/stable/). Try updating that - but ensure you have first updated the snakemake version.   

Check all log files, found in ```./snakemake/log``` and the log files generated as output. Most helpful are the output files from each slurm job, ```slurm-XXX.out```. Look at the most recent slurm out files after a job fails to find more informative errors.    

See ```snakemake -h``` for additional commands to clean up working snakemake environment, list steps, or restarting attempts, etc. When running, snakemake "locks" the working directory, so only one snakemake command can be run at a time. If you have to cancel, make sure to run ```snakemake --unlock``` to clear it. See other flags to clean up your working environment after updated conda, snakemake, or other environment versions (```--cleanup-shadow```, ```--cleanup-conda```).
To look for additional code error that may result in Syntax errors, try adding this to snakemake execution:
* ```--summary``` or ```--detailed-summary```
* ```--printshellcmds```
* ```--debug```


### to do
* include shell script to download only the test data?
* determine if updated snakemake version would be better
* figure out what is going on with multiqc
* list output files

### Notes to generalize Tara-specific pipeline

* See 2 remove directories where older files are placed - check this and remove. ```tara-euk-metaG/input/ENA_tables/rm-dir```
* separate directories for metaT and metaG - rather than PRNJ IDs
* environments? versioning eukheist?
