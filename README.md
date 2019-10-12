# Snakemake workflow for analysis of Tara Euk Metagenomes

## Introduction
Molecular and genomic approaches, particularly those applied to whole, mixed communities (e.g. metagenomics, metatranscriptomics), have shed light on the ecological roles, evolutionary histories, and physiological capabilities of these organisms. We developed a scalable and reproducible pipeline to facilitate the retrieval, taxonomic assignment, and annotation of eukaryotic metagenome assembled genomes (MAGs) from mixed community metagenomes. The below pipeline uses **EukHeist** to retrieve eukaryotic (and other) metagenome assembled genomes from the *Tara Oceans* global dataset.

## Setup

Create a conda environment to run pipeline
```
conda env create --name snakemake-tara-euk --file environmentv0.2.yaml
```

### Working directory
Provided tables & scripts include:
* ```tara-euk-metaG/input/ENA_tables```- ENA tables which are merged sample lists of the ENA download information and the Tara Ocean metadata from Pangea. [See Tara Ocean data download repo here](https://github.com/AlexanderLabWHOI/tara-download-snakemake).
* ```tara-euk-metaG/input```- List of samples to be combined for the assemblies. In the case of the Tara Oceans data this is by ocean region, province, and size fraction.


### Notes to generalize Tara-specific pipeline

* See 2 remove directories where older files are placed - check this and remove. ```tara-euk-metaG/input/ENA_tables/rm-dir```
* environments? versioning eukheist?
