# Snakemake workflow for analysis of Tara Euk Metagenomes

## Setup

Create a conda environment to run pipeline
```
conda env create --name snakemake-tara-euk --file environment.yaml
```
Enter environment- change above to newest version, etc.


### Working directory
Provided tables & scripts include:
* ```tara-euk-metaG/input/ENA_tables```- ENA tables which are merged sample lists of the ENA download information and the Tara Ocean metadata from Pangea. [See Tara Ocean data download repo here](https://github.com/AlexanderLabWHOI/tara-download-snakemake).
* ```tara-euk-metaG/input```- List of samples to be combined for the assemblies. In the case of the Tara Oceans data this is by ocean region, province, and size fraction.

