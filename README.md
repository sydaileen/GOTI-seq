# GOTI-seq

GOTI-seq is the bioinformatic pipeline for analyzing genome-wide sequencing data from GOTI experiments.

## Install

###### Requirements

```
The implement of GOTI-seq requires the following softwares installed on the system:

fastQC (v0.11.3)
Trimmomatic (v0.36) 
BWA (0.7.12) 
picard-tools (v2.3.0) 
GATK (v3.5) 
Strelka (v2.7.1) 
Lofreq (v2.1.2) 
Scalpel (v0.5.3) 
Annovar (version 2016-02-01) 
NCBI BLAST+ (v2.2.29)

Make sure these sofwares are in your PATH when running the pipeline.

Or you can install all the softwares conda package manager (https://conda.io/miniconda.html) to install required bioinformatics tools and packages in Bioconda (https://bioconda.github.io/).

conda install -c bioconda fastqc trimmomatic bwa picard gatk strelka lofreq scalpel blast

```

###### Installing

```

GOTI-seq is ready to use without any further setups in install.

cd GOTI-seq/ 
chmod +x GOTI-seq/*

```

## Usage

sh off-target-pipeline.sh work_dir sample1 sample2

###### arguments 

```
work_dir the working directorty 
sample1 the name of tdTomato+ sequencing data 
sample2 the name of tdTomato- sequencing data

```