# GOTI-seq
Detect off-target variants of genome editing tools

## Install

##### Dependencies


```

The implement of GOTI-seq requires the following sofwares pre-installed on your systems:

fastQC (v0.11.3)
Trimmomatic (v0.36)
BWA (v0.7.12)
Samtools (v1.3)
Java SE environment(v1.8.0)
Picard-tools (v2.3.0)
GATK (v3.5)
Lofreq (v2.1.2)
Strelka (v2.7.1)
Scalpel (v0.5.3)
Annovar (version 2016-02-01)
NCBI BLAST+ (v2.2.29)

Make sure these softwares are in the PATH environment before the implement.

Or you can install the required softwares using conda package manager [https://conda.io/miniconda.html] to install required bioinformatics tools and packages in [Bioconda](https://bioconda.github.io/).

conda install -c bioconda fastqc trimmomatic bwa samtools picard gatk lofreq strelka scalpel blast

```

##### Installing

```
GOTI-seq is ready to use without any further setups in install.

cd GOTI-seq
chmod +x *

```

## Usage

```
Before running the script, make sure the raw data are located in the "raw" sub-directory of working directory.

mkdir fastQC mapping picard mutect lofreq strelka scalpel

sh off-target-pipeline.sh work_dir sample1 sample2

work_dir the working directory
sample1 the name of tdTomato+ sequencing data
sample2 the name of tdTomato- sequencing data
```
