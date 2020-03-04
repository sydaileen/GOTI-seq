# GOTI-seq
Detect off-target variants of genome editing tools

## Install

##### Dependencies

The implement of GOTI-seq requires the following sofwares pre-installed on your systems:

[fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.11.3)
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (v0.36)
[BWA](http://bio-bwa.sourceforge.net/) (v0.7.12)
[Samtools](http://www.htslib.org/) (v1.3)
[Java SE environment](https://www.oracle.com/java/technologies/javase-jre8-downloads.html) (v1.8.0)
[Picard-tools](https://broadinstitute.github.io/picard/) (v2.3.0)
[GATK](https://gatk.broadinstitute.org/hc/en-us) (v3.5)
[Lofreq](https://csb5.github.io/lofreq/) (v2.1.2)
[Strelka](https://github.com/Illumina/strelka) (v2.7.1)
[Scalpel](http://scalpel.sourceforge.net/) (v0.5.3)
[Annovar](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/download/) (version 2016-02-01)
[NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (v2.2.29)

Make sure these softwares are in the PATH environment before the implement.

Or you can install the required softwares using [conda](https://conda.io/miniconda.html) package manager to install required bioinformatics tools and packages in [Bioconda](https://bioconda.github.io/).

conda install -c bioconda fastqc trimmomatic bwa samtools picard gatk lofreq strelka scalpel blast


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
The readers can generate files in the example_data folder by running the demo.sh script directly.

To get the results in the anticipated_results folder, please download the raw sequencing file from the following links:

Cas9-Tyr-C_pos: https://sra-download.ncbi.nlm.nih.gov/traces/sra77/SRR/008398/SRR8599743
Cas9-Tyr-C_neg: https://sra-download.ncbi.nlm.nih.gov/traces/sra77/SRR/008398/SRR8599744
BE3-Tyr-C-#1_pos: https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR8169133/SRR8169133.1
BE3-Tyr-C-#1_neg: https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR8169132/SRR8169132.1

Put these files into the "raw" directory and then run anticipated_results.sh.
