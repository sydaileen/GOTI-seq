#!/bin/bash
#$ -S /bin/bash
#$ -cwd

WORK_DIR=$1

REFFILE=${WORK_DIR}/ref/mm10.fa
BED=${WORK_DIR}/ref/mm10.bed
SOFT=${WORK_DIR}/soft
LOG=${WORK_DIR}/log

pos=$2
neg=$3

# Processing of raw reads

# 1. Quality control

# 1.1 pretrim

echo "`date` : ${pos} start fastqc checking of raw data" >> ${LOG}/${pos}.pipeline.log
${SOFT}/fastqc ${WORK_DIR}/raw/${pos}_R1.fastq.gz ${WORK_DIR}/raw/${pos}_R2.fastq.gz -o ${WORK_DIR}/fastQC/pretrim
echo "`date` : ${pos} finish fastqc checking of raw data" >> ${LOG}/${pos}.pipeline.log

echo "`date` : ${neg} start fastqc checking of raw data" >> ${LOG}/${neg}.pipeline.log
${SOFT}/fastqc ${WORK_DIR}/raw/${neg}_R1.fastq.gz ${WORK_DIR}/raw/${neg}_R2.fastq.gz -o ${WORK_DIR}/fastQC/pretrim
echo "`date` : ${neg} finish fastqc checking of raw data" >> ${LOG}/${neg}.pipeline.log

# 1.2 trim (Alternative, parameters should be adjusted based on the report of previous step)

echo "`date` : ${pos} start Trimmomatic trimming" >> ${LOG}/${pos}.pipeline.log
java -jar ${SOFT}/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 ${WORK_DIR}/raw/${pos}_R1.fastq.gz ${WORK_DIR}/raw/${pos}_R2.fastq.gz \
	${WORK_DIR}/fastQC/trim/${pos}_1_paired.fastq.gz ${WORK_DIR}/fastQC/trim/${pos}_1_unpaired.fastq.gz \
	${WORK_DIR}/fastQC/trim/${pos}_2_paired.fastq.gz ${WORK_DIR}/fastQC/trim/${pos}_2_unpaired.fastq.gz \
	ILLUMINACLIP:${SOFT}/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 MINLEN:70
echo "`date` : ${pos} finish Trimmomatic trimming" >> ${LOG}/${pos}.pipeline.log

echo "`date` : ${neg} start Trimmomatic trimming" >> ${LOG}/${neg}.pipeline.log
java -jar ${SOFT}/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 ${WORK_DIR}/raw/${neg}_R1.fastq.gz ${WORK_DIR}/raw/${neg}_R2.fastq.gz \
	${WORK_DIR}/fastQC/trim/${neg}_1_paired.fastq.gz ${WORK_DIR}/fastQC/trim/${neg}_1_unpaired.fastq.gz \
	${WORK_DIR}/fastQC/trim/${neg}_2_paired.fastq.gz ${WORK_DIR}/fastQC/trim/${neg}_2_unpaired.fastq.gz \
	ILLUMINACLIP:${SOFT}/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 MINLEN:70
echo "`date` : ${neg} finish Trimmomatic trimming" >> ${LOG}/${neg}.pipeline.log

# 1.3 posttrim (Alternative)

echo "`date` : ${pos} start fastqc checking of clean data" >> ${LOG}/${pos}.pipeline.log
${SOFT}/fastqc ${WORK_DIR}/fastQC/trim/${pos}_1_paired.fastq.gz ${WORK_DIR}/fastQC/trim/${pos}_2_paired.fastq.gz -o ${WORK_DIR}/fastQC/posttrim
echo "`date` : ${pos} finish fastqc checking of clean data" >> ${LOG}/${pos}.pipeline.log

echo "`date` : ${neg} start fastqc checking of clean data" >> ${LOG}/${neg}.pipeline.log
${SOFT}/fastqc ${WORK_DIR}/fastQC/trim/${neg}_1_paired.fastq.gz ${WORK_DIR}/fastQC/trim/${neg}_2_paired.fastq.gz -o ${WORK_DIR}/fastQC/posttrim
echo "`date` : ${neg} finish fastqc checking of clean data" >> ${LOG}/${neg}.pipeline.log

# 2. mapping

echo "`date` : start building index " >> ${LOG}/${pos}.pipeline.log
${SOFT}/bwa index ${REFFILE}
echo "`date` : finish building index " >> ${LOG}/${pos}.pipeline.log


echo "`date` : ${pos} start mapping " >> ${LOG}/${pos}.pipeline.log
${SOFT}/bwa mem -t 16 -M ${REFFILE} ${WORK_DIR}/raw/${pos}_R1.fastq.gz ${WORK_DIR}/raw/${pos}_R2.fastq.gz | ${SOFT}/samtools view -bS -o ${WORK_DIR}/mapping/${pos}.bam
echo "`date` : ${pos} finish mapping" >> ${LOG}/${pos}.pipeline.log

echo "`date` : ${neg} start mapping " >> ${LOG}/${neg}.pipeline.log
${SOFT}/bwa mem -t 16 -M ${REFFILE} ${WORK_DIR}/raw/${neg}_R1.fastq.gz ${WORK_DIR}/raw/${neg}_R2.fastq.gz | ${SOFT}/samtools view -bS -o ${WORK_DIR}/mapping/${neg}.bam
echo "`date` : ${neg} finish mapping" >> ${LOG}/${neg}.pipeline.log

source ~/.bashrc

# 3. picard processing

mkdir ${WORK_DIR}/picard/${pos}
mkdir ${WORK_DIR}/picard/${pos}/tmp

# picard reorder
echo "`date` : ${pos} start reordering bam" >> ${LOG}/${pos}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar ReorderSam INPUT=${WORK_DIR}/mapping/${pos}.bam  OUTPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.bam REFERENCE=${REFFILE} TMP_DIR=${WORK_DIR}/picard/${pos}/tmp
echo "`date` : ${pos} finish reordering sam" >> ${LOG}/${pos}.pipeline.log

echo "`date` : ${neg} start reordering bam" >> ${LOG}/${neg}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar ReorderSam INPUT=${WORK_DIR}/mapping/${neg}.bam  OUTPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.bam REFERENCE=${REFFILE} TMP_DIR=${WORK_DIR}/picard/${neg}/tmp
echo "`date` : ${neg} finish reordering sam" >> ${LOG}/${neg}.pipeline.log

# picard sort
echo "`date` : ${pos} start sorting bam" >> ${LOG}/${pos}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar SortSam INPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.bam OUTPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.bam SORT_ORDER=coordinate TMP_DIR=${WORK_DIR}/picard/${pos}/tmp
echo "`date` : ${pos} finish sorting sam" >> ${LOG}/${pos}.pipeline.log

echo "`date` : ${neg} start sorting bam" >> ${LOG}/${neg}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar SortSam INPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.bam OUTPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.bam SORT_ORDER=coordinate TMP_DIR=${WORK_DIR}/picard/${neg}/tmp
echo "`date` : ${neg} finish sorting sam" >> ${LOG}/${neg}.pipeline.log

# picard add read groups
echo "`date` : ${pos} start assigning read group" >> ${LOG}/${pos}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.bam OUTPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.bam RGLB=WGS RGPL=Illumina RGPU=HiSeq RGSM=${pos}
echo "`date` : ${pos} finish assigning read group" >> ${LOG}/${pos}.pipeline.log

echo "`date` : ${neg} start assigning read group" >> ${LOG}/${neg}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.bam OUTPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.bam RGLB=WGS RGPL=Illumina RGPU=HiSeq RGSM=${neg}
echo "`date` : ${neg} finish assigning read group" >> ${LOG}/${neg}.pipeline.log

# picard mark duplicated reads
echo "`date` : ${pos} start marking duplicates" >> ${LOG}/${pos}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.bam OUTPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam METRICS_FILE=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.metrics TMP_DIR=${WORK_DIR}/picard/${pos}/tmp
echo "`date` : ${pos} finish marking duplicates" >> ${LOG}/${pos}.pipeline.log

echo "`date` : ${neg} start marking duplicates" >> ${LOG}/${neg}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.bam OUTPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam METRICS_FILE=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.metrics TMP_DIR=${WORK_DIR}/picard/${neg}/tmp
echo "`date` : ${neg} finish marking duplicates" >> ${LOG}/${neg}.pipeline.log

# build index for bam file
echo "`date` : ${pos} start building bam index" >> ${LOG}/${pos}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar BuildBamIndex VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam
echo "`date` : ${pos} finish building bam index" >> ${LOG}/${pos}.pipeline.log

echo "`date` : ${neg} start building bam index" >> ${LOG}/${neg}.pipeline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar BuildBamIndex VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam
echo "`date` : ${neg} finish building bam index" >> ${LOG}/${neg}.pipeline.log

# Variant calling

# 4. Mutect variant calling

mkdir ${WORK_DIR}/mutect/${pos}

echo "`date` : ${pos} start mutect2 variant calling" >> ${LOG}/${pos}.pipeline.log
java -Xmx20g -jar ${SOFT}/GenomeAnalysisTK.jar -R ${REFFILE} -T MuTect2 -I:tumor ${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam -I:normal ${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam -o ${WORK_DIR}/mutect/${pos}/output.vcf
echo "`date` : ${pos} finish mutect2 variant calling" >> ${LOG}/${pos}.pipeline.log

# 5. Lofreq variant calling

mkdir ${WORK_DIR}/lofreq/${pos}

echo "`date` : ${pos} start lofreq variant calling" >> ${LOG}/${pos}.pipeline.log
${SOFT}/lofreq_star-2.1.2/bin/lofreq somatic -n ${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam -t ${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam -f ${REFFILE} --threads 8 -o ${WORK_DIR}/lofreq/${pos}/
echo "`date` : ${pos} finish lofreq variant calling" >> ${LOG}/${pos}.pipeline.log

# 6. Strelka variant calling

mkdir ${WORK_DIR}/strelka/${pos}

echo "`date` : ${pos} start strelka variant calling" >> ${LOG}/${pos}.pipeline.log
${SOFT}/strelka-2.7.1.centos5_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam ${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam --tumorBam ${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam --referenceFasta ${REFFILE} --runDir ${WORK_DIR}/strelka/${pos}
${WORK_DIR}/strelka/${pos}/runWorkflow.py -m local -j 8
echo "`date` : ${pos} finish strelka variant calling" >> ${LOG}/${pos}.pipeline.log


# 7. Scalpel variant calling

mkdir ${WORK_DIR}/scalpel/${pos}
echo "`date` : ${pos} start scalpel variant calling" >> ${LOG}/${pos}.pipeline.log
${SOFT}/scalpel-0.5.3/scalpel-discovery --somatic --normal ${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam --tumor ${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam --bed ${BED} --window 600 --numprocs 12 --ref ${REFFILE} --dir ${WORK_DIR}/scalpel/${pos}
echo "`date` : ${pos} finish scalpel variant calling" >> ${LOG}/${pos}.pipeline.log


# 8. Filter and summarize variant calling results

awk '$7=="PASS" {print $0}' ${WORK_DIR}/mutect/${pos}/output.vcf > ${WORK_DIR}/mutect/${pos}/output.pass.vcf
awk 'length($4)==1 && length($5)==1 {print $0}' ${WORK_DIR}/mutect/${pos}/output.pass.vcf > ${WORK_DIR}/mutect/${pos}/output.pass.snv.vcf
awk 'length($4)>1 || length($5)>1 {print $0}' ${WORK_DIR}/mutect/${pos}/output.pass.vcf > ${WORK_DIR}/mutect/${pos}/output.pass.indel.vcf

gunzip ${WORK_DIR}/strelka/${pos}/results/variants/somatic.snvs.vcf.gz
gunzip ${WORK_DIR}/strelka/${pos}/results/variants/somatic.indels.vcf.gz

awk '$7=="PASS" {print $0}' ${WORK_DIR}/strelka/${pos}/results/variants/somatic.snvs.vcf > ${WORK_DIR}/strelka/${pos}/results/variants/somatic.snvs.pass.vcf
awk '$7=="PASS" {print $0}' ${WORK_DIR}/strelka/${pos}/results/variants/somatic.indels.vcf > ${WORK_DIR}/strelka/${pos}/results/variants/somatic.indels.pass.vcf

gunzip ${WORK_DIR}/lofreq/${pos}/somatic_final.snvs.vcf.gz

perl filter_overlap.pl ${WORK_DIR}/lofreq/${pos}/somatic_final.snvs.vcf ${WORK_DIR}/strelka/${pos}/results/variants/somatic.snvs.pass.vcf ${WORK_DIR}/mutect/${pos}/output.pass.snv.vcf ${WORK_DIR}/${pos}.snv.overlap.vcf
perl filter_overlap.pl ${WORK_DIR}/scalpel/${pos}/somatic_final.indels.vcf ${WORK_DIR}/strelka/${pos}/results/variants/somatic.indels.pass.vcf ${WORK_DIR}/mutect/${pos}/output.pass.indel.vcf ${pos}.indel.overlap.vcf

# Filtration

# ${SOFT}/annotate_variation.pl -buildver mm10 -downdb -webfrom annovar refGene ${WORK_DIR}/ref/mousedb/
${SOFT}/annovar/table_annovar.pl  ${WORK_DIR}/${pos}.snv.overlap.vcf ${WORK_DIR}/ref/mousedb -buildver mm10 -out ${WORK_DIR}/${pos}.snv.overlap.vcf.anno -remove -protocol refGene -operation g -nastring . -vcfinput

perl awk_anno.pl ${WORK_DIR}/${pos}.snv.overlap.vcf.anno.mm10_multianno.txt ${WORK_DIR}/${pos}.anno.tsv
awk '$10>0.1 {print $0}' ${WORK_DIR}/${pos}.anno.tsv > ${WORK_DIR}/${pos}.anno.0.1.tsv

# Sequence comparison with on-target site

cat ${WORK_DIR}/${pos}.anno.0.1.tsv | while read line
do
	chr=$(echo $line | cut -d" " -f 1)
	loci=$(echo $line | cut -d" " -f 2)
	start=`expr $loci - 17`
	end=`expr $loci + 5`
	${SOFT}/samtools faidx ${REFFILE} $chr:${start}-${end} >> ${WORK_DIR}/${pos}.anno.0.1.tsv.fasta
done

${SOFT}/makeblastdb -in ${WORK_DIR}/sgRNA.fasta -dbtype nucl -parse_seqids
${SOFT}/blastn -db ${WORK_DIR}/sgRNA.fasta -query ${WORK_DIR}/${pos}.anno.0.1.tsv.fasta -dust no -outfmt 6 -word_size=7 -out ${WORK_DIR}/${pos}.sgRNA.blast.out
