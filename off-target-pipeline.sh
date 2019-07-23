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

# 2. mapping

echo "`date` : start building index " >> ${LOG}/${pos}.pipleline.log
bwa index ${REFFILE}
echo "`date` : finish building index " >> ${LOG}/${pos}.pipleline.log

echo "`date` : ${pos} start mapping " >> ${LOG}/${pos}.pipleline.log
bwa mem -t 16 -M ${REFFILE} ${WORK_DIR}/raw/${pos}_R1.fastq.gz ${WORK_DIR}/raw/${pos}_R2.fastq.gz | samtools view -bS -o ${WORK_DIR}/mapping/${pos}.bam
echo "`date` : ${pos} finish mapping" >> ${LOG}/${pos}.pipleline.log

echo "`date` : ${neg} start mapping " >> ${LOG}/${neg}.pipleline.log
bwa mem -t 16 -M ${REFFILE} ${WORK_DIR}/raw/${neg}_R1.fastq.gz ${WORK_DIR}/raw/${neg}_R2.fastq.gz | samtools view -bS -o ${WORK_DIR}/mapping/${neg}.bam
echo "`date` : ${neg} finish mapping" >> ${LOG}/${neg}.pipleline.log

# 3. picard processing

mkdir ${WORK_DIR}/picard/${pos}
mkdir ${WORK_DIR}/picard/${pos}/tmp

# picard reorder
echo "`date` : ${pos} start reordering bam" >> ${LOG}/${pos}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar ReorderSam INPUT=${WORK_DIR}/mapping/${pos}.bam  OUTPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.bam REFERENCE=${REFFILE} TMP_DIR=${WORK_DIR}/picard/${pos}/tmp
echo "`date` : ${pos} finish reordering sam" >> ${LOG}/${pos}.pipleline.log

echo "`date` : ${neg} start reordering bam" >> ${LOG}/${neg}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar ReorderSam INPUT=${WORK_DIR}/mapping/${neg}.bam  OUTPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.bam REFERENCE=${REFFILE} TMP_DIR=${WORK_DIR}/picard/${neg}/tmp
echo "`date` : ${neg} finish reordering sam" >> ${LOG}/${neg}.pipleline.log

# picard sort
echo "`date` : ${pos} start sorting bam" >> ${LOG}/${pos}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar SortSam INPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.bam OUTPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.bam SORT_ORDER=coordinate TMP_DIR=${WORK_DIR}/picard/${pos}/tmp
echo "`date` : ${pos} finish sorting sam" >> ${LOG}/${pos}.pipleline.log

echo "`date` : ${neg} start sorting bam" >> ${LOG}/${neg}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar SortSam INPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.bam OUTPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.bam SORT_ORDER=coordinate TMP_DIR=${WORK_DIR}/picard/${neg}/tmp
echo "`date` : ${neg} finish sorting sam" >> ${LOG}/${neg}.pipleline.log

# picard add read groups
echo "`date` : ${pos} start assigning read group" >> ${LOG}/${pos}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.bam OUTPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.bam RGLB=WGS RGPL=Illumina RGPU=HiSeq RGSM=${pos}
echo "`date` : ${pos} finish assigning read group" >> ${LOG}/${pos}.pipleline.log

echo "`date` : ${neg} start assigning read group" >> ${LOG}/${neg}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.bam OUTPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.bam RGLB=WGS RGPL=Illumina RGPU=HiSeq RGSM=${neg}
echo "`date` : ${neg} finish assigning read group" >> ${LOG}/${neg}.pipleline.log

# picard mark duplicated reads
echo "`date` : ${pos} start marking duplicates" >> ${LOG}/${pos}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.bam OUTPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam METRICS_FILE=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.metrics TMP_DIR=${WORK_DIR}/picard/${pos}/tmp
echo "`date` : ${pos} finish marking duplicates" >> ${LOG}/${pos}.pipleline.log

echo "`date` : ${neg} start marking duplicates" >> ${LOG}/${neg}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.bam OUTPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam METRICS_FILE=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.metrics TMP_DIR=${WORK_DIR}/picard/${neg}/tmp
echo "`date` : ${neg} finish marking duplicates" >> ${LOG}/${neg}.pipleline.log

# build index for bam file
echo "`date` : ${pos} start building bam index" >> ${LOG}/${pos}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar BuildBamIndex VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam
echo "`date` : ${pos} finish building bam index" >> ${LOG}/${pos}.pipleline.log

echo "`date` : ${neg} start building bam index" >> ${LOG}/${neg}.pipleline.log
java -Xmx20g -jar ${SOFT}/picard-tools-2.3.0/picard.jar BuildBamIndex VALIDATION_STRINGENCY=SILENT INPUT=${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam
echo "`date` : ${neg} finish building bam index" >> ${LOG}/${neg}.pipleline.log

# Variant calling

# 4. Mutect variant calling

mkdir ${WORK_DIR}/mutect/${pos}

echo "`date` : ${pos} start mutect2 variant calling" >> ${LOG}/${pos}.pipleline.log
java -Xmx20g -jar ${SOFT}/GenomeAnalysisTK.jar -R ${REFFILE} -T MuTect2 -I:tumor ${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam -I:normal ${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam -o ${WORK_DIR}/mutect/${pos}/output.vcf
echo "`date` : ${pos} finish mutect2 variant calling" >> ${LOG}/${pos}.pipleline.log

# 5. Lofreq variant calling

mkdir ${WORK_DIR}/lofreq/${pos}

echo "`date` : ${pos} start lofreq variant calling" >> ${LOG}/${pos}.pipleline.log
${SOFT}/lofreq_star-2.1.2/bin/lofreq somatic -n ${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam -t ${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam -f ${REFFILE} --threads 8 -o ${WORK_DIR}/lofreq/${pos}
echo "`date` : ${pos} finish lofreq variant calling" >> ${LOG}/${pos}.pipleline.log

# 6. Strelka variant calling

mkdir ${WORK_DIR}/strelka/${pos}

echo "`date` : ${pos} start strelka variant calling" >> ${LOG}/${pos}.pipleline.log
${SOFT}/strelka-2.7.1.centos5_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam ${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam --tumorBam ${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam --referenceFasta ${REFFILE} --runDir ${WORK_DIR}/strelka/${pos}
${WORK_DIR}/strelka/${pos}/runWorkflow.py -m local -j 8
echo "`date` : ${pos} finish strelka variant calling" >> ${LOG}/${pos}.pipleline.log

# 7. Scalpel variant calling

mkdir ${WORK_DIR}/scalpel/${pos}
echo "`date` : ${pos} start scalpel variant calling" >> ${LOG}/${pos}.pipleline.log
${SOFT}/scalpel-0.5.3/scalpel-discovery --somatic --normal ${WORK_DIR}/picard/${neg}/${neg}.reorder.sort.add.mkdup.bam --tumor ${WORK_DIR}/picard/${pos}/${pos}.reorder.sort.add.mkdup.bam --bed ${BED} --window 600 --numprocs 12 --ref ${REFFILE} --dir ${WORK_DIR}/scalpel/${pos}
echo "`date` : ${pos} finish scalpel variant calling" >> ${LOG}/${pos}.pipleline.log

# 8. Filter and summarize variant calling results

awk '$7=="PASS" {print $0}' ${WORK_DIR}/mutect/${pos}/output.vcf > ${WORK_DIR}/mutect/${pos}/output.pass.vcf
awk 'length($4)==1 && length($5)==1 {print $0}' ${WORK_DIR}/mutect/${pos}/output.pass.vcf > ${WORK_DIR}/mutect/${pos}/output.pass.snv.vcf
awk 'length($4)>1 || length($5)>1 {print $0}' ${WORK_DIR}/mutect/${pos}/output.pass.vcf > ${WORK_DIR}/mutect/${pos}/output.pass.indel.vcf

awk '$7=="PASS" {print $0}' ${WORK_DIR}/strelka/${pos}/results/variants/somatic.snvs.vcf > ${WORK_DIR}/strelka/${pos}/results/variants/somatic.snvs.pass.vcf
awk '$7=="PASS" {print $0}' ${WORK_DIR}/strelka/${pos}/results/variants/somatic.indels.vcf > ${WORK_DIR}/strelka/${pos}/results/variants/somatic.indels.pass.vcf

perl ${SOFT}/filter_overlap.pl ${WORK_DIR}/mutect/${pos}/output.pass.snv.vcf ${WORK_DIR}/lofreq/${pos}/somatic_final.snvs.vcf ${WORK_DIR}/strelka/${pos}/results/variants/somatic.snvs.pass.vcf ${WORK_DIR}/${pos}.snv.overlap.vcf
perl ${SOFT}/filter_overlap.pl ${WORK_DIR}/mutect/${pos}/output.pass.indel.vcf ${WORK_DIR}/scalpel/${pos}/somatic_final.indels.vcf ${WORK_DIR}/strelka/${pos}/results/variants/somatic.indels.pass.vcf ${pos}.indel.overlap.vcf

# Filtration

${SOFT}/annotate_variation.pl -buildver mm10 -downdb -webfrom annovar refGene mousedb/
${SOFT}/table_annovar.pl  ${WORK_DIR}/${pos}.snv.overlap.vcf mousedb -buildver mm10 -out ${WORK_DIR}/${pos}.snv.overlap.vcf.anno -remove -protocol refGene -operation g -nastring . -vcfinput

perl ${SOFT}/awk_anno.pl ${WORK_DIR}/${pos}.snv.overlap.vcf.anno.mm10_multianno.txt ${WORK_DIR}/${pos}.anno.tsv
awk '$10>0.1 {print $0}' ${WORK_DIR}/${pos}.anno.tsv > ${WORK_DIR}/${pos}.anno.0.1.tsv

# Sequence comparison with on-target site

cat ${WORK_DIR}/${pos}.anno.0.1.tsv | while read line
do
	chr=$(echo $line | cut -d" " -f 1)
	pos=$(echo $line | cut -d" " -f 2)
	start=`expr $pos - 5`
	end=`expr $pos + 17`
	samtools faidx ${REFFILE} $chr:${start}-${end} >> ${WORK_DIR}/${pos}.anno.0.1.tsv.fasta
done

makeblastdb -in ${WORK_DIR}/sgRNA.fasta -dbtype nucl -parse_seqids
blastn -db ${WORK_DIR}/sgRNA.fasta -query ${WORK_DIR}/${pos}.anno.0.1.tsv.fasta -dust no -outfmt 6 -word_size=7 -out ${WORK_DIR}/${pos}.sgRNA.blast.out
