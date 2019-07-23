# sgRNA offtarget analysis pipeline
# Usage: sh final_somatic.sh sample

dbSNP=/picb/bigdata/project/CPDR/metagenome/reference/mm10/mgp.v3.snps.rsIDdbSNPv137.format.vcf
dbINDEL=/picb/bigdata/project/CPDR/metagenome/reference/mm10/mgp.v3.indels.rsIDdbSNPv137.format.vcf
REPEAT=/picb/bigdata/project/CPDR/metagenome/sgRNA/lib/Repeats/mm10.Repeats.all.bed
species=mm10

sample=$1

# strelka
gunzip strelka/${sample}/results/variants/somatic.indels.vcf.gz
gunzip strelka/${sample}/results/variants/somatic.snvs.vcf.gz

awk '$7=="PASS" {print $0}' strelka/${sample}/results/variants/somatic.indels.vcf > strelka/${sample}/results/variants/somatic.indels.pass.vcf
awk '$7=="PASS" {print $0}' strelka/${sample}/results/variants/somatic.snvs.vcf > strelka/${sample}/results/variants/somatic.snvs.pass.vcf

# lofreq
gunzip lofreq/${sample}/out_somatic_final.snvs.vcf.gz

# mutect
cat mutect/${sample}*/output.vcf > mutect/${sample}.output.vcf
awk '$7=="PASS" {print $0}' mutect/${sample}.output.vcf > mutect/${sample}.pass.vcf
awk 'length($4)==1 && length($5)==1 {print $0}' mutect/${sample}.pass.vcf > mutect/${sample}.pass.snv.vcf
awk 'length($4)>1 || length($5)>1 {print $0}' mutect/${sample}.pass.vcf > mutect/${sample}.pass.indel.vcf

# scalpel
cat scalpel/${sample}*/main/somatic.indel.vcf > scalpel/${sample}.somatic.indel.vcf
awk '$7=="PASS" {print $0}' scalpel/${sample}.somatic.indel.vcf > scalpel/${sample}.somatic.indel.pass.vcf

# overlap between three methods

# snv
perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_overlap.pl strelka/${sample}/results/variants/somatic.snvs.pass.vcf lofreq/${sample}/out_somatic_final.snvs.vcf mutect/${sample}.pass.snv.vcf result/${sample}.snv.overlap.vcf
perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_overlap1.pl lofreq/${sample}/out_somatic_final.snvs.vcf mutect/${sample}.pass.snv.vcf result/${sample}.snv.overlap12.vcf
perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_overlap1.pl strelka/${sample}/results/variants/somatic.snvs.pass.vcf mutect/${sample}.pass.snv.vcf result/${sample}.snv.overlap13.vcf
perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_overlap1.pl strelka/${sample}/results/variants/somatic.snvs.pass.vcf lofreq/${sample}/out_somatic_final.snvs.vcf result/${sample}.snv.overlap23.vcf

# indel
perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_overlap.pl strelka/${sample}/results/variants/somatic.indels.pass.vcf scalpel/${sample}.somatic.indel.pass.vcf mutect/${sample}.pass.indel.vcf result/${sample}.indel.overlap.vcf
perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_overlap1.pl scalpel/${sample}.somatic.indel.pass.vcf mutect/${sample}.pass.indel.vcf result/${sample}.indel.overlap12.vcf
perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_overlap1.pl strelka/${sample}/results/variants/somatic.indels.pass.vcf mutect/${sample}.pass.indel.vcf result/${sample}.indel.overlap13.vcf
perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_overlap1.pl strelka/${sample}/results/variants/somatic.indels.pass.vcf scalpel/${sample}.somatic.indel.pass.vcf result/${sample}.indel.overlap23.vcf

# filter variants in known database
# perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_known.pl ${dbINDEL} result/${sample}.indel.overlap.vcf result/${sample}.indel.overlap.fil.vcf
# perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_known.pl ${dbSNP} result/${sample}.snv.overlap.vcf result/${sample}.snv.overlap.fil.vcf

# statistics
# snv
sed -n '$=' mutect/${sample}.pass.snv.vcf >> result/${sample}.snv.stat
sed -n '$=' lofreq/${sample}/out_somatic_final.snvs.vcf >> result/${sample}.snv.stat
sed -n '$=' strelka/${sample}/results/variants/somatic.snvs.pass.vcf >> result/${sample}.snv.stat
sed -n '$=' result/${sample}.snv.overlap.vcf >> result/${sample}.snv.stat
sed -n '$=' result/${sample}.snv.overlap12.vcf >> result/${sample}.snv.stat
sed -n '$=' result/${sample}.snv.overlap13.vcf >> result/${sample}.snv.stat
sed -n '$=' result/${sample}.snv.overlap23.vcf >> result/${sample}.snv.stat
# sed -n '$=' result/${sample}.snv.overlap.fil.vcf >> result/${sample}.snv.stat
# indel
sed -n '$=' mutect/${sample}.pass.indel.vcf >> result/${sample}.indel.stat
sed -n '$=' scalpel/${sample}.somatic.indel.pass.vcf >> result/${sample}.indel.stat
sed -n '$=' strelka/${sample}/results/variants/somatic.indels.pass.vcf >> result/${sample}.indel.stat
sed -n '$=' result/${sample}.indel.overlap.vcf >> result/${sample}.indel.stat
sed -n '$=' result/${sample}.indel.overlap12.vcf >> result/${sample}.indel.stat
sed -n '$=' result/${sample}.indel.overlap13.vcf >> result/${sample}.indel.stat
sed -n '$=' result/${sample}.indel.overlap23.vcf >> result/${sample}.indel.stat
# sed -n '$=' result/${sample}.indel.overlap.fil.vcf >> result/${sample}.indel.stat

:<<BLOCK
# filter repeats
qsub /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_repeat.sh /picb/bigdata/project/CPDR/metagenome/sgRNA/lib/Repeats/Human.Repeats.all.bed result/${sample}.snv.overlap.vcf result/${sample}.snv.overlap.norep.vcf
qsub /picb/bigdata/project/CPDR/metagenome/sgRNA/script/filter_repeat.sh /picb/bigdata/project/CPDR/metagenome/sgRNA/lib/Repeats/Human.Repeats.all.bed result/${sample}.indel.overlap.vcf result/${sample}.indel.overlap.norep.vcf
sed -n '$=' result/${sample}.snv.overlap.norep.vcf >> result/${sample}.snv.stat
sed -n '$=' result/${sample}.indel.overlap.norep.vcf >> result/${sample}.indel.stat

# annotation
sed 's/ /\t/g' result/${sample}.snv.overlap.norep.vcf > result/${sample}.snv.overlap.norep.vcf1
mv result/${sample}.snv.overlap.norep.vcf1 result/${sample}.snv.overlap.norep.vcf
sh /picb/bigdata/project/CPDR/metagenome/sgRNA/script/anno.sh result/${sample}.snv.overlap.norep.vcf
awk '$0~/exonic/ {print $0}' result/${sample}.snv.overlap.norep.vcf.anno.hg19_multianno.txt > result/${sample}.snv.overlap.exonic.vcf
awk '$0~/nonsynonymous\ SNV/ {print $0}' result/${sample}.snv.overlap.exonic.vcf > result/${sample}.snv.overlap.exonic.nonsyn.vcf
sed -n '$=' result/${sample}.snv.overlap.exonic.vcf >> result/${sample}.snv.stat
sed -n '$=' result/${sample}.snv.overlap.exonic.nonsyn.vcf >> result/${sample}.snv.stat

sed 's/ /\t/g' result/${sample}.indel.overlap.norep.vcf > result/${sample}.indel.overlap.norep.vcf1
mv result/${sample}.indel.overlap.norep.vcf1 result/${sample}.indel.overlap.norep.vcf
sh /picb/bigdata/project/CPDR/metagenome/sgRNA/script/anno.sh result/${sample}.indel.overlap.norep.vcf
awk '$0~/exonic/ {print $0}' result/${sample}.indel.overlap.norep.vcf.anno.hg19_multianno.txt > result/${sample}.indel.overlap.exonic.vcf
awk '$0~/frameshift/ {print $0}' result/${sample}.indel.overlap.exonic.vcf > result/${sample}.indel.overlap.exonic.nonsyn.vcf
sed -n '$=' result/${sample}.indel.overlap.exonic.vcf >> result/${sample}.indel.stat
sed -n '$=' result/${sample}.indel.overlap.exonic.nonsyn.vcf >> result/${sample}.indel.stat

# overlap between samples
# perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/overlap.pl result/${sample1}.snv.overlap.exonic.nonsyn.vcf result/${sample2}.snv.overlap.exonic.nonsyn.vcf result/snv.overlap.exonic.nonsyn.vcf
# perl /picb/bigdata/project/CPDR/metagenome/sgRNA/script/overlap.pl result/${sample1}.indel.overlap.exonic.nonsyn.vcf result/${sample2}.indel.overlap.exonic.nonsyn.vcf result/indel.overlap.exonic.nonsyn.vcf
BLOCK
