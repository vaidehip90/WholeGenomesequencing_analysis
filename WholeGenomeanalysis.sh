#!/bin/bash

STEP 1:Download the data 
echo "Downloading the data from 1000 genome project(Phase 3)!"

wget -P /home/vaidpatel575/WholeExomesequencing/data http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/sequence_read/ERR242939_1.filt.fastq.gz
wget -P /home/vaidpatel575/WholeExomesequencing/data http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/sequence_read/ERR242939_2.filt.fastq.gz
wget -P /home/vaidpatel575/WholeExomesequencing/data http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/sequence_read/ERR242943_1.filt.fastq.gz
wget -P /home/vaidpatel575/WholeExomesequencing/data http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/sequence_read/ERR242943_2.filt.fastq.gz
wget -P /home/vaidpatel575/WholeExomesequencing/data http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/sequence_read/ERR242947_1.filt.fastq.gz
wget -P /home/vaidpatel575/WholeExomesequencing/data http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/sequence_read/ERR242947_2.filt.fastq.gz

Download the reference genome
echo "Downloading the reference genome!"

wget -P /home/vaidpatel575/WholeExomesequencing/ref https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /home/vaidpatel575/WholeExomesequencing/ref/hg38.fa.gz

echo "Creating a Dictionary file for using gatk tools!"

gatk CreateSequenceDictionary R=home/vaidpatel575/WholeExomesequencing/ref/hg38.fa O=home/vaidpatel575/WholeExomesequencing/ref/hg38.dict

Directories
ref="/home/vaidpatel575/WholeExomesequencing/ref/hg38.fa"
rwdata="/home/vaidpatel575/WholeExomesequencing/data/"
alignreads="/home/vaidpatel575/WholeExomesequencing/aligned_reads"
results_all="/home/vaidpatel575/WholeExomesequencing/results"
dbsnps="/home/vaidpatel575/WholeExomesequencing/ref"
final_results="/home/vaidpatel575/WholeExomesequencing/final_result"

STEP 2: Checking the Quality of the fastq files
echo "Running fastqc"
fastqc ${rwdata}/*.fastq -o ${rwdata}/
All the quality check metrics passed so no trimming needed



STEP 3 :Alignment to reference genome
echo "Running alignment minimap2"
bwa mem index ${ref}
bwa mem -t 4 -R "@RG\tID:ERR242939\tPL:ILLUMINA\tSM:ERR242939" ${ref} ${rwdata}/ERR242939_1.filt.fastq ${rwdata}/ERR242939_2.filt.fastq > ${alignreads}/ERR242939_bwa.sam
bwa mem -t 4 -R "@RG\tID:ERR242943\tPL:ILLUMINA\tSM:ERR242943" ${ref} ${rwdata}/ERR242943_1.filt.fastq ${rwdata}/ERR242943_2.filt.fastq > ${alignreads}/ERR242943_bwa.sam
bwa mem -t 4 -R "@RG\tID:ERR242947\tPL:ILLUMINA\tSM:ERR242947" ${ref} ${rwdata}/ERR242947_1.filt.fastq ${rwdata}/ERR242947_2.filt.fastq > ${alignreads}/ERR242947_bwa.sam



STEP 4: Removing any duplicate reads using Mark duplicate
echo "Running Mark duplicate tool"
gatk MarkDuplicatesSpark -I ${alignreads}/ERR242939_bwa.sam -O ${alignreads}/ERR242939_deduped.bam
gatk MarkDuplicatesSpark -I ${alignreads}/ERR242943_bwa.sam -O ${alignreads}/ERR242943_deduped.bam
gatk MarkDuplicatesSpark -I ${alignreads}/ERR242947_bwa.sam -O ${alignreads}/ERR242947_deduped.bam



STEP 5 : Running gatk Mutect2 variant caller to call somatic variants
echo "Running mutect2 variant caller!"
gatk Mutect2 -R ${ref} -I ${alignreads}/ERR242939_deduped.bam -O ${results_all}/ERR242939_deduped_unfiltered.vcf
gatk FilterMutectCalls -R ${ref} -V ${results_all}/ERR242939_deduped_unfiltered.vcf -O ${results_all}/ERR242939_filtered.vcf

gatk Mutect2 -R ${ref} -I ${alignreads}/ERR242943_deduped.bam -O ${results_all}/ERR242943_deduped_unfiltered.vcf
gatk FilterMutectCalls -R ${ref} -V ${results_all}/ERR242943_deduped_unfiltered.vcf -O ${results_all}/ERR242943_filtered.vcf

gatk Mutect2 -R ${ref} -I ${alignreads}/ERR242947_deduped.bam -O ${results_all}/ERR242947_deduped_unfiltered.vcf
gatk FilterMutectCalls -R ${ref} -V ${results_all}/ERR242947_deduped_unfiltered.vcf -O ${results_all}/ERR242947_filtered.vcf

echo "Downloading pre-packageed data sources available in gatk Funcotator!"

Installation of pre-pacakaged data sources of gatk Funcotator
In the home directory,gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download



STEP 6: Variant Annotation using gatk Funcotator
echo "Running gatk Funcotator!"
gatk Funcotator --variant ${results_all}/ERR242939_filtered.vcf --reference ${ref} --ref-version hg38 --remove-filtered-variants --data-sources-path /home/vaidpatel575/funcotator_dataSources.v1.7.20200521s --output ${results_all}/ERR242939_annotated.vcf --output-file-format VCF

gatk Funcotator --variant ${results_all}/ERR242943_filtered.vcf --reference ${ref} --ref-version hg38 --remove-filtered-variants  --data-sources-path /home/vaidpatel575/funcotator_dataSources.v1.7.20200521s  --output ${results_all}/ERR242943_annotated_cosmic.vcf --output-file-format VCF

gatk Funcotator --variant ${results_all}/ERR242947_filtered.vcf --reference ${ref} --ref-version hg38 --remove-filtered-variants --data-sources-path /home/vaidpatel575/funcotator_dataSources.v1.7.20200521s --output ${results_all}/ERR242947_annotated.vcf --output-file-format VCF



STEP 7:VCF to table by extracting some field
echo "Using gatk VariantsToTable for extracting field from vcf to table format!"

gatk VariantsToTable -V ${results_all}/ERR242939_annotated.vcf  -F DP -F ECNT -F FUNCOTATION -O ${results_all}/ERR242939.table
gatk VariantsToTable -V ${results_all}/ERR242943_annotated.vcf -F DP -F ECNT  -F FUNCOTATION -O ${results_all}/ERR242943.table
gatk VariantsToTable -V ${results_all}/ERR242947_annotated.vcf -F DP -F ECNT  -F FUNCOTATION -O ${results_all}/ERR242947.table



STEP 8:Assigning the header to tables
echo "Adding headers to the table to check variants present in the file!"
All variants present
cat ${results_all}/ERR242939_annotated.vcf  | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${final_results}/ERR242939_curatedvariants.txt
cat ${results_all}/ERR242943_annotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${final_results}/ERR242943_curatedvariants.txt
cat ${results_all}/ERR242947_annotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${final_results}/ERR242947_curatedvariants.txt

cat ${results_all}/ERR242939.table | cut -f 3 |sed 's/|/\t/g' >> ${final_results}/ERR242939_curatedvariants.txt
cat ${results_all}/ERR242943.table | cut -f 3 |sed 's/|/\t/g' >> ${final_results}/ERR242943_curatedvariants.txt
cat ${results_all}/ERR242947.table | cut -f 3 |sed 's/|/\t/g' >> ${final_results}/ERR242947_curatedvariants.txt

Extacting somatic variants of the interest 
echo "Extacting the header to concentate into file with somatic variants!"
cat ${results_all}/ERR242939_annotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${final_results}/ERR242939_somaticvariants.txt
cat ${results_all}/ERR242943_annotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${final_results}/ERR242943_somaticvariants.txt
cat ${results_all}/ERR242947_annotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${final_results}/ERR242947_somaticvariants.txt

echo "Extacting only somatic mutations realted to malignant solid tumor!"
cat ${results_all}/ERR242939.table | cut -f 3 | grep "AKT1" |sed 's/|/\t/g' >> ${final_results}/ERR242939_somaticvariants.txt
cat ${results_all}/ERR242943.table | cut -f 3 | grep "AKT1" |sed 's/|/\t/g' >> ${final_results}/ERR242943_somaticvariants.txt
cat ${results_all}/ERR242947.table | cut -f 3 | grep "AKT1" |sed 's/|/\t/g' >> ${final_results}/ERR242947_somaticvariants.txt
