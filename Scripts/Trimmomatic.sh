#!/bin/bash

#### 2. Quality trimming and adapter removal using Trimmomatic v0.39
mkdir /data1/WGS_Aspat_GBR/Trimmed_files
mkdir Postqfilt_quality_check/

cd /home/hugo/PhD/Genomics/Raw_data_processing/
INDIR=/data1/WGS_Aspat_GBR/230404-A00151A_L001/
OUTDIR=/data1/WGS_Aspat_GBR/Trimmed_files/

#List R1 files only 
FILES=($INDIR/*_1.fq.gz)

#Test scripts on a subset of 20 files 
FILES=("${FILES[@]:0:20}")

#Run Trimmmomatic across files 
start=`date +%s`
for FILE in ${FILES[@]}; do
	BASE=$(basename $FILE)
	BASE=${BASE%%_1*} 
	java -jar /home/gael/Programmes/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 25 -phred33 -trimlog "${OUTDIR}${BASE}_trim.log" -summary "${OUTDIR}${BASE}_sum.txt" "${INDIR}${BASE}_1.fq.gz" "${INDIR}${BASE}_2.fq.gz" "${OUTDIR}${BASE}_R1_paired.fastq.gz" "${OUTDIR}${BASE}_R1_unpaired.fastq.gz" "${OUTDIR}${BASE}_R2_paired.fastq.gz" "${OUTDIR}${BASE}_R2_unpaired.fastq.gz" ILLUMINACLIP:/home/hugo/PhD/Genomics/Raw_data_processing/Illumina_adapters_Iva_version.fa:2:30:10:2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
done
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.


# Keep bases with phred-score quality > 20 in sliding window of 4 bp (average)
# Remove adapter sequences in user specified file
# Remove reads with length < 50bp following trimming'''

#### 3. Check quality and adapter trimming

fastqc --noextract --outdir "Postqfilt_quality_check/" $(ls /data1/WGS_Aspat_GBR/Trimmed_files/*_paired*.gz) --threads 20

#Generate one general html with multiqc
''' Warning : multiqc will crash if some fastqc reports are empty (0 sequences)'''

#Include one step that deteletes empty fastq file 
multiqc Postqfilt_quality_check/ -o Postqfilt_quality_check -f 



