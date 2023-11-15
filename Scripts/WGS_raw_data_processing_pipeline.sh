#!/bin/bash

#To make the file executable run the command below in terminal 
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/WGS_raw_data_processing_pipeline.sh '''

#To unzip folder 
'''unzip -o /data1/WGS_Aspat_GBR/230404-A00199A_L001.zip'''

#To check file copy integrity
'''md5sum -c /data1/WGS_Aspat_GBR/230404-A00151A_L001/md5.md5'''

#Set path of raw sequencing data (fastq files)
FOLDER_PATH=($(ls -d  /data1/WGS_Aspat_GBR/*/))

#Change working dirctory and create directory to store initial quality checks 
cd /home/hugo/PhD/Genomics/Raw_data_processing/
mkdir Initial_quality_check/

#### 1.Run fastqc on each individual fastq file ###

l="${#FOLDER_PATH[@]}"

for ((i=4; i<l; i++)); do
	#Run fastqc over all fastq files  
	fastqc --noextract --outdir "Initial_quality_check/" $(ls ${FOLDER_PATH[i]}*.gz) --threads 20
	
done

#Generate one general html with multiqc
multiqc Initial_quality_check/ -o Initial_quality_check -f 



#### 2. Quality trimming and adapter removal 

#Try both trimmomatic and bbduk on a subset of samples 

#Trimmomatic 
'''Check fastq to know quality score encoding'''
'''The thresholds used are a simplified log-likelihood approach. Each matching base adds just over 0.6, while each mismatch reduces the alignment score by Q/10. Therefore, a perfect match of a 12 base sequence will score just over 7, while 25 bases are needed to score 15. As such we recommend values between 7 - 15 for this parameter. For palindromic matches, a longer alignment is possible - therefore this threshold can be higher, in the range of 30. The seed mismatch parameter is used to make alignments more efficient, specifying the maximum base mismatch count in the seed (16 bases). Typical values here are 1 or 2'''

#Print head of fastq files to check adapter presence
'''zcat /data1/WGS_Aspat_GBR/230404-A00151A_L001/RRAP-ECT01-2022-Aspat-CBHE-1718_L1_1.fq.gz | head -n 8'''

file=
java -jar /home/gael/Programmes/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 20 -phred33 -trimlog "$file/Trimmomatic.log" -summary "$file/Trimmomatic_summary.txt" "$file/_1" "$file/_2" "R1_paired.fastq.gz" "R1_unpaired.fastq.gz" "R2_paired.fastq.gz" "R2_unpaired.fastq.gz" ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100





