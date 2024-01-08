#!/bin/bash

#To make the file executable run the command below in terminal 
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/Quality_check_raw_fasta.sh '''

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
