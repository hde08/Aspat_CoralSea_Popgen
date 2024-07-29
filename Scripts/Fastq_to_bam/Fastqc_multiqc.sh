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

#Generate multiple multiqc after running Trimmomatic (for visualization purpose)
#List fastqc files 
mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq1/
mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq2/
mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq3/
mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq4/
mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq5/

FASTQC=(/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/*_fastqc.zip)

mv -t /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq1/ "${FASTQC[@]:0:300}"
mv -t /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq2/ "${FASTQC[@]:300:300}"
mv -t /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq3/ "${FASTQC[@]:600:300}"
mv -t /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq4/ "${FASTQC[@]:900:300}"
mv -t /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq5/ "${FASTQC[@]:1200:334}"

#Run multiqc in each folder 
multiqc "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq1/" -o "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq1/" -f -d
multiqc "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq2/" -o "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq2/" -f -d
multiqc "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq3/" -o "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq3/" -f -d
multiqc "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq4/" -o "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq4/" -f -d
multiqc "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq5/" -o "/nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/Fastq5/" -f -d