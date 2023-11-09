#!/bin/bash

#To make the file executable run the command below in terminal 
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/WGS_raw_data_processing_pipeline.sh '''

#Set path of raw sequencing data (fastq files)
FOLDER_PATH="/data1/WGS_Aspat_GBR/230404-A00151A_L001"


#Change working dirctory and create directory to store initial quality checks 
cd /home/hugo/PhD/Genomics/Raw_data_processing/
mkdir Initial_quality_check/

#### 1.Run fastqc on each individual fastq file ###

#Run fastqc over all fastq files  
fastqc --noextract --outdir "Initial_quality_check/" $(ls $FOLDER_PATH/*.gz) --threads 15

#Generate one general html with multiqc
multiqc Initial_quality_check/ -o Initial_quality_check


#### 2. 
