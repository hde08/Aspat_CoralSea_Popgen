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

#Trimmomatic  (version 0.39)
'''Check fastq to know quality score encoding'''
'''The thresholds used are a simplified log-likelihood approach. Each matching base adds just over 0.6, while each mismatch reduces the alignment score by Q/10. Therefore, a perfect match of a 12 base sequence will score just over 7, while 25 bases are needed to score 15. As such we recommend values between 7 - 15 for this parameter. For palindromic matches, a longer alignment is possible - therefore this threshold can be higher, in the range of 30. The seed mismatch parameter is used to make alignments more efficient, specifying the maximum base mismatch count in the seed (16 bases). Typical values here are 1 or 2'''

#Print head of fastq files to check adapter presence
'''zcat /data1/WGS_Aspat_GBR/230404-A00151A_L001/RRAP-ECT01-2022-Aspat-CBHE-1718_L1_1.fq.gz | head -n 8'''

'''lenght filtering is not needed if aligned to regerence genome with soft clipping'''

file="RRAP-ECT01-2022-Aspat-CBHE-1718_L1"
outpath="/data1/WGS_Aspat_GBR/Trimmed_files/"

#Test on most problematic files and a few good ones 
list_files=("/data1/WGS_Aspat_GBR/230404-A00151A_L001/RRAP-ECT01-2022-Aspat-HICK-481_L1" "/data1/WGS_Aspat_GBR/230404-A00151A_L002/RRAP-ECT01-2022-Aspat-CBHE-1724_L2" "/data1/WGS_Aspat_GBR/230404-A00151A_L002/RRAP-ECT01-2022-Aspat-HICK-435_L2" "/data1/WGS_Aspat_GBR/230404-A00151A_L004/RRAP-ECT01-2022-Aspat-ONMO-888_L4" "/data1/WGS_Aspat_GBR/230404-A00151A_L001/RRAP-ECT01-2022-Aspat-MART-550_L1")
l="${#list_files[@]}"

#Check adapter files that are used in trimmomatic 

for ((i=0; i<l; i++)); do
	file=${list_files[i]}
	shortfile=$(basename $file)
	java -jar /home/gael/Programmes/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 20 -phred33 -trimlog "${outpath}${shortfile}_trim.log" -summary "${outpath}${shortfile}_sum.txt" "${file}_1.fq.gz" "${file}_2.fq.gz" "${outpath}${shortfile}_R1_paired.fastq.gz" "${outpath}${shortfile}_R1_unpaired.fastq.gz" "${outpath}${shortfile}_R2_paired.fastq.gz" "${outpath}${shortfile}_R2_unpaired.fastq.gz" ILLUMINACLIP:/home/hugo/PhD/Genomics/Raw_data_processing/Illumina_adapters.fa:2:30:10:2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
done


'''start=`date +%s`
java -jar /home/gael/Programmes/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 20 -phred33 -trimlog "${outpath}${file}_trim.log" -summary "${outpath}${file}_sum.txt" "${FOLDER_PATH[i]}${file}_1.fq.gz" "${FOLDER_PATH[i]}${file}_2.fq.gz" "${outpath}${file}_R1_paired.fastq.gz" "${outpath}${file}_R1_unpaired.fastq.gz" "${outpath}${file}_R2_paired.fastq.gz" "${outpath}${file}_R2_unpaired.fastq.gz" ILLUMINACLIP:/home/hugo/PhD/Genomics/Raw_data_processing/Illumina_adapters.fa:2:30:10:2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.'''

#Compare IVA and I trimming files 
file="/data1/WGS_Aspat_GBR/230404-A00151A_L001/RRAP-ECT01-2022-Aspat-HICK-481_L1"
shortfile=$(basename $file)
outpath="/data1/WGS_Aspat_GBR/Trimmed_files/"

java -jar /home/gael/Programmes/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 20 -phred33 -trimlog "${outpath}${shortfile}_trim.log" -summary "${outpath}${shortfile}_sum.txt" "${file}_1.fq.gz" "${file}_2.fq.gz" "${outpath}${shortfile}_R1_paired.fastq.gz" "${outpath}${shortfile}_R1_unpaired.fastq.gz" "${outpath}${shortfile}_R2_paired.fastq.gz" "${outpath}${shortfile}_R2_unpaired.fastq.gz" ILLUMINACLIP:/home/hugo/PhD/Genomics/Raw_data_processing/Illumina_adapters.fa:2:30:10:2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50


java -jar /home/gael/Programmes/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 20 -phred33 -trimlog "${outpath}${shortfile}_Iva_trim.log" -summary "${outpath}${shortfile}_Iva_sum.txt" "${file}_1.fq.gz" "${file}_2.fq.gz" "${outpath}${shortfile}_Iva_R1_paired.fastq.gz" "${outpath}${shortfile}_Iva_R1_unpaired.fastq.gz" "${outpath}${shortfile}_Iva_R2_paired.fastq.gz" "${outpath}${shortfile}_Iva_R2_unpaired.fastq.gz" ILLUMINACLIP:/home/hugo/PhD/Genomics/Raw_data_processing/Illumina_adapters_Iva_version.fa:2:30:10:2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

fastqc --noextract --outdir "Postqfilt_quality_check/" "${outpath}${shortfile}_R1_paired.fastq.gz" --threads 20






#BBDuk
file="RRAP-ECT01-2022-Aspat-HICK-481_L1"
outpath="/data1/WGS_Aspat_GBR/Trimmed_bbduk/"

export PATH="/Directory1:$PATH"
 
 i=0
 #Adapter trimming
 start=`date +%s`
bbduk.sh -Xmx10g in="${FOLDER_PATH[i]}${file}_1.fq.gz" in2="${FOLDER_PATH[i]}${file}_2.fq.gz" out1="${outpath}${file}_R1_clean.fastq.gz" out2="${outpath}${file}_R2_clean.fastq.gz" ref=adapters.fa stats="${outpath}${file}_stats.txt" ktrim=r k=31 mink=15 hdist=1 tpe tbo
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.

 #Quality trimming
bbduk.sh -Xmx10g in="${outpath}${file}_R1_clean.fastq.gz" in="${outpath}${file}_R2_clean_q.fastq.gz" out1="${outpath}${file}_R1_clean.fastq.gz" out2="${outpath}${file}_R2_clean_q.fastq.gz" qtrim=rl trimq=20


#### 3. Check quality and adapter filtering

mkdir Postqfilt_quality_check/
fastqc --noextract --outdir "Postqfilt_quality_check/" $(ls /data1/WGS_Aspat_GBR/Trimmed_files/*.gz) --threads 20
	
#Generate one general html with multiqc
''' Warning : multiqc will crash if some fastqc reports are empty (0 sequences)'''
multiqc Postqfilt_quality_check/ -o Postqfilt_quality_check -f 

#### 4. Map to reference genome 

cd /data1/WGS_Aspat_GBR/
mkdir Aligned_files/
outpath="/data1/WGS_Aspat_GBR/Aligned_files/"
inpath="/data1/WGS_Aspat_GBR/Trimmed_files/"

#Bwa mem (version 0.7.17-r1188)

#1 map to A.millepora v2 reference genome 

#Index database 
bwa index /home/hugo/PhD/Genomics/Ref_genomes/Amillepora_ncbi_dataset/data/GCA_013753865.1/GCA_013753865.1_Amil_v2.1_genomic.fna
REF="/home/hugo/PhD/Genomics/Ref_genomes/Amillepora_ncbi_dataset/data/GCA_013753865.1/GCA_013753865.1_Amil_v2.1_genomic.fna"

#Map to reference genome 

bwa mem ${REF} -t 20 -o .sam "${outpath}${shortfile}_pe_aln.sam.gz" "${inpath}${shortfile}_R1_paired.fastq.gz" "${inpath}${shortfile}_R2_paired.fastq.gz"


# map paired reads to reference genome
bwa mem -t 2 ${REF} ${INDIR}/${BASE}_R1_paired.fastq.gz ${INDIR}/${BASE}_R2_paired.fastq.gz > ${OUTDIR}/${BASE}.sam

# create sorted bam
samtools view -b ${OUTDIR}/${BASE}.sam > ${OUTDIR}/${BASE}_UNSORTED.bam
samtools sort ${OUTDIR}/${BASE}_UNSORTED.bam -O BAM -o ${OUTDIR}/${BASE}_UNDEDUP.bam

# removed temp files
rm ${OUTDIR}/${BASE}_UNSORTED.bam
rm ${OUTDIR}/${BASE}.sam


outpath="/data1/WGS_Aspat_GBR/Trimmed_files/"


for ((i=0; i<l; i++)); do
	file=${list_files[i]}
	shortfile=$(basename $file)
	java -jar /home/gael/Programmes/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 20 -phred33 -trimlog "${outpath}${shortfile}_trim.log" -summary "${outpath}${shortfile}_sum.txt" "${file}_1.fq.gz" "${file}_2.fq.gz" "${outpath}${shortfile}_R1_paired.fastq.gz" "${outpath}${shortfile}_R1_unpaired.fastq.gz" "${outpath}${shortfile}_R2_paired.fastq.gz" "${outpath}${shortfile}_R2_unpaired.fastq.gz" ILLUMINACLIP:/home/hugo/PhD/Genomics/Raw_data_processing/Illumina_adapters.fa:2:30:10:2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
done


#2 map to A.spathulata contig-level reference genome 

