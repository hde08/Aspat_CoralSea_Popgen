#!/bin/bash

pid=$(grep ^Pid /proc/self/status)
corelist=$(grep Cpus_allowed_list: /proc/self/status | awk '{print $2}')
host=$(hostname | sed 's/.gadi.nci.org.au//g')
echo subtask $1 running in $pid using cores $corelist on compute node $host

########################################################## PREPROCESS Symbiont Reads ###################################################################
#Script to process non coral reads to separate Symbiodiniacea reads from other DNA sources (bacteria, archae, viruses...)

#Load samtools v1.19 & bedtools v2.31.0
module load samtools
module load bedtools 



#1. There seems to have some symbiont contamination (scaffold 382) in A.millepora ref genomes thus we will extract reads mapping to this scaffold and add them to the non coral reads, prior to run any other analyses 
#This is done separately for NC and GBR samples
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Unmapped_fastq/"
BASENAME_FILES=($(find "$INDIR" -maxdepth 1 -name "*_1.fq.gz"))

INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Unmapped_fastq/"
BASENAME_FILES=($(find "$INDIR" -maxdepth 1 -name "*_1.fq.gz"))

for FILE in "${BASENAME_FILES[@]}"; do
#Find location of marked dup files between GBR and NC batch 
    BASE=$(basename $FILE)
    BASE=${BASE%%_GBR*} 
    BASE=${BASE%%_NC*} 
    
    if  echo ${BASE} | grep -q '_L'; then
      DIR_PATH="/tgt_bck2/data_lecellier/Aligned_files/"
      echo $BASE in $DIR_PATH Processing 
#Extract reads mapping to scaffold_382 in ref genome and add them to the unmapped files
      samtools view -b "${DIR_PATH}${BASE}_MARKED_DUP.bam" "scaffold_382" > "${DIR_PATH}${BASE}.scaffold382.bam"
#Concatenate the two bam files
      samtools merge -o "${DIR_PATH}${BASE}.unmapped_scaffold382.bam" "${DIR_PATH}${BASE}.unmapped.bam" "${DIR_PATH}${BASE}.scaffold382.bam" 
    else
      DIR_PATH="/tgt_bck2/data_lecellier/WGS_NC_data/Aligned_files/"
      echo $BASE in $DIR_PATH Already processed 
    fi
done 

#2. Edit bam file basenames to include the batch ID (to separate technical replicates in the different batches)
FILES=(/scratch/d85/hd9321/Unmapped_files/Unmapped_GBR/*unmapped_scaffold382.bam)
for FILE in ${FILES[@]}; do
    BASE=$(basename $FILE)
    BASE=${BASE%%.un*}
    NEW_BASE="${BASE}_GBR_batch.unmapped_scaffold382.bam"
    mv ${FILE} "/scratch/d85/hd9321/Unmapped_files/Unmapped_GBR/${NEW_BASE}"
done
#
FILES=(/scratch/d85/hd9321/Unmapped_files/Unmapped_NC/*unmapped_scaffold382.bam)
for FILE in ${FILES[@]}; do
    BASE=$(basename $FILE)
    BASE=${BASE%%.un*}
    NEW_BASE="${BASE}_NC_batch.unmapped_scaffold382.bam"
    mv ${FILE} "/scratch/d85/hd9321/Unmapped_files/Unmapped_NC/${NEW_BASE}"
done

#3. Convert bam files back to fasq format 

#GBR samples 
OUTDIR="/scratch/d85/hd9321/Unmapped_files/Unmapped_Fastq/"
INDIR="/scratch/d85/hd9321/Unmapped_files/Unmapped_GBR/"
readarray -t BAM_FILES <  $INDIR/gbr_bam.filelist.txt
for FILE in ${BAM_FILES[@]};do
  BASE=$(basename $FILE)
  BASE=${BASE%%.unmapped*}
  bedtools bamtofastq -i "${INDIR}${BASE}.unmapped_scaffold382.bam" -fq "${OUTDIR}${BASE}_1.fq" -fq2 "${OUTDIR}${BASE}_2.fq"
  
done

gzip ${OUTDIR}*.fq

#NC samples 
OUTDIR="/scratch/d85/hd9321/Unmapped_files/Unmapped_Fastq/"
INDIR="/scratch/d85/hd9321/Unmapped_files/Unmapped_NC/"
readarray -t BAM_FILES <  $INDIR/nc_bam.filelist.txt
for FILE in ${BAM_FILES[@]:261:55};do
  BASE=$(basename $FILE)
  BASE=${BASE%%.unmapped*}
  bedtools bamtofastq -i "${INDIR}${BASE}.unmapped_scaffold382.bam" -fq "${OUTDIR}${BASE}_1.fq" -fq2 "${OUTDIR}${BASE}_2.fq"
  
done

gzip ${OUTDIR}*.fq


#3. Align non coral reads against a custom database of symbiont reference genomes (one per genus)

#Genomes used for the alignment are :
#-Symbiodinium microadriaticum CCMP2467 (Nand et al., 2021)
#-Breviolum: Breviolum minutum Mf1.05b (Shoguchi et al., 2013). Re-annotated in (Chen et al., 2022)
#-Cladocopium: Cladocopium proliferum SCF055-01 (Chen et al., 2022)
#-Durusdinium: Durusdinium trenchii CCMP2556 (Dougan et al ., 2024)
#-Effrenium: Effrenium voratum RCC1521 (Shah et al., 2023)
#-Fugacium: Fugacium kawagutii CCMP2468 (Li et al., 2019)

#Concatenate reference genomes in single file 
touch "${REFDIR}Symb_ref_genome.fa"
REF_GENOMES=($REFDIR*.fa.gz)
for REF in ${REF_GENOMES[@]};
  do
  zcat $REF >> "${REFDIR}Symb_ref_genome.fa"
done
  
#Index fasta file with reference genomes 
bwa index "${REFDIR}Symb_ref_genome.fa"

#List of non coral reads file 
FASTQ_FILES=($(find "$INDIR" -maxdepth 1 -name "*_1.fq.gz"))

#Match slurm array ID to file 
FILE=${FASTQ_FILES[$((${SLURM_ARRAY_TASK_ID}-1))]}
#FILE=${FASTQ_FILES[$((${SLURM_ARRAY_TASK_ID}+400-1))]}
BASE=$(basename $FILE)
BASE=${BASE%%_1*}

#Align reads, extract mapped reads only and convert back to fasta format 
if [ ! -s "${OUTDIR}${BASE}.fa" ] 
then
  
  start=`date +%s`
  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} : start processing
  
  bwa mem "${REFDIR}Symb_ref_genome.fa" "${INDIR}${BASE}_1.fq.gz" "${INDIR}${BASE}_2.fq.gz" -o "${OUTDIR}${BASE}.sam"    
  samtools view "${OUTDIR}${BASE}.sam" -F 4 -b | samtools sort > "${OUTDIR}${BASE}.bam"
  samtools fasta "${OUTDIR}${BASE}.bam" > "${OUTDIR}${BASE}.fa"
  
  #Remove intermediate files
  rm "${OUTDIR}${BASE}.sam"
  rm "${OUTDIR}${BASE}.bam"
 
  end=`date +%s`
  echo ${BASE} : Execution time was `expr $(( ($end - $start) / 60))` minutes.
   
else

  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} already processed

fi




