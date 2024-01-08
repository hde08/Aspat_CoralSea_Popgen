#!/bin/bash
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/BAM_statistics.sh '''

#7. Generate BAM file statistics using samtools version 1.10 

cd /data1/WGS_Aspat_GBR

INDIR="/data1/WGS_Aspat_GBR/Aligned_files/"
OUTDIR="/home/hugo/PhD/Genomics/Raw_data_processing/BAM_statistics/"

FILES=($(ls -1 ${INDIR}*UNDEDUP_RG.bam))
BASE=$(basename ${FILES[0]})
BASE=${BASE%%_U*} 

#7.1 Generate index bai and statistics from BAM aligned reads 
samtools index "${INDIR}${BASE}_UNDEDUP_RG.bam"

cd ${INDIR}

#7.2 Output statistics 

samtools idxstats "${INDIR}${BASE}_UNDEDUP_RG.bam" > "${OUTDIR}${BASE}-index_stats.txt"
#Columns in output file are : reference sequence name, sequence lenght, mapped reads segment, unmapped reads segment 


samtools coverage "${INDIR}${BASE}_UNDEDUP_RG.bam" > "${OUTDIR}${BASE}-coverage.txt"
samtools flagstat -O tsv "${INDIR}${BASE}_UNDEDUP_RG.bam" > "${OUTDIR}${BASE}-flagstats.txt"


