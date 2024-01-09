#!/bin/bash
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/BAM_statistics.sh '''

#7. Generate BAM file statistics using samtools version 1.10 

cd /data1/WGS_Aspat_GBR

INDIR="/data1/WGS_Aspat_GBR/Aligned_files/"
OUTDIR="/home/hugo/PhD/Genomics/Raw_data_processing/BAM_statistics/"

FILES=($INDIR*UNDEDUP_RG.bam)

#Loop across files 
for FILE in "${FILES[@]}"; do
    BASE=$(basename $FILE)
    BASE=${BASE%%_U*} 
    echo ${BASE} >> ${OUTDIR}/sample_name.txt		# generate sample names
    samtools view -c --threads 25 "${INDIR}${BASE}_MARKED_DUP.bam" >> ${OUTDIR}allCounts.txt		# generate read counts
    samtools view -b -f 12 -F 256 --threads 25 "${INDIR}${BASE}_MARKED_DUP.bam" > "${INDIR}${BASE}.unmapped.bam"		# generate bam with unmapped reads = symbiont, microbes and other reads 
    samtools view -c --threads 25 "${INDIR}${BASE}.unmapped.bam" >> "${OUTDIR}unmappedCounts.txt"		# generate unmapped read counts
	
done


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


