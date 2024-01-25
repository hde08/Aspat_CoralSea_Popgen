#!/bin/bash
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/BAM_statistics.sh '''

#7. Generate BAM file statistics using samtools version 1.10 

cd /data1/WGS_Aspat_GBR

INDIR="/data1/WGS_Aspat_GBR/Aligned_files/"
OUTDIR="/home/hugo/PhD/Genomics/Raw_data_processing/BAM_statistics/"

FILES=($INDIR*MARKED_DUP.bam)

#Loop across files 
start=`date +%s`
for FILE in "${FILES[@]}"; do
    BASE=$(basename $FILE)
    BASE=${BASE%%_M*} 
    #7.1 Generate index bai and statistics from BAM aligned reads 
    samtools index -@ 25 "${INDIR}${BASE}_MARKED_DUP.bam"
    #7.2 Output statistics 
    samtools idxstats --threads 25 "${INDIR}${BASE}_MARKED_DUP.bam" > "${OUTDIR}${BASE}-index_stats.txt"
    #Columns in output file are : reference sequence name, sequence lenght, mapped reads segment, unmapped reads segment 
    samtools coverage "${INDIR}${BASE}_MARKED_DUP.bam" > "${OUTDIR}${BASE}-coverage.txt"
    samtools flagstat --threads 25 -O tsv "${INDIR}${BASE}_MARKED_DUP.bam" > "${OUTDIR}${BASE}-flagstats.txt"
done
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.





