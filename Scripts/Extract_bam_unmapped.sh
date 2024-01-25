#!/bin/bash
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/Extract_bam_unmapped.sh '''

#6. Extract unmapped reads (non host reads) from bam files to separate files 

cd /data1/WGS_Aspat_GBR

INDIR="/data1/WGS_Aspat_GBR/Aligned_files/"
OUTDIR="/home/hugo/PhD/Genomics/Raw_data_processing/BAM_statistics/"

#List files for which duplicates have been previoulsy removed 
FILES=($INDIR*MARKED_DUP.bam)

#Loop across files 

start=`date +%s`
for FILE in "${FILES[@]}"; do
    BASE=$(basename $FILE)
    BASE=${BASE%%_M*} 
    # generate sample names
    echo ${BASE} >> ${OUTDIR}/sample_name.txt		
    # generate read counts
    samtools view -c --threads 25 "${INDIR}${BASE}_MARKED_DUP.bam" >> ${OUTDIR}allCounts.txt		
    # generate bam with unmapped reads = symbiont, microbes and other reads 
    samtools view -b -f 12 -F 256 --threads 25 "${INDIR}${BASE}_MARKED_DUP.bam" > "${INDIR}${BASE}.unmapped.bam"
    # generate unmapped read counts
    samtools view -c --threads 25 "${INDIR}${BASE}.unmapped.bam" >> "${OUTDIR}unmappedCounts.txt"		
done
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

# Notes:
# -f 4: extract all unmapped reads
# -f12: extract only the reads where read 1 is unmapped AND read 2 is unmapped (= both mates are unmapped). Applies only to paired reads.


