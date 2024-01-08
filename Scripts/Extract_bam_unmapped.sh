#!/bin/bash
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/Extract_bam_unmapped.sh '''

#6. Extract unmapped reads (non host reads) from bam files to separate files 

cd /data1/WGS_Aspat_GBR

INDIR="/data1/WGS_Aspat_GBR/Aligned_files/"
OUTDIR="/home/hugo/PhD/Genomics/Raw_data_processing/BAM_statistics/"

FILES=($(ls -1 ${INDIR}*MARKED_DUP.bam))
BASE=$(basename ${FILES[0]})
BASE=${BASE%%_M*} 

echo ${BASE} >> ${OUTDIR}/sample_name.txt		# generate sample names
samtools view -c "${INDIR}${BASE}_MARKED_DUP.bam" >> ${OUTDIR}allCounts.txt		# generate read counts

samtools view -b -f 12 -F 256 "${INDIR}${BASE}_MARKED_DUP.bam" > "${INDIR}${BASE}.unmapped.bam"		# generate bam with unmapped reads = symbiont, microbes and other reads 

samtools view -c "${INDIR}${BASE}.unmapped.bam" >> "${OUTDIR}unmappedCounts.txt"		# generate unmapped read counts

# Notes:
# -f 4: extract all unmapped reads
# -f12: extract only the reads where read 1 is unmapped AND read 2 is unmapped (= both mates are unmapped). Applies only to paired reads.
