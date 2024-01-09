#!/bin/bash
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/picard.sh '''

#5. Mark and delete PCR duplicates using picard version 2.27.5 

#Recommanded to generate before after quality plot to check base recalibration 
cd /data1/WGS_Aspat_GBR
mkdir /home/hugo/PhD/Genomics/Raw_data_processing/BAM_statistics 

INDIR="/data1/WGS_Aspat_GBR/Aligned_files/"
OUTDIR="/home/hugo/PhD/Genomics/Raw_data_processing/BAM_statistics/"

FILES=($INDIR*_UNDEDUP.bam)

#Create file IDs
>/home/hugo/PhD/Genomics/Raw_data_processing/ids.txt
for FILE in ${FILES[@]}; do
	BASE=$(basename $FILE)
	BASE=${BASE%%_U*} 
	echo ${BASE} >> /home/hugo/PhD/Genomics/Raw_data_processing/ids.txt
done 

#5.1 Add read group information (parallel mode)
#Besoin de régler le problème des echo | cut 
start=`date +%s`
cat /home/hugo/PhD/Genomics/Raw_data_processing/ids.txt | parallel --jobs 25 "echo /home/gael/Programmes/gatk-4.3.0.0/gatk AddOrReplaceReadGroups INPUT="${INDIR}{}_UNDEDUP.bam" OUTPUT="${INDIR}{}_UNDEDUP_RG.bam" RGID=$(echo {} | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/") RGLB=$(echo {} | cut -d"-" -f1,2,3) RGPL=ILLUMINA RGPU=$(echo {} | cut -d"_" -f2) RGSM=$(echo {} | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/")"
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.


#Single file mode 
/home/gael/Programmes/gatk-4.3.0.0/gatk AddOrReplaceReadGroups INPUT="${INDIR}${BASE}_UNDEDUP.bam" OUTPUT="${INDIR}${BASE}_UNDEDUP_RG.bam" RGID=$(echo $BASE | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/") RGLB=$(echo $BASE | cut -d"-" -f1,2,3) RGPL=ILLUMINA RGPU=$(echo $BASE | cut -d"_" -f2) RGSM=$(echo $BASE | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/")

'''View read group : samtools view -H "${INDIR}${BASE}_UNDEDUP.bam" | grep '^@RG''''


#RGID : unique read group identifier
#RGLB : library index
#RGPL : platform name
#RGPU : flowcell, lane, sample barcode 
#RGSM : sample name (=RGID in this case)

#5.2 Mark and remove duplicates 

/home/gael/Programmes/gatk-4.3.0.0/gatk MarkDuplicates INPUT="${INDIR}${BASE}_UNDEDUP_RG.bam" OUTPUT="${INDIR}${BASE}_MARKED_DUP.bam" METRICS_FILE="${OUTDIR}${BASE}_marked_dup_metrics.txt" REMOVE_DUPLICATES=true TMP_DIR="${INDIR}" VALIDATION_STRINGENCY=LENIENT


      
