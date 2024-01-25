#!/bin/bash
'''chmod u+x /home/hugo/PhD/Genomics/Scripts/Bwa_mem.sh '''

#### 4. Map to reference genome using Bwa mem (version 0.7.17-r1188) with default parameters

cd /data1/WGS_Aspat_GBR/
#mkdir Aligned_files/
INDIR="/data1/WGS_Aspat_GBR/Trimmed_files/"
OUTDIR="/data1/WGS_Aspat_GBR/Aligned_files/"

#List reference genomes to be tested 
REF_1="/home/hugo/PhD/Genomics/Ref_genomes/Amillepora_ncbi_dataset/data/GCA_013753865.1/GCA_013753865.1_Amil_v2.1_genomic.fna"
REF_2="/home/hugo/PhD/Genomics/Ref_genomes/Aspathulata_ncbi_dataset/data/GCA_031770025.1/GCA_031770025.1_AGI_CSIRO_Aspa_v1_genomic.fna"
REFS=(${REF_1} ${REF_2})

#### 4.1. Index reference genomes 
'''bwa index $REF_1 #A.millepora'''

'''bwa index $REF_2 #A.spathulata'''

#List files 
FILES=($INDIR/*_R1_paired.fastq.gz)
FILES=("${FILES[@]:0:20}")


start=`date +%s`
#Loop across reference genomes that we want to align to 
for REF in ${REFS[@]}; do
    REF_NAME=${REF%%_ncbi*}
    REF_NAME=${REF_NAME#*genomes/}
    #Loop across input files 
    for FILE in "${FILES[@]}"; do
        BASE=$(basename $FILE)
        BASE=${BASE%%_R*} 
	#### 4.2. Map to reference genome 
	bwa mem -t 25 ${REF} -o "${OUTDIR}${BASE}_pe_aln_${REF_NAME}.sam" "${INDIR}${BASE}_R1_paired.fastq.gz" "${INDIR}${BASE}_R2_paired.fastq.gz"
	#### 4.3. Create sorted bam  
	samtools view -b --threads 25 "${OUTDIR}${BASE}_pe_aln_${REF_NAME}.sam" > "${OUTDIR}${BASE}_pe_aln_${REF_NAME}_UNSORTED.bam"
	samtools sort --threads 25 "${OUTDIR}${BASE}_pe_aln_${REF_NAME}_UNSORTED.bam" -O BAM -o "${OUTDIR}${BASE}_pe_aln_${REF_NAME}_UNDEDUP.bam"
	##### 4.4. removed temporary files
	rm "${OUTDIR}${BASE}_pe_aln_${REF_NAME}_UNSORTED.bam"
	rm "${OUTDIR}${BASE}_pe_aln_${REF_NAME}.sam"	
    done
done
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.
