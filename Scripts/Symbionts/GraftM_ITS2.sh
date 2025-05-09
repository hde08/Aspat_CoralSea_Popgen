#!/bin/bash

#Load conda environment with graftM package
#Moved to the Submit_nci_parallel.sh script 

pid=$(grep ^Pid /proc/self/status)
corelist=$(grep Cpus_allowed_list: /proc/self/status | awk '{print $2}')
host=$(hostname | sed 's/.gadi.nci.org.au//g')
echo subtask $1 running in $pid using cores $corelist on compute node $host

#GraftM creation with custom ITS2 HMM profile was previously done by Take Ishida
#see https://github.com/hisatakeishida/Symb-SHIN/blob/main/C_Read-based_alignment.md 
#graftM create --output ITS2_graftm.gpkg --sequences symportal_ITS2.fa --taxonomy ITS2_taxonomy.txt

#ITS2_taxonomy.txt looks like
#C1	k__Dinophyceae,p__Suessiales,c__Symbiodiniaceae,o__Cladocopium,f__C1,g__C1
#C1b	k__Dinophyceae,p__Suessiales,c__Symbiodiniaceae,o__Cladocopium,f__C1,g__C1b
#C1au	k__Dinophyceae,p__Suessiales,c__Symbiodiniaceae,o__Cladocopium,f__C1,g__C1au

#Apply graftM to all unmapped fastq files (recover ITS2 sequences)
INDIR="/scratch/d85/hd9321/Unmapped_files/Unmapped_Fastq/"
OUTDIR="/scratch/d85/hd9321/ITS2_GraftM/"
GRAFTM_PACK_DIR="/home/600/hd9321/Symb-SHIN/"

#mkdir ${OUTDIR}

#Create file with list of samples 
#FILES=($INDIR*_1.fq.gz)
#printf "%s\n" "${FILES[@]}" > $INDIR/fq.filelist.txt

readarray -t FQ_FILES <  $INDIR/fq.filelist.txt

#Match the subtask index $1 to file 
FILE=${FQ_FILES[$(($1-1))]}
BASE=$(basename $FILE)
BASE=${BASE%%_1.fq.gz}

#Apply graftM package for ITS2
if [ ! -s "${OUTDIR}${BASE}" ] 
then
  
  start=`date +%s`
  echo Array Id : ${1} File : ${BASE} : start processing
  
  /home/600/hd9321/Programs/graftM/bin/graftM graft --forward "${INDIR}${BASE}_1.fq.gz" --reverse "${INDIR}${BASE}_2.fq.gz" --graftm_package "${GRAFTM_PACK_DIR}ITS2_graftm_final.gpkg" --input_sequence_type nucleotide --output_directory "${OUTDIR}${BASE}" --threads 4

  end=`date +%s`
  echo ${BASE} : Execution time was `expr $(( ($end - $start) / 60))` minutes.
   
else

  echo Array Id : ${1} File : ${BASE} already processed

fi

#Copy each combined_count_table.txt file into a common directory 
COUNT_TAB_FILES=($(find "$OUTDIR" -maxdepth 2 -name "*combined_count_table.txt"))
COUNT_TAB_DIR="/scratch/d85/hd9321/ITS2_GraftM/Count_tabs/"
mkdir $COUNT_TAB_DIR
for FILE in ${COUNT_TAB_FILES[@]}; do
    FILENAME="${FILE%/combined*}"
    FILENAME=$(basename $FILENAME)
    mv $FILE "${COUNT_TAB_DIR}${FILENAME}.txt"
done

#Run at the end to get a summary across all samples 
python /home/600/hd9321/Scripts/graftM_result_summary.py




