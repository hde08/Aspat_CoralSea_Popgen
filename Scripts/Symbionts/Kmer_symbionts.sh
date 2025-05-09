#!/bin/bash

#Load conda environment with graftM package
#Moved to the Submit_nci_parallel.sh script 

########################################################## METAGENOMIC DISTANCE ###################################################################
#Script to count kmer from Symbiodiniaceae reads found in each sample and used to compute a metagenomic dissimilarity matrix (D2S distance)

pid=$(grep ^Pid /proc/self/status)
corelist=$(grep Cpus_allowed_list: /proc/self/status | awk '{print $2}')
host=$(hostname | sed 's/.gadi.nci.org.au//g')
echo subtask $1 running in $pid using cores $corelist on compute node $host

cd /scratch/d85/hd9321/
INDIR="/scratch/d85/hd9321/Kmer_files/"
OUTDIR="/scratch/d85/hd9321/Kmer_files/"

#1. Count kmer in each samples 
#As samples contain uneven number of reads, each sample is subsampled without replacement to 50Mbp (corresponding to the 0.1 quantile of the distribution across samples). 
#Samples originally containing <50Mbp are then excluded from the analysis 
source activate base
conda activate d2ssect

FA_FILES=(${OUTDIR}*.fa)

TASK_ID=0

START=$((${TASK_ID}*100))
END=$((${TASK_ID}+100))

for ((i = $START; i < $END; i++)); do 

  FILE=${FA_FILES[i]}
  BASE=$(basename $FILE)
  BASE=${BASE%%_100M*}

  if [ ! -s "${OUTDIR}${BASE}_50M.jf" ] 
  then
    
    start=`date +%s`
    echo  File : ${BASE} : start processing
    
    #Downsample the file 
    #samplebasestarget=100M Extract 100M bp per sample
    #supsample=f prevents duplicating reads when the sample size is lower than 100M bp 
    /home/600/hd9321/Programs/bbmap/reformat.sh in="${OUTDIR}${BASE}.fa" out="${OUTDIR}${BASE}_50M.fa" samplebasestarget=50M upsample=f overwrite=t
    
    #Count kmer in sequencing reads using jellyfish
    jellyfish count -m 21 -s 500M "${OUTDIR}${BASE}_50M.fa" -o "${OUTDIR}${BASE}_50M.jf"
    
    #Remove suffix automatically outputed by jellyfish
    mv "${OUTDIR}${BASE}_50M.jf_0" "${OUTDIR}${BASE}_50M.jf"
    
    end=`date +%s`
    echo ${BASE} : Execution time was `expr $(( ($end - $start) / 60))` minutes.
     
  else
  
    echo File : ${BASE} already processed
  
  fi


done

#2. Compute D2S distance 
#This can be done for all samples at once. 
#However due to the large number of samples including in this analysis, the original script was modified to compute individual pairwise distances separately. 
#Each pairwise distance is then aggregated into a single matrix using a custom R script 

#Create file with non-directional pairs of files 
JF_FILES=(${OUTDIR}*50M.jf)
FA_FILES=(${OUTDIR}*50M.fa)

mkdir "${OUTDIR}Dist_sample_pairs/"
touch "${OUTDIR}_JF_sample_pairs.txt"
touch "${OUTDIR}_FA_sample_pairs.txt"
for ((i = 0; i < ${#JF_FILES[@]}; i++)); do 
      for ((j = i + 1; j < ${#JF_FILES[@]}; j++)); do 
        printf "%s\t%s\n" "${JF_FILES[i]}" "${JF_FILES[j]}" >> "${OUTDIR}_JF_sample_pairs.txt"
        printf "%s\t%s\n" "${FA_FILES[i]}" "${FA_FILES[j]}" >> "${OUTDIR}_FA_sample_pairs.txt"
      done;
done 

#Match the subtask index $1 to file 
TASK_ID=$1

for ((i = 0; i < 750; i++)); do 

  LINE_ID="$((1+(TASK_ID-1)*750+i))"
  
  IFS=$'\t' JF_PAIR=($(awk "NR==$LINE_ID" "${OUTDIR}_JF_sample_pairs.txt"))
  IFS=$'\t' FA_PAIR=($(awk "NR==$LINE_ID" "${OUTDIR}_FA_sample_pairs.txt"))
    
  NAME1=$(basename ${JF_PAIR[0]})
  NAME1=${NAME1%%_pe*}
  NAME1=${NAME1%%_L*}
    
  NAME2=$(basename ${JF_PAIR[1]})
  NAME2=${NAME2%%_pe*}
  NAME2=${NAME2%%_L*}

  if [ ! -s "${OUTDIR}Dist_sample_pairs/${NAME1}_${NAME2}_50M_pairwise_dist.txt" ] 
  then
    
    
    start=`date +%s`
    echo Array Id : ${TASK_ID} Line Id : ${LINE_ID} Pair ${NAME1} ${NAME2} : start processing
    
    #Use d2ssect to compute the overall distance matrix between the pair 
    d2ssect -l ${JF_PAIR[@]} -f ${FA_PAIR[@]} -o "${OUTDIR}Dist_sample_pairs/${NAME1}_${NAME2}_50M_pairwise_dist.txt" -t 1
    
    end=`date +%s`
    echo  Execution time was `expr $(( ($end - $start) / 60))` minutes.
  else

    echo Pair ${NAME1} ${NAME2} already processed
  fi
        
done;