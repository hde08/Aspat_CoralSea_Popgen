#!/bin/bash

#Load conda environment with graftM package
#Moved to the Submit_nci_parallel.sh script 

pid=$(grep ^Pid /proc/self/status)
corelist=$(grep Cpus_allowed_list: /proc/self/status | awk '{print $2}')
host=$(hostname | sed 's/.gadi.nci.org.au//g')
echo subtask $1 running in $pid using cores $corelist on compute node $host

cd /scratch/d85/hd9321/
INDIR="/scratch/d85/hd9321/Unmapped_files/Unmapped_Fastq/"
OUTDIR="/scratch/d85/hd9321/psbAncr/"

SEED="/scratch/d85/hd9321/0_symb_ref/2_markers/Cladocopium_psbA_frontier.fa"

#Read number of PBS task and associate to sample ID 
FILES=($INDIR*_1.fq.gz)

PSBANCR_FILES="${OUTDIR}psbAncr_gbr_nc_samples.list.txt"

TASK_ID=$1
NAME=${FILES[$((${TASK_ID}-1))]}
BASE=$(basename $NAME)
BASE=${BASE%%_1*}

#Index reference file 
#bwa index $SEED

#1. Map reads against psbAncr seed
#Using bwa mem v0.7.17
#if [ ! -s "${OUTDIR}${FILE}.filtered100-90.bam" ] 
#  then
#  
#  #Only process files that are in the psbAncr list
#  if cat $PSBANCR_FILES | grep -q $BASE; then 
#    
#    start=`date +%s`
#    echo Array Id : ${TASK_ID} Sample : ${BASE} : start processing
#    
#    #Map reads against psbAncr seed
#    bwa mem $SEED "${INDIR}${BASE}_1.fq.gz" "${INDIR}${BASE}_2.fq.gz" | samtools view -b -@ 4 -F 4 - -o "${OUTDIR}${BASE}.bam"
#    
#    # keep only full-length reads aligned with more than 90% of identity
#    perl /home/600/hd9321/Scripts/BamFiltration.pl -in "${OUTDIR}${BASE}.bam" -out "${OUTDIR}${BASE}.filtered100-90.bam" -minPCaligned 100 -minPCidentity 90    
#    
#    end=`date +%s`
#    echo  Execution time was `expr $(( ($end - $start) / 60))` minutes.
#  
#  fi
#   
#  else
#
#  echo Sample : ${BASE} already processed
#fi
  
  
##2. Identify reference with highest number of reads mapped 
#for FILE in ${OUTDIR}*.filtered100-90.bam
#do
#    echo -ne "$FILE\t" ;samtools view $FILE | awk '{print $3}' | sort | uniq -c | sort -nrk1,1 |head -1 | awk '{print $2"\t"$1}'
#done > ${OUTDIR}Summary_mapping_Cladocopium_psbA_100-90.tab

 #echo -ne "$FILE\t" ;samtools view $FILE | awk '{print $3}' | sort | uniq -c | sort -nrk1,1 | awk '{print $2"\t"$1}'
#done > ${OUTDIR}Summary_mapping_Cladocopium_psbA_100-90_allref.tab

##2.2 Count uniquely mapped reads as proxy of abundance
for FILE in ${OUTDIR}*.filtered100-90.bam
do
    echo -ne "$FILE\t" ;samtools view $FILE | awk '{print $3}' | sort | uniq -c | sort -nrk1,1 |head -1 | awk '{print $2"\t"$1}'
done > ${OUTDIR}Summary_mapping_Cladocopium_psbA_100-90.tab

 echo -ne "$FILE\t" ;samtools view $FILE | awk '{print $3}' | sort | uniq -c | sort -nrk1,1 | awk '{print $2"\t"$1}'
done > ${OUTDIR}Summary_mapping_Cladocopium_psbA_100-90_allref.tab

#-F 4: Excludes unmapped reads.
#-F 256: Excludes secondary alignments.
#-F 2048: Excludes supplementary alignments.


#3. Sort bam files 
#cat ${OUTDIR}Summary_mapping_Cladocopium_psbA_100-90.tab | while read a b c ;do samtools sort -o ${a%.filtered100-90.bam}.sorted.bam $a; done
#
##4. Index sorted bam files 
#for i in ${OUTDIR}*.sorted.bam
#do
#    samtools index $i
#done
##
##
##5. Generate text pileup output
#cat ${OUTDIR}Summary_mapping_Cladocopium_psbA_100-90.tab | while read a b c
#do
#    samtools mpileup -Aa -f $SEED -r $b ${a%.filtered100-90.bam}.sorted.bam -o ${a%.filtered100-90.bam}.sorted.mpileup
#done
#
##6. building consensus (minimum coverage of 2) 
#for i in ${OUTDIR}*sorted.mpileup
#do
#    perl /home/600/hd9321/Scripts/MpileuptoConsensus.pl -mpileup $i -Mincov 2 -out ${i}.consensus.fa
#done
#
#7; keep consensus sequences in one fasta file 
#for i in ${OUTDIR}*.sorted.mpileup.consensus.fa
#do
#    base=$(basename $i)
#    echo -e ">${base%_*}"; tail -n +2 $i
#done > ${OUTDIR}All_psbA-Cladocopium100-90_psbA.aln.mpileup.consensus.fa


######### Create a phylogenetic tree with additional psbAncr reference sequences ########

#Create joint fasta file with consensus sequences and ref
#cat ${OUTDIR}All_psbA-Cladocopium100-90_psbA.aln.mpileup.consensus.fa "/scratch/d85/hd9321/0_symb_ref/2_markers/psbA_ref_marker_sequences.fa" > ${OUTDIR}psbA_consensus_and_ref.fa
#
##Create Multi sequence alignment using mafft v7.525
#/home/600/hd9321/bin/mafft --maxiterate 1000 --localpair --thread 24 ${OUTDIR}psbA_consensus_and_ref.fa > ${OUTDIR}psbA_consensus_and_ref_alignment.fa
##
### convert fasta to nexus
#/home/600/hd9321/.local/bin/seqmagick convert --output-format nexus --alphabet dna ${OUTDIR}psbA_consensus_and_ref_alignment.fa ${OUTDIR}psbA_consensus_and_ref_alignment.nexus
#
##Replace quotation marks that have been introduced in nexus file so that MrBayes can run 
#sed -e "s/'//g" ${OUTDIR}psbA_consensus_and_ref_alignment.nexus > ${OUTDIR}psbA_consensus_and_ref_alignment_formated.nexus

#Create Bayesian Tree using Mr Bayes 
cd /scratch/d85/hd9321/psbAncr/MrBayes/
/home/600/hd9321/bin/bin/mb /home/600/hd9321/Scripts/MrbayesCommands.nex


