#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Trimmomatic_2nd_round_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J Trimmomatic

### WALLTIME
#SBATCH -t 56:00:00

#MEMORY
#SBATCH --mem=10G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)

#Output and error directory
#SBATCH -o /home/hdenis/Slurm/trim_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/trim_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=1		# number of cores per job

#SBATCH --array=1-350%20        	# job array

### Ressources

# This job's working directory :
echo Working Directory is $SLURM_SUBMIT_DIR
export PBS_O_WORKDIR=$SLURM_SUBMIT_DIR

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo NTASKS $SLURM_NTASKS

echo This job runs on the following processors :
##PBS_NODEFILE=`srun hostname`
PBS_NODEFILE=$SLURM_SUBMIT_DIR/node_list_$SLURM_JOB_ID.txt
printf '%s\n' `srun hostname` > $PBS_NODEFILE
export PBS_NODEFILE=$PBS_NODEFILE

###echo `cat $SLURM_JOB_NODELIST`
#`scontrol show hostnames $SLURM_JOB_NODELIST`

NPROCS=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`
echo This job has allocated $NPROCS cpus
echo This job has allocated $NNODES servers

# This job's working directory :
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

ulimit -s unlimited

########################################################## READS TRIMMING 2nd round ###################################################################
#This scripts serves to remove adapters and trim low quality bases from raw reads 

#1. Running Fastqc on trimmed files detected rare occurence of sequencing artifacs with sequences containing only A or T
#We Perform trimmomatic again using those sequences as adapter removal 
#TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

#mkdir /data1/WGS_Aspat_GBR/Trimmed_files
mkdir /nvme/disk0/lecellier_data/WGS_NC_data/Postqfilt_quality_check/QC_2nd_round/
mkdir /nvme/disk0/lecellier_data/WGS_NC_data/Trimmed_files/Trim_2nd_round/

cd /nvme/disk0/lecellier_data/WGS_NC_data/
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Trimmed_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Trimmed_files/Trim_2nd_round/"

#List R1 files only 
FILES=($INDIR/*_R1_paired.fastq.gz)

#Associate filename with array index 
FILENAME=${FILES[$((${SLURM_ARRAY_TASK_ID}-1))]}
BASE=$(basename $FILENAME)
BASE=${BASE%%_*}

#If statement to avoid reprocessing file
if [ ! -s "${OUTDIR}{BASE}_R1_paired.fastq.gz" ] 
then
  
  start=`date +%s`
  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} : start trimming 
  
  #Run Trimmomatic 
  #Trim log removed as it takes too much storage can be added with -trimlog "${OUTDIR}{}_trim.log"
  start=`date +%s`
  java -jar /home/hdenis/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 -summary "${OUTDIR}${BASE}_sum.txt" "${INDIR}${BASE}_R1_paired.fastq.gz" "${INDIR}${BASE}_R2_paired.fastq.gz" "${OUTDIR}${BASE}_R1_paired.fastq.gz" "${OUTDIR}${BASE}_R1_unpaired.fastq.gz" "${OUTDIR}${BASE}_R2_paired.fastq.gz" "${OUTDIR}${BASE}_R2_unpaired.fastq.gz" ILLUMINACLIP:/nvme/disk0/lecellier_data/WGS_NC_data/artifacts.fa:2:30:10:4:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
  end=`date +%s`
  echo Execution time was `expr $(( ($end - $start) / 60))` minutes.
  
  #Run fastq on trimmed files 
  fastqc --noextract --outdir "Postqfilt_quality_check/QC_2nd_round/" "${OUTDIR}${BASE}_R1_paired.fastq.gz" --threads 1
  fastqc --noextract --outdir "Postqfilt_quality_check/QC_2nd_round/" "${OUTDIR}${BASE}_R2_paired.fastq.gz" --threads 1
  
else

  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} already processed

fi


#2. Run multiqc on all samples (2n trimmed round)
multiqc "/nvme/disk0/lecellier_data/WGS_NC_data/Postqfilt_quality_check/QC_2nd_round/" -o "/nvme/disk0/lecellier_data/WGS_NC_data/Postqfilt_quality_check/QC_2nd_round/" -f -d


