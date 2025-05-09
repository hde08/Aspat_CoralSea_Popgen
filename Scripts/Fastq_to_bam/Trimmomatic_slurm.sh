#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Trimmomatic_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J Trimmomatic

### WALLTIME
#SBATCH -t 56:00:00

#MEMORY
#SBATCH --mem=150G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)

#Output and error directory
#SBATCH -o /home/hdenis/Slurm/trim_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/trim_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=10		# number of cores per job

#SBATCH --array=1%1        	# job array

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

########################################################## READS TRIMMING ###################################################################
#This scripts serves to remove adapters and trim low quality bases from raw reads 

####1. Quality trimming and adapter removal using Trimmomatic v0.39
#mkdir /data1/WGS_Aspat_GBR/Trimmed_files
mkdir /nvme/disk0/lecellier_data/WGS_NC_data/Postqfilt_quality_check/
mkdir /nvme/disk0/lecellier_data/WGS_NC_data/Trimmed_files/

cd /nvme/disk0/lecellier_data/WGS_NC_data/
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/HN00216654_hdd1/RawData/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Trimmed_files/"

#List R1 files only 
FILES=($INDIR/*_1.fastq.gz)

#Associate filename with array index 
FILENAME=${FILES[$((${SLURM_ARRAY_TASK_ID}-1))]}
BASE=$(basename $FILENAME)
BASE=${BASE%%_*}

#Trimming parameters
#Keep bases with phred-score quality > 20 in sliding window of 4 bp (average)
#Keep reads with a minimum length of 50bp after trimming
#Remove leading and trailing low quality bases (3)
#Remove adapters using the illuminaclip option in 'palindrome mode'

#If statement to avoid reprocessing file
if [ ! -s "${OUTDIR}{BASE}_R1_paired.fastq.gz" ] 
then
  
  start=`date +%s`
  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} : start trimming 
  
  #Run Trimmomatic 
  #Trim log removed as it takes too much storage can be added with -trimlog "${OUTDIR}{}_trim.log"
  start=`date +%s`
  java -jar /home/hdenis/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 -summary "${OUTDIR}${BASE}_sum.txt" "${INDIR}${BASE}_1.fastq.gz" "${INDIR}${BASE}_2.fastq.gz" "${OUTDIR}${BASE}_R1_paired.fastq.gz" "${OUTDIR}${BASE}_R1_unpaired.fastq.gz" "${OUTDIR}${BASE}_R2_paired.fastq.gz" "${OUTDIR}${BASE}_R2_unpaired.fastq.gz" ILLUMINACLIP:/nvme/disk0/lecellier_data/WGS_NC_data/Illumina_adapters_Iva_version.fa:2:30:10:4:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
  end=`date +%s`
  echo Execution time was `expr $(( ($end - $start) / 60))` minutes.
  
  #Run fastq on trimmed files to confirm adapters have been removed
  fastqc --noextract --outdir "Postqfilt_quality_check/" "${OUTDIR}${BASE}_R1_paired.fastq.gz" --threads 1
  fastqc --noextract --outdir "Postqfilt_quality_check/" "${OUTDIR}${BASE}_R2_paired.fastq.gz" --threads 1
  
else

  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} already processed

fi


#2. Identify empty files after trimmomatic and delete them to avoid crashing multiqc
Store their ID in a file to record files that have been eliminated
TRIMMED_FILES=(/nvme/disk0/lecellier_data/WGS_NC_data/Trimmed_files/*_paired*.gz)
for FILE in ${TRIMMED_FILES[@]}; do
    NREADS=$(awk '{s++}END{print s/4}' $FILE)
    if [ "$NREADS" = 0 ]; then
        echo $(basename ${FILE}) >> /nvme/disk0/lecellier_data/WGS_NC_data/Trimmed_files/Trimming_empty_ids.txt
        #rm $FILE
    fi   
done

#3. Run multiqc on all samples  
multiqc "/nvme/disk0/lecellier_data/WGS_NC_data/Postqfilt_quality_check/" -o "/nvme/disk0/lecellier_data/WGS_NC_data/Postqfilt_quality_check/" -f -d
#''' Warning : multiqc will crash if some fastqc reports are empty (0 sequences)'''






