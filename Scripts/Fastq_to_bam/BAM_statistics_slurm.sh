#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/BAM_statistics_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J BAM_stat

### WALLTIME
#SBATCH -t 56:00:00

#MEMORY
#SBATCH --mem=15G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)

#Output and error directory
#SBATCH -o /home/hdenis/Slurm/bam_stat_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/bam_stat_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=2		# number of cores per job

#SBATCH --array=7-350%15        	# job array

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

#7. Generate BAM file statistics using samtools version 1.10 

cd /nvme/disk0/lecellier_data/WGS_NC_data/

#GBR
#INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"
#OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/BAM_statistics/"

#NC
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/BAM_statistics/"

#List files for which duplicates have been previoulsy removed 
FILES=($INDIR*MARKED_DUP.bam)

#Associate filename with array index 
FILENAME=${FILES[$((${SLURM_ARRAY_TASK_ID}-1))]}
BASE=$(basename $FILENAME)
BASE=${BASE%%_M*}

if [ ! -s "${OUTDIR}${BASE}-flagstats.txt" ] 
then
  
  start=`date +%s`
  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} : start processing
  
    #7.1 Generate index bai and statistics from BAM aligned reads 
    samtools index -@ 2 "${INDIR}${BASE}_MARKED_DUP.bam"
    #7.2 Output statistics 
    samtools idxstats --threads 2 "${INDIR}${BASE}_MARKED_DUP.bam" > "${OUTDIR}${BASE}-index_stats.txt"
    #Columns in output file are : reference sequence name, sequence lenght, mapped reads segment, unmapped reads segment 
    samtools coverage "${INDIR}${BASE}_MARKED_DUP.bam" > "${OUTDIR}${BASE}-coverage.txt"
    samtools flagstat --threads 2 -O tsv "${INDIR}${BASE}_MARKED_DUP.bam" > "${OUTDIR}${BASE}-flagstats.txt"	
  
  end=`date +%s`
  echo ${BASE} : Execution time was `expr $(( ($end - $start) / 60))` minutes.

else

  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} already processed

fi


#Count uniquely mapped reads to A.millepora chr 
#UNIQ_DIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Uniquely_mapped_reads_files/"
#
#FILES=($INDIR*Amillepora_MARKED_DUP.bam)
#FILES=("${FILES[@]:13:40}")
#
#start=`date +%s`
#for FILE in "${FILES[@]}"; do
#    BASE=$(basename $FILE)
#    BASE=${BASE%%_M*} 
#    #7.1 Generate uniquely mapped reads file 
#    samtools view -F 2048 -b "${INDIR}${BASE}_MARKED_DUP.bam" --threads 20 > "${UNIQ_DIR}${BASE}.uniquemapped.bam"
#    #samtools view -b -e '!([XA] | [SA])' --threads 20 "${INDIR}${BASE}_MARKED_DUP.bam" > "${UNIQ_DIR}${BASE}.uniquemapped.bam"
#    #7.2 Output statistics  
#    samtools coverage "${UNIQ_DIR}${BASE}.uniquemapped.bam" > "${UNIQ_DIR}${BASE}-coverage.txt"
#done
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.