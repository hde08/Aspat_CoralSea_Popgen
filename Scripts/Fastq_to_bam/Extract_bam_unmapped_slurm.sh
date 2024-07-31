#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Extract_bam_unmapped_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J Extract_bam

### WALLTIME
#SBATCH -t 56:00:00

#MEMORY
#SBATCH --mem=15G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)

#Output and error directory
#SBATCH -o /home/hdenis/Slurm/extr_bam_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/extr_bam_%A_%a.e     # standard error

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

#6. Extract unmapped reads (non host reads) from bam files to separate files 

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

if [ ! -s "${INDIR}${BASE}.unmapped.bam" ] 
then
  
  start=`date +%s`
  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} : start processing
  
    # generate sample names
    echo ${BASE} >> ${OUTDIR}/sample_name.txt		
    # generate read counts
    samtools view -c --threads 2 "${INDIR}${BASE}_MARKED_DUP.bam" >> ${OUTDIR}allCounts.txt		
    # generate bam with unmapped reads = symbiont, microbes and other reads 
    samtools view -b -f 12 -F 256 --threads 2 "${INDIR}${BASE}_MARKED_DUP.bam" > "${INDIR}${BASE}.unmapped.bam"
    # generate unmapped read counts
    samtools view -c --threads 2 "${INDIR}${BASE}.unmapped.bam" >> "${OUTDIR}unmappedCounts.txt"		
  
  end=`date +%s`
  echo ${BASE} : Execution time was `expr $(( ($end - $start) / 60))` minutes.

else

  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} already processed

fi

# Notes:
# -f 4: extract all unmapped reads
# -f12: extract only the reads where read 1 is unmapped AND read 2 is unmapped (= both mates are unmapped). Applies only to paired reads.

