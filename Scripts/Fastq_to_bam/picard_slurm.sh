#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/picard_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J Picard

### WALLTIME
#SBATCH -t 56:00:00

#MEMORY
#SBATCH --mem=10G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)

#Output and error directory
#SBATCH -o /home/hdenis/Slurm/picard_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/picard_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=1		# number of cores per job

#SBATCH --array=1-6%20        	# job array

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


#5. Mark and delete PCR duplicates using picard version 2.27.5 

#Recommanded to generate before after quality plot to check base recalibration 
cd /nvme/disk0/lecellier_data/WGS_NC_data/
mkdir /nvme/disk0/lecellier_data/WGS_NC_data/BAM_statistics/

#GBR
#INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"
#OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/BAM_statistics/"

#NC
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/BAM_statistics/"

#List files that have been aligned and indexed with bwa mem / samtools 
FILES=($INDIR*_UNDEDUP.bam)

#Associate filename with array index 
FILENAME=${FILES[$((${SLURM_ARRAY_TASK_ID}-1))]}
BASE=$(basename $FILENAME)
BASE=${BASE%%_U*}

#Load correct java version (>= java v17)
#ml load gcc-11.4.1/jdk/17.0.2_openjdk-rplh5ry

if [ ! -s "${INDIR}${BASE}_MARKED_DUP.bam" ] 
then
  
  start=`date +%s`
  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} : start processing
  
  #5.1 Add read group information
  if echo ${BASE} | grep -q 'RRAP'; then
  #GBR
  #Get file info from other repertory
  BASE_REP="${BASE#*-}"
  BASE_REP=$(ls /tgt_bck2/data_lecellier/BAM_files/  | grep $BASE)
  BASE_REP=${BASE%%.u*}
  
  singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif AddOrReplaceReadGroups INPUT="${INDIR}${BASE}_UNDEDUP.bam" OUTPUT="${INDIR}${BASE}_UNDEDUP_RG.bam" RGID=$(echo ${BASE_REP} | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/") RGLB=$(echo ${BASE_REP} | cut -d"-" -f1,2,3) RGPL=ILLUMINA RGPU=$(echo ${BASE_REP} | cut -d"_" -f2) RGSM=$(echo ${BASE_REP} | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/")
  
  else if echo ${BASE} | grep -q 'RA-'; then
    #NC - Reef Adapt samples
    singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif AddOrReplaceReadGroups INPUT="${INDIR}${BASE}_UNDEDUP.bam" OUTPUT="${INDIR}${BASE}_UNDEDUP_RG.bam" RGID=${BASE} RGLB="RECOVER_2023" RGPL=ILLUMINA RGPU="L1" RGSM=${BASE}
      
  else
    #NC - Recover samples 
    singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif AddOrReplaceReadGroups INPUT="${INDIR}${BASE}_UNDEDUP.bam" OUTPUT="${INDIR}${BASE}_UNDEDUP_RG.bam" RGID=${BASE} RGLB="REEF_ADAPT_2021" RGPL=ILLUMINA RGPU="L1" RGSM=${BASE}
  fi

  #RGID : unique read group identifier
  #RGLB : library index
  #RGPL : platform name
  #RGPU : flowcell, lane, sample barcode 
  #RGSM : sample name (=RGID in this case)
  
  #5.2 Mark and remove duplicate
  singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif MarkDuplicates INPUT="${INDIR}${BASE}_UNDEDUP_RG.bam" OUTPUT="${INDIR}${BASE}_MARKED_DUP.bam" METRICS_FILE="${OUTDIR}${BASE}_marked_dup_metrics.txt" REMOVE_DUPLICATES=true TMP_DIR="${INDIR}" VALIDATION_STRINGENCY=LENIENT
  
  end=`date +%s`
  echo ${BASE} : Execution time was `expr $(( ($end - $start) / 60))` minutes.
  
  rm "${INDIR}${BASE}_UNDEDUP_RG.bam"

else

  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} already processed

fi


#Trash 
##Remove UNDEDUP files that have been processed
#INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"
#
##List files that have been aligned and indexed with bwa mem / samtools 
#MARKED_FILES=($INDIR*_MARKED_DUP.bam)
#for FILE in ${MARKED_FILES[@]}; do
#  BASE=$(basename $FILE)
#  BASE=${BASE%%_M*}
#  rm "${INDIR}${BASE}_UNDEDUP.bam"
#done


##Create file IDs
#>/nvme/disk0/lecellier_data/WGS_GBR_data/Raw_data_processing/aln_ids.txt
#for FILE in ${FILES[@]}; do
#	BASE=$(basename $FILE)
#	BASE=${BASE%%_U*} 
#	echo ${BASE} >> /nvme/disk0/lecellier_data/WGS_GBR_data/Raw_data_processing/aln_ids.txt
#done 

#######################################################################################

#View read group 
#samtools view -H "${INDIR}${BASE}_UNDEDUP.bam" | grep '^@RG