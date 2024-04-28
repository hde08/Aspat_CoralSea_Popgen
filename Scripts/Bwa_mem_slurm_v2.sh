#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Bwa_mem_slurm_v2.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J BWA_MEM

### WALLTIME
#SBATCH -t 56:00:00

#MEMORY
#SBATCH --mem=150G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)

#Output and error directory
#SBATCH -o /home/hdenis/Slurm/bwa_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/bwa_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=10		# number of cores per job

#SBATCH --array=1-250%2        	# job array

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

#### 4. Map to reference genome using Bwa mem (version 0.7.17-r1188) with default parameters


cd /nvme/disk0/lecellier_data/WGS_GBR_data/
#mkdir Aligned_files/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Trimmed_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"

#List reference genomes to be tested 
REF_1="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amillepora_ncbi_dataset/data/GCA_013753865.1/GCA_013753865.1_Amil_v2.1_genomic.fna"
REF_2="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Aspathulata_ncbi_dataset/data/GCA_031770025.1/GCA_031770025.1_AGI_CSIRO_Aspa_v1_genomic.fna"
REFS=(${REF_1} ${REF_2})
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"

#Align only against A.millepora v3 reference genome 
REF_NAME="Amilleporav3"
#REF_NAME=${REF_3%%_ncbi*}
#REF_NAME=${REF_NAME#*genomes/}

#### 4.1. Index reference genomes 
#bwa index $REF_1 #A.millepora

#bwa index $REF_2 #A.spathulata

#bwa index $REF_3 #A.spathulata v3

FILES=($INDIR/*_R1_paired.fastq.gz)
#FILES=("${FILES[@]:75:500}")

#Associate filename with array index 
FILENAME=${FILES[$((${SLURM_ARRAY_TASK_ID}-1))]}
BASE=$(basename $FILENAME)
BASE=${BASE%%_R*}

if [ ! -s "${OUTDIR}${BASE}_pe_aln_${REF_NAME}_UNDEDUP.bam" ] 
then
  
  start=`date +%s`
  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} : start processing
  #Filename is associated with array index 
  
  #### 4.2. Map to reference genome
  bwa mem -t 10 ${REF_3} -o "${OUTDIR}${BASE}_pe_aln_${REF_NAME}.sam" "${INDIR}${BASE}_R1_paired.fastq.gz" "${INDIR}${BASE}_R2_paired.fastq.gz" 
  
  #### 4.3. Create sorted bam  
  samtools view -b --threads 10 "${OUTDIR}${BASE}_pe_aln_${REF_NAME}.sam" | samtools sort --threads 10  -O BAM -o "${OUTDIR}${BASE}_pe_aln_${REF_NAME}_UNDEDUP.bam"
  
  end=`date +%s`
  echo ${BASE} : Execution time was `expr $(( ($end - $start) / 60))` minutes.
  
  ##### 4.4. removed temporary files
  rm "${OUTDIR}${BASE}_pe_aln_${REF_NAME}.sam"	
  
  #Remove the trimmed input files (temporarly saved on tgt_bck2)
  #rm "${INDIR}${BASE}_R1_paired.fastq.gz" 
  #rm "${INDIR}${BASE}_R2_paired.fastq.gz"
  #rm "${INDIR}${BASE}_R1_unpaired.fastq.gz"
  #rm "${INDIR}${BASE}_R2_unpaired.fastq.gz"
  #rm "${INDIR}${BASE}_sum.txt"
  
  
    
else

  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} already processed

fi


