#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/GATK_HaplotypeCaller_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J GATK

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=10G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/gatk_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/gatk_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=1		# number of cores per job

#SBATCH --array=1-823%30        	# job array

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

### Ressources : Load compatible java version
module load gcc-11.4.1/jdk/17.0.2_openjdk-rplh5ry


# This job's working directory :
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

ulimit -s unlimited


#### 11. Perform genotype hard callig variants using gatk v4.5.0.0
#mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/
cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/"

#Get filelist of clean aspat files (good quality samples, no bad taxID)
readarray -t BAM_FILES <  $INDIR/aspat_bam_clean.filelist.txt


#Match slurm array ID to file 
FILE=${BAM_FILES[$((${SLURM_ARRAY_TASK_ID}-1))]}
BASE=$(basename $FILE)
BASE=${BASE%%_M*}

#A.millepora v3 reference genome 
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"
REF_NAME="Amilleporav3"

#List with chromosome names 
readarray -t CHROMOSOMES < "/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/chromosomes_header.txt"

#10.1 Create samtools genome index and GATK reference dictionary

#cd /nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/
#
#gatk CreateSequenceDictionary -R $REF_3
#
#samtools faidx $REF_3


#10.2 Call genotypes per sample using Haplotype Caller gatk v4.5.0.0

#Genotype calling on complete reference (chr + scaffolds)
OPTIONS="-ERC GVCF"

if [ ! -s "$OUTDIR$BASE.g.vcf.gz" ] 
then
  
  start=`date +%s`
  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} : start processing
  
  singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx8g" HaplotypeCaller --input $FILE --output "$OUTDIR$BASE.g.vcf.gz" --reference $REF_3 $OPTIONS
  
  end=`date +%s`
  echo ${BASE} : Execution time was `expr $(( ($end - $start) / 60))` minutes.
   
else

  echo Array Id : ${SLURM_ARRAY_TASK_ID} File : ${BASE} already processed

fi

#Validate format 
#gatk ValidateVariants -V /nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/RRAP-ECT01-2022-Aspat-CBHE-1718_L1_pe_aln_Amilleporav3.g.vcf.gz --validation-type-to-exclude ALL


  

