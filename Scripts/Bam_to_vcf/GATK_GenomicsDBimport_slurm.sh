#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/GATK_GenomicsDBimport_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J GATK

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=20G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/gatk_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/gatk_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=2		# number of cores per job

#SBATCH --array=1-36%15        	# job array

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
#module load gcc-11.4.1/jdk/17.0.2_openjdk-rplh5ry


# This job's working directory :
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

ulimit -s unlimited

########################################################## CREATE GENOMICDB ###################################################################
#This scripts consolidates GVCF using GenomicsDBImport

cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/Vcf_files/"

#A.millepora v3 reference genome 
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"
REF_NAME="Amilleporav3"


#1. Create map file of following format
#sample1      sample1.vcf.gz
#sample2      sample2.vcf.gz
#sample3      sample3.vcf.gz

GVCF_FILES=($INDIR*.g.vcf.gz)

echo -n > "${INDIR}aspat_gvcf_clean.sample_map"
for FILE in ${GVCF_FILES[@]}; do
    NAME="$(basename $FILE)"
    NAME="${NAME%.g.vcf*}"
    printf "%s\t%s\n" "$NAME" "$FILE" >> "${INDIR}aspat_gvcf_clean.sample_map"
done


#2. Read from a list of intervals  (total of 36 intervals)
#Chromosomes 1-4 divided in 4 intervals
#Chromosomes 5-14 divided in 2 inverals 
INTERVALS_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/chromosome_intervals.txt"
INTERVALS_NAMES_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/chromosome_interval_names.txt"

readarray -t INTERVALS <  $INTERVALS_FILE
readarray -t INTERVALS_NAMES <  $INTERVALS_NAMES_FILE

#Match slurm array ID to interval 
CONTIG=${INTERVALS[$((${SLURM_ARRAY_TASK_ID}-1))]}
CONTIG_NAME=${INTERVALS_NAMES[$((${SLURM_ARRAY_TASK_ID}-1))]}

#3. Create GenomicDB gatk v4.5.0.0

start=`date +%s`
echo Array Id : ${SLURM_ARRAY_TASK_ID} Chromosome : ${CONTIG} : start processing

singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx15g -Xms4g" GenomicsDBImport --genomicsdb-workspace-path "${INDIR}GenomicDB/${CONTIG_NAME}" --batch-size 50 -L $CONTIG --sample-name-map "${INDIR}aspat_gvcf_clean.sample_map" --tmp-dir /nvme/disk0/lecellier_data/WGS_GBR_data/tmp --reader-threads 7 --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader true 

end=`date +%s`
echo ${CONTIG} : Execution time was `expr $(( ($end - $start) / 60))` minutes.


