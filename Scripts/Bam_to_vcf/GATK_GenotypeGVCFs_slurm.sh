#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/GATK_GenotypeGVCFs_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J GATK

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=70G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/gatk_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/gatk_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=4		# number of cores per job

#SBATCH --array=26-36%4        	# job array

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

##10.4 Perform joint genotyping with Genotype GVCF

cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/Vcf_files/"

#A.millepora v3 reference genome 
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"
REF_NAME="Amilleporav3"


#10.41 Read from a list of intervals  (total of 36 intervals)
#Chromosomes 1-4 divided in 4 intervals
#Chromosomes 5-14 divided in 2 inverals 
INTERVALS_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/chromosome_intervals.txt"
INTERVALS_NAMES_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/chromosome_interval_names.txt"

readarray -t INTERVALS <  $INTERVALS_FILE
readarray -t INTERVALS_NAMES <  $INTERVALS_NAMES_FILE

#Match slurm array ID to file 
CONTIG=${INTERVALS[$((${SLURM_ARRAY_TASK_ID}-1))]}
CONTIG_NAME=${INTERVALS_NAMES[$((${SLURM_ARRAY_TASK_ID}-1))]}

#OPTIONS : 
#--include-non-variant-sites true : output non-variable sites
#--only-output-calls-starting-in-intervals : makes sure that only sites starting within the specified intervals are reported
#However sites whose end is located outside the interval are reported as long as they start within

start=`date +%s`
echo Array Id : ${SLURM_ARRAY_TASK_ID} Chromosome : ${CONTIG} : start joint genotyping

singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx58g" GenotypeGVCFs -R $REF_3 -V "gendb://${INDIR}GenomicDB/${CONTIG_NAME}" -O "${OUTDIR}aspat_clean_${CONTIG_NAME}.vcf.gz" --tmp-dir /nvme/disk0/lecellier_data/WGS_GBR_data/tmp --include-non-variant-sites true -L $CONTIG --only-output-calls-starting-in-intervals true

end=`date +%s`
echo ${CONTIG} : Execution time was `expr $(( ($end - $start) / 60))` minutes.


#Check changes in variants processed over time 
#cat /home/hdenis/Slurm/gatk_7338_1.e | grep 'ProgressMeter' > /nvme/disk0/lecellier_data/WGS_GBR_data/Vcf_performance.txt


