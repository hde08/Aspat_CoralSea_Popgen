#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Analyses_Scripts/IBS_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J IBS

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=400G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/IBS_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/IBS_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=15		# number of cores per job

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

### Ressources
#module load intel/2022.0.2
#module load openmpi/1.10.7
#module load netcdf-fortran-4.5.3-intel-2021.5.0-gzuyinc
#module load netcdf-c-4.8.1-intel-2021.5.0-j7rukz5

# This job's working directory :
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

ulimit -s unlimited

########################################################## SAMPLES IBS ###################################################################
#Identify clones and related individuals using IBS method in ANGSD (Identity by state)
#The analysis is done separately for each genomic cluster identified based on PCA and Admixture results

#GBR
cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/"

#NC
cd /nvme/disk0/lecellier_data/WGS_NC_data/
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/ANGSD_files/"

#A.millepora v3 reference genome 
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"
REF_NAME="Amilleporav3"

mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/Analyses_outputs/IBS_results/

#Match slurm array index with group number 
FILE="${INDIR}aspat_bam_group${SLURM_ARRAY_TASK_ID}.filelist.txt"

#Get number of files
N_FILES=$(wc -l < $FILE)
MIN_N=$((95*N_FILES/100))

#1. Run IBS analysis on each cluster (angsd) separately

FILTERS="-minMapQ 30 -minQ 30 -minInd MIN_N -setMinDepthInd 3 -uniqueOnly 1 -remove_bads 1 -doSNPstat -SNP_pval 1e-6 -minMaf 0.05 -only_proper_pairs 1"

TODO="-doIBS 1 -makeMatrix 1 -GL 1 -doMaf 2 -doMajorMinor 1 -doCounts 1"

start=`date +%s`
echo Start computing IBS 
angsd -bam $FILE  -out "Analyses_outputs/IBS_results/NC_aspat_group${SLURM_ARRAY_TASK_ID}_ibs05" -nThreads 5 $FILTERS $TODO
angsd -bam $FILE  -out "Analyses_outputs/IBS_results/GBR_aspat_group${SLURM_ARRAY_TASK_ID}_ibs05" -nThreads 15 $FILTERS $TODO
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.





