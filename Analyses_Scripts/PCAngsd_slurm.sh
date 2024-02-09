#!/bin/bash
#chmod u+x /home/hdenis/Coral-Genomics/Analyses_Scripts/PCAngsd_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J PCAngsd

### WALLTIME
#SBATCH -t 56:00:00

### MPI TASKS (cores)
#SBATCH -n 8

#MEMORY
#SBATCH --mem=10G

#Output and error directory
#SBATCH -o /home/hdenis/Slurm/outFile_%j.out

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


#### 9. Perform PCA on ANGSD output (beagle file) 
#mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/Analyses_outputs/

cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Analyses_outputs/"

#python /home/hdenis/Programs/pcangsd/pcangsd/pcangsd.py -h

start=`date +%s`
pcangsd --beagle "${INDIR}GBR_sub_test.beagle.gz" --out "${OUTDIR}pcangsd" --maf 0.05 --threads $NPROCS
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Parameters 
#--maf 0.05 (minimum minor allele frequency)

#Perform Bayesian hiearchical clustering admixture analyses
#start=`date +%s`
#python /home/hdenis/Programs/pcangsd/pcangsd/pcangsd.py -beagle "${INDIR}GBR_sub_test.input.beagle.gz" -out "${OUTDIR}admix" --maf 0.05 --admix -threads $NPROCS
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#To run this script 
#sbatch /home/hdenis/Coral-Genomics/Analyses_Scripts/PCAngsd_slurm.sh --output=/home/hdenis/Slurm/output.

