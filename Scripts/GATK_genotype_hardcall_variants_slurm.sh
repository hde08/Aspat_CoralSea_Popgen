#!/bin/bash

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J GATK

### WALLTIME
#SBATCH -t 56:00:00

### MPI TASKS (cores)
#SBATCH -n 8

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


#### 10. Perform genotype hard callig variants using gatk v4.5.0.0
#mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/
cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/"

#Create filelist with full path of each BAM file (one filepath per line)
FILES=($INDIR*Amillepora_MARKED_DUP.bam)
printf "%s\n" "${FILES[@]}" > $INDIR/bam.filelist.txt
#Get number of files
N_FILES="${#FILES[@]}"
MIN_N=$((95*N_FILES/100))

#A.millepora reference genome 
REF_1="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amillepora_ncbi_dataset/data/GCA_013753865.1/GCA_013753865.1_Amil_v2.1_genomic.fna"

#10.1 Call genotypes per sample using Haplotype Caller 

gatk HaplotypeCaller --h

gatk HaplotypeCaller --input bam --output $OUTDIR 

#Parameters 
--input,-I <GATKPath>         BAM/SAM/CRAM file containing reads  This argument must be specified at least once.
                              Required.

--output,-O <GATKPath>        File to which variants should be written  Required.

--reference,-R <GATKPath>     Reference sequence file  Required.

- generates intermediate GVCF file (-ERC GVCF mode


#All in one single command line

#Edit the outdir path
start=`date +%s`
angsd -bam $INDIR/bam.filelist.txt -out "${OUTDIR}GBR_sub_test" -ref $REF_1 -uniqueOnly 0 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd MIN_N -setMinDepth 3  -doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -nThreads $NPROCS
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Note on beagle file coding
#0=A, 1=C, 2=G, 3=T
