#!/bin/bash

### memento:
###           sbatch POE_BEACH_slurm.sh
###           squeue  (qstat)
###           scancel JOB_ID  (qdel)
###           sinfo           (equivalent qstat -q avec reservation en cours)
###           scontrol show node  (equi pbsnodes -a)
###           scontrol show job JOB_ID  (detail ton job)

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J ANGSD

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

### Ressources
#module load intel/2022.0.2
#module load openmpi/1.10.7
#module load netcdf-fortran-4.5.3-intel-2021.5.0-gzuyinc
#module load netcdf-c-4.8.1-intel-2021.5.0-j7rukz5

# This job's working directory :
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

ulimit -s unlimited


#### 8. Perform SNP and genotype calling in ANGSD v0.1.17

cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/"
STATDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Raw_data_processing/BAM_statistics/"

#Create filelist with full path of each BAM file (one filepath per line)
FILES=($INDIR*Amillepora_MARKED_DUP.bam)
printf "%s\n" "${FILES[@]}" > $INDIR/bam.filelist.txt
#Get number of files
N_FILES="${#FILES[@]}"
MIN_N=$((95*N_FILES/100))

#Get median and sd depth value across bam files 
#Store genome-wide mean depth for each fle 
> $STATDIR/bamfile_depth.txt
> $STATDIR/bamfile_sd_depth.txt
DEPTH_FILES=($STATDIR*Amillepora*coverage.txt)
for FILE in "${DEPTH_FILES[@]}"; do
    MEAN_DEPTH=$(awk '{sum+=$7} END { print sum/NR}' $FILE)
    SD_DEPTH=$(awk '{delta = $7 - avg; avg += delta / NR; mean2 += delta * ($7 - avg); } END { print sqrt(mean2 / (NR-1)); }' $FILE)
    echo $MEAN_DEPTH >> $STATDIR/bamfile_depth.txt
    echo $SD_DEPTH>> $STATDIR/bamfile_sd_depth.txt
done

#Set max depth as the median genome-wide depth + 1 SD (following Lou 2021 guidelines)
MEDIAN_DEPTH=$(sort -k1 -n $STATDIR/bamfile_depth.txt | awk '{arr[NR]=$1}
   END { if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}')
MEDIAN_SD=$(sort -k1 -n $STATDIR/bamfile_sd_depth.txt | awk '{arr[NR]=$1}
   END { if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}')
MAX_DEPTH=`echo "$MEDIAN_DEPTH + $MEDIAN_SD" | bc -l`


#Reference genomes -> Pick only one of the two
REF_1="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amillepora_ncbi_dataset/data/GCA_013753865.1/GCA_013753865.1_Amil_v2.1_genomic.fna"
REF_2="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Aspathulata_ncbi_dataset/data/GCA_031770025.1/GCA_031770025.1_AGI_CSIRO_Aspa_v1_genomic.fna"

#8.1 Re-index reference using samtools
#samtools faidx $REF_1

#8.2 Filtering of polymorphic sites

#Parameters :
#--minQ 30 base quality >30
#--minMapQ mapping quality >30
#--minInd MIN_N sites with data from at least 95% of individuals
#--setMinIndDepth 3 minimum of 3 reads
#--uniqueOnly 0 keep reads with multiple hits
#--only_proper_pairs 1 only paired reads (already filtered)
#--remove_bads 1 (remove reads flags with 256)

#Not used parameters
# -setMaxDepth 30
#-trim 0 -C 50 -baq 1

#8.3 Infer genotypes likelihood and major/minor allele frequencies

#Parameters
#--GL 1 (SAMtools model)
#--doGlf 2 (output in beagle format)
#--doMaf 2 (fixed major, unkown minor, consider switching to 4 )
#---doMajorMinor (infer from genotype likelihoods)
#--SNP_pval 1e-6 (variant likelihood ratio <0.0000001)

#All in one single command line

#Edit the outdir path
#See the difference between setMinDepth and setMinDepthInd
#Look to set a max depth as well 
#-setMaxDepthInd 
#-minInd

#Check the difference between the different version or major minor inference 

start=`date +%s`
angsd -bam $INDIR/bam.filelist.txt -out "${OUTDIR}GBR_sub_test" -ref $REF_1 -uniqueOnly 0 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd MIN_N -setMinDepthInd 3  -doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -nThreads $NPROCS -setMaxDepthInd $MAX_DEPTH
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Note on beagle file coding
#0=A, 1=C, 2=G, 3=T


