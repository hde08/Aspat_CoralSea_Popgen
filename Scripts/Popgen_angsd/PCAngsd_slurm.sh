#!/bin/bash
#chmod u+x /home/hdenis/Coral-Genomics/Analyses_Scripts/PCAngsd_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J PCANGSD

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=50G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/pcangsd_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/pcangsd_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=4		# number of cores per job

#SBATCH --array=1-5%5        	# job array

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

########################################################## POPGEN ANGSD ###################################################################
#Perform PCA and ADMIXTURE analyses using genotype likelihoods (beagle file)

#1. GBR samples only 
cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Analyses_outputs/"

#Associate slurm array index with K value 
K=${SLURM_ARRAY_TASK_ID}

#Parameters 
#--maf 0.05 (minimum minor allele frequency)

## 1.1 Perform PCA and ADMIXTURE on all GBR samples (including Ahya and Amil outgroups) to assess possible misID 
Remove some Amil samples that were uncorrectly identified 


if (( $K == 1 )); then
  #Compute covariance matrix 
  start=`date +%s`
  echo Start computing covariance matrix 
  
  pcangsd --beagle "${INDIR}GBR_allsamples_filt_uniq_all_chr.beagle.gz" --out "${OUTDIR}GBR_allsamples_filt_uniq_all_chr" --maf 0.05 --threads 4 --filter "/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/aspat_amil_ahya_subset.filelist.txt"
  
  end=`date +%s`
  echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 
  
  else
  #Perform Bayesian hiearchical clustering admixture analyses
    start=`date +%s`
    echo Start computing admixture K=$K
    
    pcangsd --beagle "${INDIR}GBR_allsamples_filt_uniq_all_chr.beagle.gz" --out "${OUTDIR}GBR_allsamples_filt_uniq_all_chr_K${K}" --maf 0.05 --admix --admix_K $K --threads 4 --filter "/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/aspat_amil_ahya_subset.filelist.txt"
    
    end=`date +%s`
    echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 
  fi

#1.2 Re-run PCA and ADMIXTURE on Aspat samples only after removal of 4 mis-identifications and 49 putative clones (based on IBS)
if (( $K == 1 )); then
  #Compute covariance matrix 
  start=`date +%s`
  echo Start computing covariance matrix 
  
  pcangsd --beagle "${INDIR}GBR_aspatsamples_filt_05mis_uniq_all_chr.beagle.gz" --out "${OUTDIR}GBR_aspatsamples_filt_05mis_uniq_all_chr_noclones" --maf 0.05 --threads 4 --filter "/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/aspat_bam_noclones.filelist.txt"
  
  end=`date +%s`
  echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 
  
  else
  #Perform Bayesian hiearchical clustering admixture analyses
    start=`date +%s`
    echo Start computing admixture K=$K
    
    pcangsd --beagle "${INDIR}GBR_aspatsamples_filt_05mis_uniq_all_chr.beagle.gz" --out "${OUTDIR}GBR_aspatsamples_filt_05mis_uniq_all_chr_noclones_K${K}" --maf 0.05 --admix --admix_K $K --threads 4 --filter "/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/aspat_bam_noclones.filelist.txt"
    
    end=`date +%s`
    echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 
  fi

#2. NC samples only 
cd /nvme/disk0/lecellier_data/WGS_NC_data/
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/ANGSD_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Analyses_outputs/"

#Associate slurm array index with K value 
K=${SLURM_ARRAY_TASK_ID}

#Parameters 
#--maf 0.05 (minimum minor allele frequency)

### 9.1 Perform PCA and ADMIXTURE on all NC samples
if (( $K == 1 )); then
  #Compute covariance matrix 
  start=`date +%s`
  echo Start computing covariance matrix 
  
  pcangsd --beagle "${INDIR}NC_allsamples_filt_05mis_uniq_all_chr.beagle.gz" --out "${OUTDIR}NC_allsamples_filt_05mis_uniq_all_chr" --maf 0.05 --threads 4 
  
#  --filter "/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/aspat_amil_ahya_subset.filelist.txt"
  
  end=`date +%s`
  echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 
  
  else
  #Perform Bayesian hiearchical clustering admixture analyses
    start=`date +%s`
    echo Start computing admixture K=$K
    
    pcangsd --beagle "${INDIR}NC_allsamples_filt_05mis_uniq_all_chr.beagle.gz" --out "${OUTDIR}NC_allsamples_filt_05mis_uniq_all_chr_K${K}" --maf 0.05 --admix --admix_K $K --threads 4 
    
#    --filter "/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/aspat_amil_ahya_subset.filelist.txt"
    
    end=`date +%s`
    echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 
  fi


#3. GBR+NC samples 
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/ANGSD_files/"

#3.1 PCA and ADMIXTURE on NC + Aspat/Amil GBR samples 
if (( $K == 1 )); then
  #Compute covariance matrix 
  start=`date +%s`
  echo Start computing covariance matrix 
  
  pcangsd --beagle "${INDIR}GBR_NC_aspat_amil_samples_05mis_all_chr.beagle.gz" --out "${OUTDIR}GBR_NC_aspat_amil_samples_05mis_all_chr" --maf 0.05 --threads 4 
  end=`date +%s`
  echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 
  
  else
  #Perform Bayesian hiearchical clustering admixture analyses
    start=`date +%s`
    echo Start computing admixture K=$K
    
    pcangsd --beagle "${INDIR}GBR_NC_aspat_amil_samples_05mis_all_chr.beagle.gz" --out "${OUTDIR}GBR_NC_aspat_amil_samples_05mis_all_chr_K${K}" --maf 0.05 --admix --admix_K $K --threads 4  
    end=`date +%s`
    echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 
  fi
  
  
#3.2 PCA and ADMIXTURE on NC + GBR samples excluding Amil outgroup
if (( $K == 1 )); then
#Compute covariance matrix 
start=`date +%s`
echo Start computing covariance matrix 

pcangsd --beagle "${INDIR}GBR_NC_aspat_amil_samples_05mis_all_chr.beagle.gz" --out "${OUTDIR}GBR_NC_aspat_samples_05mis_all_chr" --maf 0.05 --threads 4 --filter "${INDIR}GBR_NC_aspat_samples.filterlist.txt"


end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 

else
#Perform Bayesian hiearchical clustering admixture analyses
  start=`date +%s`
  echo Start computing admixture K=$K
  
  pcangsd --beagle "${INDIR}GBR_NC_aspat_amil_samples_05mis_all_chr.beagle.gz" --out "${OUTDIR}GBR_NC_aspat_samples_05mis_all_chr_K${K}" --maf 0.05 --admix --admix_K $K --threads 4 --filter "${INDIR}GBR_NC_aspat_samples.filterlist.txt"
 
  end=`date +%s`
  echo Execution time was `expr $(( ($end - $start) / 60))` minutes.; 
fi
  

