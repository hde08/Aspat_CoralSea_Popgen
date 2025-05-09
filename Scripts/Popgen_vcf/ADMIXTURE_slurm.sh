#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Popgen_vcf/ADMIXTURE_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J ADMIXTURE

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=10G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/admix_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/admix_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=3		# number of cores per job

#SBATCH --array=1-6%2            	# job array 1-6%2    

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

########################################################## ADMIXTURE ###################################################################
#Perform ADMIXTURE analysis 

#1. Run ADMIXTURE reduced GBR+NC Aspat/Amil dataset 
cd /nvme/disk0/lecellier_data/GBR_NC_phylo_subset/
INDIR="/nvme/disk0/lecellier_data/GBR_NC_phylo_subset/Vcf_files/"
BASE="gbr_nc_aspat_amil_phylo_subset"
NEWBASE="gbr_nc_aspat_clean_phylo_subset"
OUTDIR="/nvme/disk0/lecellier_data/GBR_NC_phylo_subset/Popgen_outputs/"

cd "/nvme/disk0/lecellier_data/GBR_NC_phylo_subset/Popgen_outputs/"

K=${SLURM_ARRAY_TASK_ID}

source activate base
conda activate admixture 

Now that it has been recoded, Run admixture for different values of K 
admixture -B -j3 --cv "${INDIR}${BASE}_all_chr_SNP_filtered_2_recoded_chr_nosim_LD_pruned0.1.bed" $K | tee "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_recoded_chr_nosim_LD_pruned0.1.log${K}.out"

#Save CV error for each K in single file 
grep -h CV ${OUTDIR}${BASE}_all_chr_SNP_filtered_2_recoded_chr_nosim_LD_pruned0.1.log*.out > "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_recoded_chr_nosim_LD_pruned0.1.Kselection.txt"

#Save Log-likelihood for each K in single file 
grep -h '^Loglikelihood' ${OUTDIR}${BASE}_all_chr_SNP_filtered_2_recoded_chr_nosim_LD_pruned0.1.log*.out > "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_recoded_chr_nosim_LD_pruned0.1.Loglikelihood.txt"

#Save Log-likelihood across each bootstrap replicates for each value of K in single file
for K in {1..10}
  do
  grep -h '^1 ' ${OUTDIR}${BASE}_all_chr_SNP_filtered_2_recoded_chr_nosim_LD_pruned0.1.log$K.out | tail -n +3 > "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_recoded_chr_nosim_LD_pruned0.1.${K}.boostrap_Loglikelihood.txt"
done

#2. Run ADMIXTURE on whole GBR+NC datasets, Aspat samples only 
#Including different amount of missing data
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Vcf_files/"
BASE="aspat_clean"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/Admixture/"

cd "/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/Admixture/"

K=${SLURM_ARRAY_TASK_ID}


MISS=(05 10 20 50)
MAF="0.05"

for M in ${MISS[@]};
    do
    admixture -B -j3 --cv "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2.bed" $K | tee "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2.log${K}.out"
    
    #Save CV error for each K in single file 
    grep -h CV ${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2.log${K}.out >> "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2.Kselection.txt"

    #Save Log-likelihood for each K in single file 
    grep -h '^Loglikelihood' ${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2.log*.out > "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2.Loglikelihood.txt"

    #Save Log-likelihood across each bootstrap replicates for each value of K in single file
    for K in {1..6}
      do
      grep -h '^1 ' ${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2.log${K}.out | tail -n +3 > "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2.log${K}.boostrap_Loglikelihood.txt"
    done
    
done