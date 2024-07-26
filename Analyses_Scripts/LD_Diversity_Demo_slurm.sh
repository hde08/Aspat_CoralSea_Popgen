#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/LD_Diversity_Demo_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J Div_Demo

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=300G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/vcftool_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/vcftool_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=20		# number of cores per job

#SBATCH --array=1        	# job array

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

#11 Merge and filter vcf files 

cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/Vcf_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Analyses_outputs/Diversity_Demo/"

#A.millepora v3 reference genome 
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"
REF_NAME="Amilleporav3"

#12.1 Estimate LD decay (set of filtered SNPs produced by GATK Variant filtration)
#Retain 1% of variants
#Compute r2 in a 100kb distance 
#Compute r2 for a maximum of 999999 variants 
#Prevent error due to unrecognized chromosome codes 
start=`date +%s`
echo `Start computing LD decay`
plink --allow-extra-chr --ld-window 999999 --ld-window-kb 100 --ld-window-r2 0 --out "${INDIR}aspat_clean_all_chr_SNP_nofilt_filtered1" --r2 --thin 0.01 --vcf "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" --threads 4 --const-fid 0  
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.  

#12.2 Performed PCA on filtered and pruned set of SNPs
./plink --allow-extra-chr --extract "${INDIR}aspat_clean_all_chr_SNP_filtered3_LD_pruned0.15.prune.in" --out "${OUTDIR}aspat_clean_all_chr_SNP_filtered3_LD_pruned0.15_pca" --vcf "${INDIR}aspat_clean_all_chr_SNP_filtered3.vcf.gz" --pca


#12.3 Estimate Fst between each pair of reef (set of filtered SNPs in LE)

## run bcftools query and create population files using R script 
bcftools query -l "${INDIR}aspat_clean_all_chr.vcf.gz" > "${INDIR}aspat_vcf_sample_names" 

#Compute pairwise Fst between reefs 
POP_FILES=(${OUTDIR}*samples*)
for POP1 in ${POP_FILES[@]};
  for POP2 in ${POP_FILES[@]};
    do
      if (( $POP1 != $POP2 )); then
        echo `Computing pairwise Fst for $POP1 and $POP2`
        start=`date +%s`
        plink --vcf "${INDIR}aspat_clean_all_chr_SNP_filtered2.vcf.gz" --extract pruned.prune.in --weir-fst-pop ${POP1} --weir-fst-pop ${POP2} --out ${OUTDIR}${POP1}_${POP2}_pairwise
        end=`date +%s`
        echo Execution time was `expr $(( ($end - $start) / 60))` minutes.
      else
        `Same population : skipping `
      fi
  done
done

#12.3 Estimate Ne, Pi and Heterozygosity per group (1 and 2)
#Need to recreate group1-2 files without clones 


#12.4 Assess changes in Ne over time -> not yet 



#12.5 Assess migration rates between reefs








