#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Popgen_vcf/RAxML_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J RAxML

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=130G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/raxml_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/raxml_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=15		# number of cores per job

#SBATCH --array=1        	# job array 1-6%2      

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

########################################################## RAXML TREE ###################################################################
#Run RAxML on host genomic data  (RAxML version 8.2.12)

#We use a subset of 70 samples per genomic cluster corresponding to the minimal sample size 

#Filtration of input file
#bi-allelic sites
#05,10 or 20 missing data
#Minor allele frequency >0.05
#LD<0.2

MISS=(05 10 20)
MAF="0.05"

INDIR="/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Vcf_files/"
BASE="aspat_clean"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/RAxML/"

cd $OUTDIR 

N=70
##Make sure we only take samples for which we have been able to compute the kmer distance matrix for the symbionts to create a tanglegram
##Make sure we also take both samples from inshore and offshore central/gbr since they host distinct symbiont communities 
SYMBIONTS_SUBSET="/nvme/disk0/lecellier_data/Symbionts/Symbionts_clean_kmer.filelist.txt"
SYMBIONTS_INSHORE_CN_GBR="/nvme/disk0/lecellier_data/Symbionts/Symbionts_cluster2_C50_samples.filelist.txt"


> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"

cat "${INDIR}${BASE}_filtered2_Group1_noclones.filelist.txt" | grep -F -f "${SYMBIONTS_SUBSET}"  | grep -F -v -f "${INDIR}${BASE}_PCA_outliers.txt" | grep -F -f "${SYMBIONTS_INSHORE_CN_GBR}" | shuf -n 35 >> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"
cat "${INDIR}${BASE}_filtered2_Group1_noclones.filelist.txt" | grep -F -f "${SYMBIONTS_SUBSET}"  | grep -F -v -f "${INDIR}${BASE}_PCA_outliers.txt" | grep -F -v -f "${SYMBIONTS_INSHORE_CN_GBR}" | shuf -n 35 >> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"
cat "${INDIR}${BASE}_filtered2_Group2_noclones.filelist.txt" | grep -F -f "${SYMBIONTS_SUBSET}"  | grep -F -v -f "${INDIR}${BASE}_PCA_outliers.txt" | shuf -n ${N} >> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"
cat "${INDIR}${BASE}_filtered2_Group3_noclones.filelist.txt" | grep -F -f "${SYMBIONTS_SUBSET}"  | grep -F -v -f "${INDIR}${BASE}_PCA_outliers.txt" | shuf -n ${N} >> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"
cat "${INDIR}${BASE}_filtered2_Group4_noclones.filelist.txt" | grep -F -f "${SYMBIONTS_SUBSET}"  | grep -F -v -f "${INDIR}${BASE}_PCA_outliers.txt" | shuf -n ${N} >> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"


cat "${INDIR}${BASE}_filtered2_Group1234_noclones.filelist.txt" | grep 'Group1' | grep -F -v -f "${INDIR}${BASE}_PCA_outliers.txt" | shuf -n ${N} >> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"
cat "${INDIR}${BASE}_filtered2_Group1234_noclones.filelist.txt" | grep 'Group2' | grep -F -v -f "${INDIR}${BASE}_PCA_outliers.txt" | grep -F -v -f "${INDIR}${BASE}_South_GBR_outliers.txt"  | shuf -n ${N} >> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"
cat "${INDIR}${BASE}_filtered2_Group1234_noclones.filelist.txt" | grep 'Group3' | grep -F -v -f "${INDIR}${BASE}_PCA_outliers.txt" | shuf -n ${N} >> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"
cat "${INDIR}${BASE}_filtered2_Group1234_noclones.filelist.txt" | grep 'Group4' | grep -F -v -f "${INDIR}${BASE}_PCA_outliers.txt" | shuf -n ${N} >> "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"

for M in ${MISS[@]};
    do
    
    #Create new vcf file with the subset of samples 
    vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2.vcf.gz"  --keep "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt" --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70.vcf.gz"
    
    #Convert vcf file to PHYLIP format 
    python /home/hdenis/Coral-Genomics/Scripts/Python_scripts/vcf2phylip.py -i "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70.vcf.gz" -o "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70.phy"
    
    #Run RAxML (GTR GAMMA model)
    #-p : random seed
    #-T : number of threads
    #-x : random seed for bootstrap replicates 
    #-N : automatic determination of number of bootstrapping replicates 
    
    #GAMMA model for ML tree (as more reliable)
    raxmlHPC-PTHREADS -s "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70.min4.phy" -n "${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70_ML_tree" -m GTRGAMMA -p 12345 -T 15
    
#    #GTRCAT model for bootstrap
    raxmlHPC-PTHREADS -s "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70.min4.phy" -n "${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70_boot_tree_GAMMA_100" -m GTRGAMMA -p 12345 -x 67890 -# 100 -T 15
#    
#    #Combine results (add bootstrap info to ML tree)
     raxmlHPC-PTHREADS -f b -t RAxML_bestTree.${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70_ML_tree -z RAxML_bootstrap.${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70_boot_tree_GAMMA_100 -n ${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF${MAF}_LD0.2_subsample70_final_tree_GAMMA -m GTRGAMMA -T 15
     
done





