#!/bin/bash
#sbatch /home/hdenis/Coral-Genomics/Scripts/Popgen_vcf/Compare_sequencing_batch.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J COMPSEQ

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=100G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/compseq_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/compseq_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=5		# number of cores per job

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

### Ressources
#module load intel/2022.0.2
#module load openmpi/1.10.7
#module load netcdf-fortran-4.5.3-intel-2021.5.0-gzuyinc
#module load netcdf-c-4.8.1-intel-2021.5.0-j7rukz5

# This job's working directory :
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

ulimit -s unlimited

#########################################################
###Check for batch effects using hard called variants ###
#########################################################

cd /nvme/disk0/lecellier_data/WGS_GBR_NC_data/Vcf_files/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Vcf_files/"
OUTDIR="/nvme/disk0/lecellier_data/Compare_sequencing_batch/"
BASE="aspat_clean"

#1. Extract resequenced individuals in separate vcf file (N = 21 samples x2)
RESEQ_INDIV="/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Vcf_files/aspat_resequenced_samples.filelist.txt"

vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis.vcf.gz" --keep "${RESEQ_INDIV}" --remove-filtered-all  --recode --recode-INFO-all --stdout | gzip -c > "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals.vcf.gz"

#Remove monomorphic sites in this subset using MAF>0.01 (retains 130344 sites)
plink --vcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals.vcf.gz" --allow-extra-chr --maf 0.01  --out "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals_MAF0.01" --const-fid --keep-allele-order --recode vcf-iid bgz --set-missing-var-ids @:#[Aspat]

#LD-pruning to conduct PCA analysis (R2<0.2, retains 9664 sites) 
plink --vcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals_MAF0.01.vcf.gz" --allow-extra-chr --indep-pairwise 200 20 0.2 --out "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals_MAF0.01_LD0.2" --const-fid --keep-allele-order --set-missing-var-ids @:#[Aspat]
      
plink --vcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals_MAF0.01.vcf.gz" --allow-extra-chr --extract "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals_MAF0.01_LD0.2.prune.in" --out "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals_MAF0.01_LD0.2" --const-fid --keep-allele-order --recode vcf-iid bgz --set-missing-var-ids @:#[Aspat] 


#2. Compute samples relatedness
#For this step we generate another vcf file with resequenced individuals plus 100 random other individuals for comparison. And compute relatedness (method of Yang et al, Nature Genetics 2010 (doi:10.1038/ng.608))
cat "${RESEQ_INDIV}" "${INDIR}aspat_random100.filelist.txt" > "${INDIR}aspat_random100_resequenced.filelist.txt"

vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis.vcf.gz" --keep "${INDIR}aspat_random100_resequenced.filelist.txt" --remove-filtered-all  --recode --recode-INFO-all --stdout | gzip -c > "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals_100random.vcf.gz"

OPTION="--relatedness"
vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals_100random.vcf.gz" --out "${OUTDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals_100random"  --temp /nvme/disk0/lecellier_data/WGS_GBR_data/tmp $OPTION

#3. Compute samples heterosigosity 
#This step is done in two several steps :

#First we compute SNP heterosigosity using no missing data, MAF>0.01 (as recommended per Schmidt 2021 DOI: 10.1111/2041-210X.13659).
MISS=(0)
for M in ${MISS[@]};
  do
  MAX_MISS=$(echo "scale=2; (100 - $M) / 100" | bc)
  
  vcftools --gzvcf "${OUTDIR}${BASE}_all_chr_SNP_filtered2_reseq_indiv_${M}mis.vcf.gz"  --maf 0.01 --het --out "${OUTDIR}Stats/${BASE}_all_chr_SNP_filtered2_reseq_indiv_${M}mis"
  
  
done

#Then we compute the number of unvariant sites for these individuals and we sum it with the previous denominator to compute autosomal heterozygosity which is more accurate. 
MISS=(0 05 10)
for M in ${MISS[@]};
  do
  MAX_MISS=$(echo "scale=2; (100 - $M) / 100" | bc)
  
  #Get number of sites that are all homozigous reference (those sites are outputed in the NOVARIATION type of GATK SelectVariants)
  vcftools --gzvcf "${INDIR}${BASE}_all_chr_NOVARIATION_filtered_2.vcf.gz" --keep "${RESEQ_INDIV}" --remove-filtered-all --max-missing $MAX_MISS --min-meanDP 5 --max-meanDP 15 --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > "/nvme/disk0/lecellier_data/Compare_sequencing_batch/${BASE}_all_chr_NOVARIATION_reseq_indiv_${M}mis.vcf.gz" 
  
  bcftools query -f '[%SAMPLE\t%GT\n]' "/nvme/disk0/lecellier_data/Compare_sequencing_batch/${BASE}_all_chr_NOVARIATION_reseq_indiv_${M}mis.vcf.gz" | \
awk '$2 != "./." {count[$1]++} END {for (s in count) print s "\t" count[s]}' > "${OUTDIR2}${BASE}_all_chr_NOVARIATION_reseq_indiv_${M}mis_counts.txt"
  
  #Get number of sites that are all heterozigous or homozigous mutant among samples (those sites are considered variant as they difer from the reference allele, although monomorphic for our datasets thus they are outputed in the SNP type of GATK SelectVariants)
  vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss.vcf.gz" --keep "${RESEQ_INDIV}" --exclude-bed "${INDIR}${BASE}_all_chr_mask_indel.bed" --max-missing $MAX_MISS --min-meanDP 5 --max-meanDP 15 --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > "/nvme/disk0/lecellier_data/Compare_sequencing_batch/${BASE}_all_chr_SNP_filtered2_reseq_indiv_${M}mis.vcf.gz" 
  
  bcftools view -i '((COUNT(GT="0/0") + COUNT(GT="./.")) == N_SAMPLES) || ((COUNT(GT="1/1") + COUNT(GT="./.")) == N_SAMPLES) || ((COUNT(GT="0/1") + COUNT(GT="./.")) == N_SAMPLES)' "/nvme/disk0/lecellier_data/Compare_sequencing_batch/${BASE}_all_chr_SNP_filtered2_reseq_indiv_${M}mis.vcf.gz" \
    -Oz -o "/nvme/disk0/lecellier_data/Compare_sequencing_batch/${BASE}_all_chr_SNP_filtered2_reseq_indiv_${M}mis_invariants.vcf.gz" 
    
  bcftools query -f '[%SAMPLE\t%GT\n]' "/nvme/disk0/lecellier_data/Compare_sequencing_batch/${BASE}_all_chr_SNP_filtered2_reseq_indiv_${M}mis_invariants.vcf.gz" | \
awk '$2 != "./." {count[$1]++} END {for (s in count) print s "\t" count[s]}' > "${OUTDIR2}${BASE}_all_chr_SNP_filtered2_reseq_indiv_${M}mis_counts.txt"
  
  #Sum both number to get the total number of invariant sites this set of samples
  #This is added to the denominatory to compute autosomal observed heterosigosity 
  awk 'NR==FNR {a[$1]=$2; next} 
     $1 in a {print $1, a[$1] + $2}' "${OUTDIR2}${BASE}_all_chr_NOVARIATION_reseq_indiv_${M}mis_counts.txt" "${OUTDIR2}${BASE}_all_chr_SNP_filtered2_reseq_indiv_${M}mis_counts.txt" > "${OUTDIR2}${BASE}_all_chr_reseq_indiv_${M}mis_unvariant_sites_counts.txt"
    
done

###NEED TO ADD HERE LINES TO COMPUTE HETEROSIGOSITY from the het script 