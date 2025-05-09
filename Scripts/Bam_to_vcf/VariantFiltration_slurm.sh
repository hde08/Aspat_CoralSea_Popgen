#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Bam_to_vcf/VariantFiltration_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH 

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=50G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/varfilt_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/varfilt_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=1		# number of cores per job

#SBATCH --array=1-4%4        	# job array

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

########################################################## VCF FILE FILTERING ###################################################################
#This scripts serves to filter vcf files outputted by the GATK pipeline prior conducting population genomic analyses 
#Other specific filtrations for different analyses (e.g. demographic inference in dadi were done in separate scripts provided along these analyses scripts)


#1. Convert genotypes with bad quality (<20) to missing data and additional hard-filter of variants as recommanded in GATK guidelines

#https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
#https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants


cd /nvme/disk0/lecellier_data/WGS_GBR_NC_data/Vcf_files/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Vcf_files/"
BASE="aspat_clean"

REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"
REF_NAME="Amilleporav3"

start=`date +%s`
echo start annotating vcf 
singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx260g"   VariantFiltration \
    -V "${INDIR}${BASE}_all_chr_SNP.vcf.gz" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "FS > 50.0" --filter-name "FS50" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -G-filter "DP < 1.0" -G-filter-name "DP1" \
    -G-filter "GQ < 20.0" -G-filter-name "GQ20" \
    -O "${INDIR}${BASE}_all_chr_SNP_annotated.vcf.gz" \
    --set-filtered-genotype-to-no-call true \
    --verbosity ERROR
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes. 

#Remove variants that don't have the PASS filter from vcf file 
#start=`date +%s`
echo start outputting filtered vcf 
singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx260g" SelectVariants --variant "${INDIR}${BASE}_all_chr_SNP_annotated.vcf.gz" --output "${INDIR}${BASE}_all_chr_SNP_filtered_1.vcf.gz" --exclude-filtered --verbosity ERROR
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes. 

##Correct RRAP sample names in file 
#Some samples labels were mispecificied in the original fasta file and were corrected at this stage to ensure compatible binding with metadata in downstream analyses 

#Get previous sample names
bcftools query -l "${INDIR}${BASE}_all_chr_SNP_filtered_1.vcf.gz" > "${INDIR}${BASE}_all_chr_SNP_filtered_1_old_sample_names.txt"
cp "${INDIR}${BASE}_all_chr_SNP_filtered_1_old_sample_names.txt" "${INDIR}${BASE}_all_chr_SNP_filtered_1_correct_sample_names.txt"
#Correct manually samples names in text editor using Correct_sample_names.xlsx file
#Correct file header
bcftools reheader -s "${INDIR}${BASE}_all_chr_SNP_filtered_1_correct_sample_names.txt" "${INDIR}${BASE}_all_chr_SNP_filtered_1.vcf.gz" > "${INDIR}${BASE}_all_chr_SNP_filtered_1_reheader.vcf.gz"

mv "${INDIR}${BASE}_all_chr_SNP_filtered_1_reheader.vcf.gz" "${INDIR}${BASE}_all_chr_SNP_filtered_1.vcf.gz"



#2. Compute individual statistics on filtered VCF to remove individuals with very low missingness, quality or depth 
#Statistics outputs : 
#Individual average depths
#Individual average misingness
#Percentage of multi-allelic SNPs

OUTPUT_OPTIONS=("--depth" "--missing-indv")
start=`date +%s`
echo start outputting stats
for OPTION in "${OUTPUT_OPTIONS[@]}";
do
vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_1.vcf.gz" --out "${INDIR}${BASE}_all_chr_SNP_filtered_1"  --temp /nvme/disk0/lecellier_data/WGS_GBR_data/tmp $OPTION
end=`date +%s`
done
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Compute percentage of multiallelic SNPs in filtered VCF 
bcftools view --no-header -G -m 2 -M 2 --types snps "${INDIR}${BASE}_all_chr_SNP_filtered_1.vcf.gz" | wc -l > "${INDIR}${BASE}_all_chr_SNP_filtered_1_perc_biallelic.txt"



#3. Remove individuals with high missingness and outlier samples 
#The distribution of missing data among individuals is visualized using using R script Vcf_statistics.Rmd and used to create files below 
keep_indiv="${INDIR}${BASE}_keep_individuals.txt"
rm_indiv="${INDIR}${BASE}_remove_individuals.txt"

vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_1.vcf.gz" --remove "${INDIR}${BASE}_remove_individuals.txt" --remove-filtered-all  --recode --recode-INFO-all --stdout | bgzip -c > "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss.vcf.gz"



#4. Compute site statistics after they have been removed 
#Statistics plot are visualized using script Vcf_statistics.Rmd
#The vcf file is thinned to 1% to increase the speed of computation 

#Create index 
gunzip -c "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss.vcf.gz" | bgzip -c > "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss.vcf.bgz"
singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx260g" IndexFeatureFile -I "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss.vcf.gz"

#Thin vcf file 
singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx200g" SelectVariants --variant "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss.vcf.gz" --output "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss_thinned0.01.vcf.gz" --select-random-fraction 0.01

#Compute different site statistics on the thinned file 
OUTPUT_OPTIONS=("--site-mean-depth" "--site-pi" "--site-quality" "--missing-site" "--freq" "--hardy")
start=`date +%s`
echo start outputting stats
for OPTION in "${OUTPUT_OPTIONS[@]}";
do
vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss_thinned0.01.vcf.gz" --out "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss_thinned0.01"  --temp /nvme/disk0/lecellier_data/WGS_GBR_data/tmp $OPTION
end=`date +%s`
done
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.


#5. General filtration
#Different datasets with different amount of missing data are generated following Hemstrom et al., 2023
#Consistency of population genomic results (such as PCA or ADMIXTURE) are checked across these datasets
MISS=(05 10 20 50)
for M in ${MISS[@]};
    do
    MAX_MISS=$(echo "scale=2; (100 - $M) / 100" | bc)
    vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_1.vcf.gz" --remove "${INDIR}${BASE}_remove_individuals.txt"  --remove-filtered-all --mac 1 --max-missing $MAX_MISS --min-meanDP 5 --max-meanDP 15 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > "${INDIR}${BASE}_all_chr_SNP_filtered_2_${MISS}mis.vcf.gz"
done



#6. Confirm the absence of batch effects
#see Script Compare_sequencing_batch.sh


#7. Compute samples relatedness (method of Yang et al, Nature Genetics 2010 (doi:10.1038/ng.608))
#This is done to confirm clones previously identified in ANGSD
#Clones are identified separately for each genomic cluster using the 20% missing dataset 
POP_FILES=(${INDIR}${BASE}_filtered2_Group?.filelist.txt)

for POP in ${POP_FILES[@]};
    do
        GROUP=$(basename $POP)
        GROUP=${GROUP%%.filelist*}
        GROUP=${GROUP#*2_}
        vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis.vcf.gz" --keep "${POP}" --remove-filtered-all  --relatedness --out "/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Clones_ID/${BASE}_all_chr_SNP_filtered_2_20mis_${GROUP}"      
done


#8. Filtration for Population genetic analyses (LD, MAF) ("Dataset 2")
#Technical replicates and one clone from each pairs of natural clones are removed 
#MAC>1, MAF>0.05, LD>0.2 (200 20 0.2)
#Test for missingness varying between 5 and 50%
#MISS=(05 10 20 50)
MISS=(20)

for M in ${MISS[@]};
    do
      MAX_MISS=$(echo "scale=2; (100 - $M) / 100" | bc)
      #Remove clones / Tech replicates and re-apply filters 
      vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis.vcf.gz" --remove "${INDIR}${BASE}_all_chr_SNP_filtered_2_clones_technical_replicates.txt" --remove-filtered-all --mac 1 --max-missing ${MAX_MISS} --min-meanDP 5 --max-meanDP 15 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones.vcf.gz"
      #In ADMIXTURE chromosomes are coded as integers thus the vcf file needs to be recoded first and then plink ran
      zcat "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones.vcf.gz" |\
      	sed s/^scaffold_//g |\
      	bgzip > "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded.vcf.gz"
      
      #Filter for LD and MAF using plink
      plink --vcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded.vcf.gz" --allow-extra-chr --maf 0.05 --indep-pairwise 200 20 0.2 --out "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF0.05_LD0.2" --const-fid --keep-allele-order --set-missing-var-ids @:#[Aspat]
      
      plink --vcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded.vcf.gz" --allow-extra-chr --extract "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF0.05_LD0.2.prune.in" --out "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF0.05_LD0.2" --const-fid --keep-allele-order --recode vcf-iid bgz --set-missing-var-ids @:#[Aspat] 
       
      #Output in bed format for ADMIXTURE 
      plink --vcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded.vcf.gz" --allow-extra-chr --extract "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF0.05_LD0.2.prune.in" --out "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF0.05_LD0.2" --const-fid --keep-allele-order --make-bed --set-missing-var-ids @:#[Aspat] 
      
      #Same filtering applied without LD pruning for comparison 
      #MAF filtering without LD filtering for comparison
      plink --vcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded.vcf.gz" --allow-extra-chr --maf 0.05 --out "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded_MAF0.05" --const-fid --keep-allele-order --recode vcf-iid bgz --set-missing-var-ids @:#[Aspat]
      
      #Remove intermediary files
      rm "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones.vcf.gz"
      rm "${INDIR}${BASE}_all_chr_SNP_filtered_2_${M}mis_noclones_recoded.vcf.gz"
done


#9. Compute heterosigosity for each cluster 
#POP_FILES=(${INDIR}${BASE}_filtered2_Group?.filelist.txt)

#for POP in ${POP_FILES[@]};
#    do
#        GROUP=$(basename $POP)
#        GROUP=${GROUP%%.filelist*}
#        GROUP=${GROUP#*2_}
#        vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis.vcf.gz" --keep "${POP}" --remove-filtered-all  --het --out "/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/${BASE}_all_chr_SNP_filtered_2_20mis_${GROUP}"      
#done


#10. Add filtering information to final VCF file header provided along the publication 


#Get number of individuals in VCF file (1128)
#cd /nvme/disk0/lecellier_data/WGS_GBR_NC_data/Vcf_files/
#zcat aspat_clean_all_chr_SNP_filtered_1.vcf.gz | grep -m1 "^#CHROM" | cut -f 10- | tr "\t" "\n"  | wc -l

#Get number of SNPs in VCF file () 
#zcat aspat_clean_all_chr_SNP_filtered_1.vcf.gz | grep -v '^[#;]' | wc -l

#zcat "${INDIR}${BASE}_all_chr_SNP_filtered_2_20mis_reseqindividuals.vcf.gz" | grep -v '^[#;]' | wc -l
