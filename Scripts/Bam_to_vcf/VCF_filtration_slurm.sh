#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/VCF_filtration_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J VCFtools

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=300G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/vcftool_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/vcftool_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=10		# number of cores per job

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

#A.millepora v3 reference genome 
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"
REF_NAME="Amilleporav3"

#Associate slurm array index with contig (chr) name 
CHROMOSOME_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/chromosomes_header.txt"
CONTIG=`sed -n ${SLURM_ARRAY_TASK_ID}p $CHROMOSOME_FILE`


#11.3 Basic site filtering using vcftools

#11.31 Output statistics of raw vcf before filtering to select parameters thresholds 
#Note : vcftools can only apply one function per output so it is required to call the command several time 
#Tune parameter for ld calculation -> might need to add a --thin parameter or subsample a lower number of sites 

###Get number of sites in vcf file 
#NSITE=$(zcat "${INDIR}aspat_clean_all_chr_SNP.vcf.gz"| grep -v '^[#;]' | wc -l)
#echo Number of sites in the vcf file
#echo $NSITE 

#Get number of bi-allelic sites in vcf file 
#bcftools view --no-header -G -m 2 -M 2 --types snps "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" | wc -l

#
### Get number of samples in vcf file 
#echo Number of samples in the vcf file 
#zcat "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" | grep -m1 "^#CHROM" | cut -f 10- | tr "\t" "\n"  | wc -l


##Statistics outputs : 
##Genotype frequencies
##Individual average depths
##Individual average misingness
##Sites average depths
##Sites average misingness
##Nucleotide diversity
##Individuals relatedness
##Sites quality
##Heterozygoty and inbreeding coefficients 
#
#OUTPUT_OPTIONS=("--freq" "--depth" "--site-mean-depth" "--site-pi" "--relatedness" "--site-quality" "--missing-indv" "--missing-site" "--het")
#OUTPUT_OPTIONS=("--missing-site" "--site-quality" "--missing-indv" "--het")
#
#start=`date +%s`
#echo start outputting stats
#for OPTION in "${OUTPUT_OPTIONS[@]}";
#do
#vcftools --gzvcf "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" --out "${INDIR}aspat_clean_all_chr_SNP_nofilt"  --temp /nvme/disk0/lecellier_data/WGS_GBR_data/tmp $OPTION
#end=`date +%s`
#done
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.
#

#Extract additional genotype quality for investigation 
start=`date +%s`
echo start create random subset file 
#
#singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx200g" SelectVariants --variant "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" --output "${INDIR}aspat_clean_all_chr_SNP_random_subset.vcf.gz" --select-random-fraction 0.001
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes. 


start=`date +%s`
echo start outputting GQ table stats
singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx58g" VariantsToTable --show-filtered true -V "${INDIR}aspat_clean_all_chr_SNP_random_subset.vcf.gz" -F CHROM -F POS -F FILTER -GF GQ -O "${INDIR}aspat_clean_all_chr_SNP_random_subset_GQ_stat.table" 
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.  

start=`date +%s`
echo start outputting DP table stats
singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx58g" VariantsToTable --show-filtered true -V "${INDIR}aspat_clean_all_chr_SNP_random_subset.vcf.gz" -F CHROM -F POS -F FILTER -GF DP -O "${INDIR}aspat_clean_all_chr_SNP_random_subset_DP_stat.table" 
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes. 

##Subset frequency file to keep only first two alleles 
#cat aspat_clean_all_chr_SNP_nofilt.frq | awk '{ print $1, $2, $3, $4, $5, $6 }' > aspat_clean_all_chr_SNP_nofilt_biallelic.frq

#Extract additional statistics using bcftools
#start=`date +%s`
#echo start outputting additional stats
#FS score 
#bcftools query "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" -f'%FS\n' > "${INDIR}aspat_clean_all_chr_SNP_nofilt.FS"
##Mapping quality rank sum
#bcftools query "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" -f'%MQRankSum\n' > "${INDIR}aspat_clean_all_chr_SNP_nofilt.MQRankSum"
##ReadPosRankSum
#bcftools query "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" -f'%ReadPosRankSum\n' > "${INDIR}aspat_clean_all_chr_SNP_nofilt.ReadPosRankSum"
##Quality by depth
#bcftools query "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" -f'%QD\n' > "${INDIR}aspat_clean_all_chr_SNP_nofilt.QD"

#Use plink to compute hwe 
#plink --allow-extra-chr --hardy --out "${INDIR}aspat_clean_all_chr_SNP_nofilt" --vcf "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" --threads 4 --const-fid 0 --memory 200000 
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.  

#11.32 First round of filtration using gatk Variant filtration 
#https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
#https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
#QD : Quality by depth : variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.
#QUAL : Variant quality (values below 30 have already been assigned a Lowqual fields by GenotypeGvcfs)
#FS : Fisher strand : Phred-scaled probability that there is strand bias at this site
#MQ : Mapping quality : Root mean square mapping quality over all the reads at the site
#MQRankSum : Mapping quality rank sum test : Compares the mapping qualities of the reads supporting the reference allele and the alternate allele
#ReadPosRankSum : Read position rank sum test : Compares whether the positions of the reference and alternate alleles are different within the reads
#GQ : Genotype quality : Quality of genotype call (for each sample)
#DP : Genotype depth : Unfiltered number of reads supporting genotype call usefull as no call genotypes are marked as 0/0 instead of ./. in gatk version prior to 4.6.0.0 : https://gatk.broadinstitute.org/hc/en-us/articles/6012243429531-GenotypeGVCFs-and-the-death-of-the-dot-obsolete-as-of-GATK-4-6-0-0

#start=`date +%s`
#echo start outputting filtered vcf 
#singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx58g" VariantFiltration \
#    -V "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" \
#    -filter "QD < 2.0" --filter-name "QD2" \
#    -filter "QUAL < 30.0" --filter-name "QUAL30" \
#    -filter "FS > 30.0" --filter-name "FS30" \
#    -filter "MQ < 20.0" --filter-name "MQ20" \
#    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#    -G-filter "GQ < 20.0" -G-filter-name "GQ20" \
#    -G-filter "DP < 1" -G-filter-name "DP1" \
#    -O "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" \
#    --set-filtered-genotype-to-no-call true
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.  

#Extract statistics after filtration
#start=`date +%s`
#echo start outputting filtered vcf stats
#singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx58g" VariantsToTable --show-filtered true -V "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" -F CHROM -F POS -F FILTER -F DP -O "${INDIR}aspat_clean_all_chr_SNP_filtered1_stat.table" 
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.  

##Compute misingness using another tool to confirm (plink) -> remove  
#start=`date +%s`
#echo Start computing missingness
#plink --allow-extra-chr --missing --out "${INDIR}aspat_clean_all_chr_SNP_nofilt_filtered1_plink" --vcf "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" --threads 10 --const-fid 0 --memory 200000  
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.  

#OPTION="--site-mean-depth"
#vcftools --gzvcf "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" --out "${INDIR}aspat_clean_all_chr_SNP_filtered1"  --temp /nvme/disk0/lecellier_data/WGS_GBR_data/tmp $OPTION --remove-filtered-all
#
#OPTION="--missing-indv"
#vcftools --gzvcf "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" --out "${INDIR}aspat_clean_all_chr_SNP_filtered1"  --temp /nvme/disk0/lecellier_data/WGS_GBR_data/tmp $OPTION --remove-filtered-all
#
#OPTION="--missing-site"
#vcftools --gzvcf "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" --out "${INDIR}aspat_clean_all_chr_SNP_filtered1"  --temp /nvme/disk0/lecellier_data/WGS_GBR_data/tmp $OPTION --remove-filtered-all

##11.33 Compute relatedness between samples and filter out clones
#start=`date +%s`
#echo start computing relatedness
#OPTION="--relatedness"
#vcftools --gzvcf "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" --out "${INDIR}aspat_clean_all_chr_SNP_filtered1"  --temp /nvme/disk0/lecellier_data/WGS_GBR_data/tmp $OPTION
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

##Filter clone individuals 
#
##11.34 Additional filtering using vcftools
#
##11.341 Mask SNPs 5bp around indels (optional)
#
##Create bed file from INDEL vcf file to mask SNPs that are +/- 5bp from indel positions
##bcftools query -f '%CHROM\t%POS0\t%END\n' "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" | awk -v s=5 'BEGIN{OFS="\t"}{print $1,$2-s,$3+s}' > "${INDIR}aspat_clean_all_chr_mask_indel.bed"
#
##Exclude SNPs within intervals specified in bed file 
#vcftools --gzvcf "${INDIR}aspat_clean_all_chr_SNP_filtered1.vcf.gz" --exclude-bed "${INDIR}aspat_clean_all_chr_mask_indel.bed" --recode --recode-INFO-all --stdout | gzip -c > "${INDIR}aspat_clean_all_chr_SNP_filtered2.vcf.gz"
#
##11.342 Filter out variants with missing genotypes, that have low depth or are not bi-allelic, MAF>0.05 
##Remove sites with hwe p.value < 0.0001
##Remove sites with >=1 missing genoype
##Remove sites that have <3 mean depth and >13 depth (median +/- sd)
##Select only bi-allelic sites 
#vcftools --gzvcf "${INDIR}aspat_clean_all_chr_SNP_filtered2.vcf.gz" --remove-filtered-all --hwe 0.0001 --max-missing 1 --min-meanDP 3 --max-meanDP 13 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --recode-INFO-all --stdout | gzip -c > "${INDIR}aspat_clean_all_chr_SNP_filtered3.vcf.gz"
#
##11.343 Prune variants by LD with plink -> Generate two files with different levels of pruning 
#plink --vcf "${INDIR}aspat_clean_all_chr_SNP_filtered3.vcf.gz" --allow-extra-chr --indep-pariwise 200 20 0.15 --out "${INDIR}aspat_clean_all_chr_SNP_filtered3_LD_pruned0.15" 
#
#
##--export vcf bgz if necessary to create new vcf file 
#
##Check --set-missing-var-ids if necessary 
#
###11.34 Re-output the same statistics after filtration 

