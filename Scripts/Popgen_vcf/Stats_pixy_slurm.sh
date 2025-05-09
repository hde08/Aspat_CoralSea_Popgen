#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Popgen_vcf/Stats_pixy_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J pixy

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=150G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/pixy_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/pixy_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=1		# number of cores per job

#SBATCH --array=1-10%10      	# job array 1-14%4

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

########################################################## DIVERSITY/DIVERGENCE METRICS ###################################################################
#Scripts to compute population statistics using pixy 1.2.10.beta2 : Fst, Fis, Pi, Dxy
#This method produces statistics unbiased by missing data by using both variant and invariant sites 

INDIR="/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Vcf_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/Pixy/"
OUTDIR2="/nvme/disk0/lecellier_data/GBR_NC_phylo_subset/Vcf_files/"
NEWBASE="Group1234"
BASE="aspat_clean"

#1. Filter invariant sites

#Annotate vcf file per contig to fasten the process
CONTIG="scaffold_${SLURM_ARRAY_TASK_ID}"
start=`date +%s`
echo start annotating vcf 
singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx35g"   VariantFiltration \
    -V "${INDIR}${BASE}_all_chr_NOVARIATION.vcf.gz" \
    -G-filter "DP < 1.0" -G-filter-name "DP1" \
    -G-filter "RGQ < 20.0" -G-filter-name "RGQ20" \
    -L $CONTIG \
    -O "${INDIR}${BASE}_all_chr_NOVARIATION_annotated_${CONTIG}.vcf.gz" \
    --set-filtered-genotype-to-no-call true \
    --verbosity ERROR
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes. 

#Combine VCF files in single VCF 
ls ${INDIR}${BASE}_all_chr_NOVARIATION_annotated_scaffold*.vcf.gz -v > "${INDIR}${BASE}_all_chr_NOVARIATION_annotated_scaffold_list.txt"

singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx150g" GatherVcfs --INPUT  "${INDIR}${BASE}_all_chr_NOVARIATION_annotated_scaffold_list.txt" --OUTPUT "${INDIR}${BASE}_all_chr_NOVARIATION_annotated.vcf.gz" --CREATE_INDEX true 

#Create index file 
singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx260g" IndexFeatureFile -I "${INDIR}${BASE}_all_chr_NOVARIATION_annotated.vcf.gz"

#Remove variants that don't have the PASS filter from vcf file 
start=`date +%s`
echo start outputting filtered vcf 
singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx120g" SelectVariants --variant "${INDIR}${BASE}_all_chr_NOVARIATION_annotated.vcf.gz" --output "${INDIR}${BASE}_all_chr_NOVARIATION_filtered_1.vcf.gz" --exclude-filtered --verbosity ERROR
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes. 

#Correct sample names in header of vcf file 
bcftools reheader -s "${INDIR}${BASE}_all_chr_SNP_filtered_1_correct_sample_names.txt" "${INDIR}${BASE}_all_chr_NOVARIATION_filtered_1.vcf.gz" > "${INDIR}${BASE}_all_chr_NOVARIATION_filtered_1_reheader.vcf.gz"
mv "${INDIR}${BASE}_all_chr_NOVARIATION_filtered_1_reheader.vcf.gz" "${INDIR}${BASE}_all_chr_NOVARIATION_filtered_1.vcf.gz"


#Filter invariant file, we restrict the dataset to 70 samples per genomic cluster
KEEP_SAMPLES="${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70.filelist.txt"

#Filters: 
#RGQ > 20.0
#--max-missing 0.8
#--min-meanDP 5
#--max-meanDP 50
#--max-maf 0

vcftools --gzvcf "${INDIR}${BASE}_all_chr_NOVARIATION_filtered_1.vcf.gz" --keep "${KEEP_SAMPLES}" --remove-filtered-all --max-missing 0.8 --min-meanDP 5 --max-meanDP 50 --max-maf 0 --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > "${OUTDIR2}${NEWBASE}_all_chr_NOVARIATION_filtered_2_20mis_subsample70.vcf.gz" 

tabix "${OUTDIR2}${NEWBASE}_all_chr_NOVARIATION_filtered_2_20mis_subsample70.vcf.gz" 


#2. Filter SNPs following recommandations of Hemstrom 2024 and Sopniewski 2024

#Filters:
#--min-meanDP 15
#--max-meanDP 50
#--min-alleles 2
#retain multi-allelic sites 
#--mac 1
#--hwe 0.0001
#--max-missing 0.8
#GQ>20

vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss.vcf.gz" --keep "${KEEP_SAMPLES}" --remove-filtered-all --max-missing 0.8 --min-meanDP 5 --max-meanDP 50 --mac 1 --hwe 0.0001 --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > "${OUTDIR2}${NEWBASE}_all_chr_SNP_filtered_2_20mis_mac1_hwe0.0001_subsample70.vcf.gz"  

tabix "${OUTDIR2}${NEWBASE}_all_chr_SNP_filtered_2_20mis_mac1_hwe0.0001_subsample70.vcf.gz" 


#3. Concatenate invariant and variant vcfs
bcftools concat --allow-overlaps "${OUTDIR2}${NEWBASE}_all_chr_SNP_filtered_2_20mis_mac1_hwe0.0001_subsample70.vcf.gz" "${OUTDIR2}${NEWBASE}_all_chr_NOVARIATION_filtered_2_20mis_subsample70.vcf.gz"  -O z -o "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample70.vcf.gz"

tabix "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample70.vcf.gz"

paste <(cat "$KEEP_SAMPLES") <(printf 'Group1\n%.0s' {1..67}; printf 'Group2\n%.0s' {1..67}; printf 'Group3\n%.0s' {1..70}; printf 'Group4\n%.0s' {1..68}) > "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70_popfile.txt"


#4. Run pixy : compute metrics using genomic clusters as populations 

#Parameters:
#Window size = 10kb 
#N samples per genomic cluster = 2-70
#10 random samples replicates 

source activate base
conda activate pixy

N=(2 5 10 20 30 40 50 70)
K=${SLURM_ARRAY_TASK_ID}

#echo -e "Replicate\tN_SAMPLES\tN_SPNs" > /nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/Pixy/v7/summary_pixy_replicates.txt

for n in "${N[@]}"; do
    
    #Select a random set of n individuals from each group and filter for missingness
    for g in 1 2 3 4; do
      cat "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70_popfile.txt" | grep "Group${g}" | shuf -n "$n" >> "${OUTDIR2}${BASE}_temp_popfile_${n}_${K}.txt"
    done
    
    awk '{print $1}' "${OUTDIR2}${BASE}_temp_popfile_${n}_${K}.txt" > "${OUTDIR2}${BASE}_temp_${n}_${K}.filelist.txt"
    
    #Filter the vcf file t 
    plink \
  --vcf "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample70.vcf.gz" \
  --allow-extra-chr \
  --geno 0.2 \
  --keep "${OUTDIR2}${BASE}_temp_${n}_${K}.filelist.txt" \
  --out "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample${n}_${K}" \
  --const-fid \
  --keep-allele-order \
  --recode vcf-iid bgz \
  --set-missing-var-ids @:#[Aspat]
  
    vcftools --gzvcf "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample70.vcf.gz" --keep "${OUTDIR2}${BASE}_temp_${n}_${K}.filelist.txt" --remove-filtered-all --max-missing 0.8 --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample${n}_${K}.vcf.gz"  

    tabix "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample${n}_${K}.vcf.gz"
    
    #Save number of individuals and SNPs in a folder 
    N_SAMPLES=$n
    N_SNPS=$(zcat "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample${n}_${K}.vcf.gz" | grep -v '^#' | wc -l)
    
    echo -e "$K\t$n\t$N_SNPS" >> /nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/Pixy/v7/summary_pixy_replicates.txt

    #Use it to compute values with pixy 
    pixy --stats pi fst dxy --vcf "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample${n}_${K}.vcf.gz"  --populations "${OUTDIR2}${BASE}_temp_popfile_${n}_${K}.txt" --window_size 10000 --n_cores 1 --output_folder /nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/Pixy/v7/ --output_prefix Group1234_${n}_${K}
    
    #remove temporary files 
    rm "${OUTDIR2}${BASE}_temp_popfile_${n}_${K}.txt" 
    rm "${OUTDIR2}${BASE}_temp_${n}_${K}.filelist.txt"
    rm "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample${n}_${K}.vcf.gz"
    rm "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample${n}_${K}.vcf.gz.tbi"
    rm "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample${n}_${K}.nosex"
    rm "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample${n}_${K}.log"
    
done


cd /nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/Pixy/v5/

pixy --stats pi fst dxy --vcf "${OUTDIR2}${NEWBASE}_all_chr_allsites_filtered_v3_subsample70.vcf.gz"  --populations "${INDIR}${BASE}_filtered2_Group1234_noclones_subsample70_popfile.txt" --window_size 10000 --n_cores 4



#4. Run pixy : compute metrics at the reef level (5 individuals per reef, 10 replicates per reef)

#Compute pi per reef 
n=5
K=${SLURM_ARRAY_TASK_ID}

#4.1 Create a random subset of the dataset of 5 samples per reef

#Find unique reefs
awk '{print $2}' "${INDIR}${BASE}_Group1234_noclones_nooutliers_popfile_perreef.txt" | sort -u > "${INDIR}${BASE}_reef_list.txt"
FITZ_N=1

> "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}_popfile.txt"
while IFS= read -r line; do
  REEF="$line"
  echo $REEF
  
  #There is an issue coming from the fact Fitzroy Reef and Fitzroy Island have the same abbreviation 
  #Fitzroy Reef : wgsID range between 1656 and 1716
  #Fitzroy Island : wgsID range between 772 and 846
  if [ "$REEF" == "FITZ" ]; then
    #First FITZ is used to pick individuals from Fitzroy Reef
    if [ "$FITZ_N" == "1" ]; then
      count=$(cat "${INDIR}${BASE}_Group1234_noclones_nooutliers_popfile_perreef.txt" | grep "$REEF" | grep -E '165[6-9]|16[6-9][0-9]|170[0-9]|171[0-6]'  | shuf -n "$n" | wc -l)
      
      if [ "$count" -ge 5 ]; then
        echo "Computing pi for this reef"
        cat "${INDIR}${BASE}_Group1234_noclones_nooutliers_popfile_perreef.txt" | grep "$REEF" | grep -E '165[6-9]|16[6-9][0-9]|170[0-9]|171[0-6]' | shuf -n "$n" | awk '{$2 = "FITZ_R"; print}' >> "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}_popfile.txt"
        FITZ_N=2
    
      fi
    
    else
      count=$(cat "${INDIR}${BASE}_Group1234_noclones_nooutliers_popfile_perreef.txt" | grep "$REEF" | grep -E '77[2-9]|7[8-9][0-9]|8[0-3][0-9]|84[0-6]'  | shuf -n "$n" | wc -l)
      
      if [ "$count" -ge 5 ]; then
    echo "Computing pi for this reef"
    cat "${INDIR}${BASE}_Group1234_noclones_nooutliers_popfile_perreef.txt" | grep "$REEF" | grep -E '77[2-9]|7[8-9][0-9]|8[0-3][0-9]|84[0-6]' | shuf -n "$n" | awk '{$2 = "FITZ_I"; print}' >> "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}_popfile.txt"
      fi
    fi
  else

    # Assign the result of the command to a variable
  count=$(cat "${INDIR}${BASE}_Group1234_noclones_nooutliers_popfile_perreef.txt" | grep "$REEF" | shuf -n "$n" | wc -l)

 # Test if the count is greater than or equal to 5
   if [ "$count" -ge 5 ]; then
     echo "Computing pi for this reef"
     cat "${INDIR}${BASE}_Group1234_noclones_nooutliers_popfile_perreef.txt" | grep "$REEF" | shuf -n "$n" >> "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}_popfile.txt"
   
   fi
  fi
done < "${INDIR}${BASE}_reef_list.txt"

awk '{print $1}' "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}_popfile.txt" > "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}.filelist.txt"


#4.2 Use this subset for the filtration and computation 
KEEP_SAMPLES="${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}.filelist.txt"
   
    
#Filter invariant sites 
vcftools --gzvcf "${INDIR}${BASE}_all_chr_NOVARIATION_filtered_1.vcf.gz" --keep "${KEEP_SAMPLES}" --remove-filtered-all --max-missing 0.8 --min-meanDP 5 --max-meanDP 50 --max-maf 0 --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > "${INDIR}${NEWBASE}_all_chr_NOVARIATION_filtered_2_20mis_5perreef_${K}.vcf.gz" 

#Filter variant sites 
vcftools --gzvcf "${INDIR}${BASE}_all_chr_SNP_filtered_indivmiss.vcf.gz" --keep "${KEEP_SAMPLES}" --remove-filtered-all --max-missing 0.8 --min-meanDP 5 --max-meanDP 50 --mac 1 --hwe 0.0001 --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > "${INDIR}${NEWBASE}_all_chr_SNP_filtered_2_20mis_mac1_hwe0.0001_5perreef_${K}.vcf.gz"  

tabix "${INDIR}${NEWBASE}_all_chr_NOVARIATION_filtered_2_20mis_5perreef_${K}.vcf.gz" 
tabix "${INDIR}${NEWBASE}_all_chr_SNP_filtered_2_20mis_mac1_hwe0.0001_5perreef_${K}.vcf.gz" 

#Concatenate files 
bcftools concat --allow-overlaps "${INDIR}${NEWBASE}_all_chr_SNP_filtered_2_20mis_mac1_hwe0.0001_5perreef_${K}.vcf.gz" "${INDIR}${NEWBASE}_all_chr_NOVARIATION_filtered_2_20mis_5perreef_${K}.vcf.gz"  -O z -o "${INDIR}${NEWBASE}_all_chr_allsites_filtered_v3_5perreef_${K}.vcf.gz"

tabix "${INDIR}${NEWBASE}_all_chr_allsites_filtered_v3_5perreef_${K}.vcf.gz"

cd /nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/Pixy/v9/

#Create popfile for pixty 
awk '{print $1 "\t" $2}' "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}_popfile.txt" > tmp_${K}.txt
mv tmp_${K}.txt "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}_popfile.txt"

#Run pixy 
pixy --stats pi fst dxy --vcf "${INDIR}${NEWBASE}_all_chr_allsites_filtered_v3_5perreef_${K}.vcf.gz"  --populations "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}_popfile.txt" --window_size 10000 --n_cores 1 --output_folder /nvme/disk0/lecellier_data/WGS_GBR_NC_data/Popgen_outputs/Pixy/v9/ --output_prefix Group1234_${n}perreef_${K}

#Remove temp files 
rm "${INDIR}${NEWBASE}_all_chr_NOVARIATION_filtered_2_20mis_5perreef_${K}.vcf.gz" 
rm "${INDIR}${NEWBASE}_all_chr_SNP_filtered_2_20mis_mac1_hwe0.0001_5perreef_${K}.vcf.gz"  
rm "${INDIR}${NEWBASE}_all_chr_allsites_filtered_v3_5perreef_${K}.vcf.gz"
rm "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}_popfile.txt"
rm "${INDIR}${BASE}_Group1234_noclones_nooutliers_5_perreef_${K}.filelist.txt"





