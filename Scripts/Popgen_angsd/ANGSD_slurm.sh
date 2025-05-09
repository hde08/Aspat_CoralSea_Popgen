#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/ANGSD_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J ANGSD

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=200G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/angsd_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/angsd_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=5		# number of cores per job

#SBATCH --array=1-14%2        	# job array

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

########################################################## SNP CALLING USING GL ###################################################################
#This scripts perform SNP and genotype calling using genotype likelihoods in ANGSD v0.1.17

#A.millepora v3 reference genome 
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"
REF_NAME="Amilleporav3"
#Create samtools ID
samtools faidx $REF_3

#1. Create separate beagle file for GBR and NC data

#1.1 GBR
cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/"
STATDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/BAM_statistics/"

#1.1.1 Create filelist with full path of each BAM file (one filepath per line)
FILES=($INDIR*MARKED_DUP.bam)
printf "%s\n" "${FILES[@]}" > $INDIR/bam.filelist.txt

#Remove one file with low number of reads (<1 M : Aspat-NDIR-489)
#cat $INDIR/bam.filelist.txt | grep 'Aspat-NDIR-489' -v > $INDIR/bam.filelist2.txt

readarray -t FILES <  $INDIR/aspat_bam_clean.filelist.txt

#Get number of files
N_FILES="${#FILES[@]}"
MIN_N=$((95*N_FILES/100))
DEPTH_FILES="${FILES[@]%%_MARKED*}"

#1.1.2 Get median and sd depth value across bam files 
#Store genome-wide mean depth for each fle 
DEPTH_FILES=($STATDIR*Amillepora*coverage.txt)
> $STATDIR/bamfile_depth.txt
> $STATDIR/bamfile_sd_depth.txt

for FILE in ${DEPTH_FILES[@]}; do
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

#1.1.3 Infer SNPs and compute genotype likelihoods 

#Filters :
#--minQ 30 base quality >30
#--minMapQ mapping quality >30
#--minInd MIN_N sites with data from at least 95% of individuals
#--setMinDepthInd 3 minimum of 3 reads for each individual 
#--uniqueOnly 0 keep reads with multiple hits
#--only_proper_pairs 1 only paired reads (already filtered)
#--remove_bads 1 (remove reads flags with 256)

FILTERS="-uniqueOnly 0 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd $MIN_N -setMinDepthInd 3 -SNP_pval 1e-6 -setMaxDepthInd $MAX_DEPTH"

#Parameters
#--GL 1 (SAMtools model)
#--doGlf 2 (output in beagle format)
#--doMaf 2 (fixed major, unkown minor, consider switching to 4 )
#---doMajorMinor 1 (inferred directly from genotype likelihoods)
#--SNP_pval 1e-6 (variant likelihood ratio <0.0000001)

TODO="-doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1"

#Keep only SNPs that matched on the first 14 scaffolds (chromosomes)
cat $REF_3 | grep -m 14  '>' | cut -c 2- > $CHROMOSOME_FILE
#Associate slurm array index with contig (chr) name 
CHROMOSOME_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/chromosomes_header.txt"
CONTIG=`sed -n ${SLURM_ARRAY_TASK_ID}p $CHROMOSOME_FILE`

#Run ANGSD separately by chromosome
start=`date +%s`
angsd -bam $INDIR/aspat_bam_clean.filelist.txt -out "${OUTDIR}GBR_allsamples_filt_05mis_uniq_${CONTIG}" -ref $REF_3 -r $CONTIG: -nThreads 5 $TODO $FILTERS
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Merge beagle files together at the end 
cat <(zcat $INDIR/GBR_allsamples_filt_05mis_uniq_scaffold_1.beagle.gz | head -n 1) <(zcat $INDIR/GBR_allsamples_filt_05mis_uniq_scaffold_*.beagle.gz | grep -v -w marker) | gzip > $INDIR/GBR_allsamples_filt_05mis_uniq_all_chr.beagle.gz
rm $INDIR/GBR_aspatsamples_filt_05mis_uniq_scaffold_*.beagle.gz

#Note on beagle file coding
#0=A, 1=C, 2=G, 3=T



#1.2 NC
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/ANGSD_files/"
STATDIR="/nvme/disk0/lecellier_data/WGS_NC_data/BAM_statistics/"
mkdir "${OUTDIR}"
mkdir "${STATDIR}"

#1.2.1Create filelist with full path of each BAM file (one filepath per line)
FILES=($INDIR*MARKED_DUP.bam)
printf "%s\n" "${FILES[@]}" > $INDIR/bam.filelist.txt
readarray -t FILES <  $INDIR/aspat_bam_clean.filelist.txt

#Get number of files
N_FILES="${#FILES[@]}"
MIN_N=$((95*N_FILES/100))

DEPTH_FILES="${FILES[@]%%_MARKED*}"

#1.2.2 Get median and sd depth value across bam files 
#Store genome-wide mean depth for each fle 
#DEPTH_FILES=($STATDIR*Amillepora*coverage.txt)
> $STATDIR/bamfile_depth.txt
> $STATDIR/bamfile_sd_depth.txt

for FILE in ${DEPTH_FILES[@]}; do
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

#1.2.3 Infer SNPs and compute genotype likelihoods 

#Filters :
#--minQ 30 base quality >30
#--minMapQ mapping quality >30
#--minInd MIN_N sites with data from at least 95% of individuals
#--setMinDepthInd 3 minimum of 3 reads for each individual 
#--uniqueOnly 0 keep reads with multiple hits
#--only_proper_pairs 1 only paired reads (already filtered)
#--remove_bads 1 (remove reads flags with 256)

FILTERS="-uniqueOnly 0 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd $MIN_N -setMinDepthInd 3 -SNP_pval 1e-6 -setMaxDepthInd $MAX_DEPTH"

#Parameters
#--GL 1 (SAMtools model)
#--doGlf 2 (output in beagle format)
#--doMaf 2 (fixed major, unkown minor, consider switching to 4 )
#---doMajorMinor 1 (inferred directly from genotype likelihoods)
#--SNP_pval 1e-6 (variant likelihood ratio <0.0000001)

TODO="-doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1"


#Keep only SNPs that matched on the first 14 scaffolds (chromosomes)
cat $REF_3 | grep -m 14  '>' | cut -c 2- > $CHROMOSOME_FILE

#Associate slurm array index with contig (chr) name 
CHROMOSOME_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/chromosomes_header.txt"
CONTIG=`sed -n ${SLURM_ARRAY_TASK_ID}p $CHROMOSOME_FILE`

#Run ANGSD separately by chromosome
start=`date +%s`
angsd -bam $INDIR/bam.filelist.txt -out "${OUTDIR}NC_allsamples_filt_05mis_uniq_${CONTIG}" -ref $REF_3 -r $CONTIG: -nThreads 5 $TODO $FILTERS
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.


#Merge beagle files together at the end 
cat <(zcat $INDIR/NC_allsamples_filt_05mis_uniq_scaffold_1.beagle.gz | head -n 1) <(zcat $INDIR/NC_allsamples_filt_05mis_uniq_scaffold_*.beagle.gz | grep -v -w marker) | gzip > $INDIR/NC_allsamples_filt_05mis_uniq_all_chr.beagle.gz
rm $INDIR/GBR_aspatsamples_filt_05mis_uniq_scaffold_*.beagle.gz


#2. Create single beagle file for GBR+NC data
INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/ANGSD_files/"
STATDIR="/nvme/disk0/lecellier_data/WGS_NC_data/BAM_statistics/"

cat "/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/nc_aspat_noclones_amil.filelist.txt" > "/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/gbr_nc_aspat_noclones_amil.filelist.txt"
cat "/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/gbr_aspat_noclones_amil.filelist.txt" >> "/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/gbr_nc_aspat_noclones_amil.filelist.txt"
readarray -t FILES <  "/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/gbr_nc_aspat_noclones_amil.filelist.txt"

#Get number of files
N_FILES="${#FILES[@]}"
MIN_N=$((95*N_FILES/100))
MAX_DEPTH=50

CHROMOSOME_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/chromosomes_header.txt"
CONTIG=`sed -n ${SLURM_ARRAY_TASK_ID}p $CHROMOSOME_FILE`

TODO="-doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1"

FILTERS="-uniqueOnly 0 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd $MIN_N -setMinDepthInd 3 -SNP_pval 1e-6 -setMaxDepthInd $MAX_DEPTH"

#Run ANGSD separately by chromosome
start=`date +%s`
angsd -bam $INDIR/gbr_nc_aspat_noclones_amil.filelist.txt -out "${OUTDIR}GBR_NC_aspat_amil_samples_05mis_${CONTIG}" -ref $REF_3 -r $CONTIG: -nThreads 4 $TODO $FILTERS
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Merge at the end 
if [[ "$SLURM_ARRAY_TASK_ID" == 14 ]] 
then
  cat <(zcat $INDIR/GBR_NC_aspat_amil_samples_05mis_scaffold_1.beagle.gz | head -n 1) <(zcat $INDIR/GBR_NC_aspat_amil_samples_05mis_scaffold_*.beagle.gz | grep -v -w marker) | gzip > $INDIR/GBR_NC_aspat_amil_samples_05mis_all_chr.beagle.gz
  
else
  echo Not last file 
fi


