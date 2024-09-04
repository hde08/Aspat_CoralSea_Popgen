#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Fastq_to_bam/ANGSD_NC_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J ANGSD

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=100G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/angsd_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/angsd_%A_%a.e     # standard error

#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=4		# number of cores per job

#SBATCH --array=1-14%1        	# job array

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
#### On NC samples ####

#INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/"
#OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/ANGSD_files/"
#STATDIR="/nvme/disk0/lecellier_data/WGS_NC_data/BAM_statistics/"
#
#mkdir "${OUTDIR}"
#mkdir "${STATDIR}"
#
##Create filelist with full path of each BAM file (one filepath per line)
#FILES=($INDIR*MARKED_DUP.bam)
#printf "%s\n" "${FILES[@]}" > $INDIR/bam.filelist.txt
#
##Remove one file with low number of reads (<1 M : Aspat-NDIR-489)
##cat $INDIR/bam.filelist.txt | grep 'Aspat-NDIR-489' -v > $INDIR/bam.filelist2.txt
#
#
#readarray -t FILES <  $INDIR/bam.filelist.txt
##readarray -t FILES <  $INDIR/aspat_bam_clean.filelist.txt
#
##Get number of files
#N_FILES="${#FILES[@]}"
#MIN_N=$((95*N_FILES/100))
#
#DEPTH_FILES="${FILES[@]%%_MARKED*}"
#
##Get median and sd depth value across bam files 
##Store genome-wide mean depth for each fle 
##DEPTH_FILES=($STATDIR*Amillepora*coverage.txt)
#> $STATDIR/bamfile_depth.txt
#> $STATDIR/bamfile_sd_depth.txt
#
#for FILE in ${DEPTH_FILES[@]}; do
#    FILE="$(basename $FILE)"
#    FILE=$STATDIR$FILE
#    FILE+="-coverage.txt"
#    MEAN_DEPTH=$(awk '{sum+=$7} END { print sum/NR}' $FILE)
#    SD_DEPTH=$(awk '{delta = $7 - avg; avg += delta / NR; mean2 += delta * ($7 - avg); } END { print sqrt(mean2 / (NR-1)); }' $FILE)
#    echo $MEAN_DEPTH >> $STATDIR/bamfile_depth.txt
#    echo $SD_DEPTH>> $STATDIR/bamfile_sd_depth.txt
#done
#
##Set max depth as the median genome-wide depth + 1 SD (following Lou 2021 guidelines)
#MEDIAN_DEPTH=$(sort -k1 -n $STATDIR/bamfile_depth.txt | awk '{arr[NR]=$1}
#   END { if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}')
#MEDIAN_SD=$(sort -k1 -n $STATDIR/bamfile_sd_depth.txt | awk '{arr[NR]=$1}
#   END { if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}')
#MAX_DEPTH=`echo "$MEDIAN_DEPTH + $MEDIAN_SD" | bc -l`

#A.millepora v3 reference genome 
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"

REF_NAME="Amilleporav3"

#8.1 Re-index reference using samtools
#samtools faidx $REF_3

#8.2 Filtering of polymorphic sites

#Parameters :
#--minQ 30 base quality >30
#--minMapQ mapping quality >30
#--minInd MIN_N sites with data from at least 95% of individuals
#--setMinDepthInd 3 minimum of 3 reads for each individual 
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
#---doMajorMinor 1 (inferred directly from genotype likelihoods)
#--SNP_pval 1e-6 (variant likelihood ratio <0.0000001)


#Keep only SNPs that matched on the first 14 scaffolds (chromosomes)

#cat $REF_3 | grep -m 14  '>' | cut -c 2- > $CHROMOSOME_FILE

#Associate slurm array index with contig (chr) name 
#CHROMOSOME_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/chromosomes_header.txt"
#CONTIG=`sed -n ${SLURM_ARRAY_TASK_ID}p $CHROMOSOME_FILE`

#Run a first time on all samples from this sequencing run 

#TODO="-doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1"
#
#FILTERS="-uniqueOnly 0 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd $MIN_N -setMinDepthInd 3 -SNP_pval 1e-6 -setMaxDepthInd $MAX_DEPTH"

#start=`date +%s`
#angsd -bam $INDIR/bam.filelist.txt -out "${OUTDIR}NC_allsamples_filt_05mis_uniq_${CONTIG}" -ref $REF_3 -r $CONTIG: -nThreads 5 $TODO $FILTERS
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Run a second time on all samples and output only sites with MAF >0.05 

#start=`date +%s`
#angsd -bam $INDIR/bam.filelist.txt -out "${OUTDIR}NC_allsamples_filt_05mis_MAF0.05_${CONTIG}" -ref $REF_3 -r $CONTIG: -nThreads 5 $TODO $FILTERS -minMaf 0.05
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Previous version with low filtering 
#start=`date +%s`
#angsd -bam $INDIR/bam.filelist.txt -out "${OUTDIR}GBR_allsamples_${CONTIG}" -ref $REF_3 -r $CONTIG: -uniqueOnly 0 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd MIN_N -setMinDepthInd 3  -doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -nThreads 5 -setMaxDepthInd $MAX_DEPTH
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#New version with more stringent filtering
#No reads with multiple hit 
#No missing data 
#Minimum 5 reads per variant

##Run a first time on the different species, allowing for some missing data
#start=`date +%s`
#angsd -bam $INDIR/bam.filelist2.txt -out "${OUTDIR}GBR_allsamples_filt_uniq_${CONTIG}" -ref $REF_3 -r $CONTIG: -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd $MIN_N -setMinDepthInd 3  -doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -nThreads 5 -setMaxDepthInd $MAX_DEPTH
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Note on beagle file coding
#0=A, 1=C, 2=G, 3=T


#Run a second time on aspat samples only with and without missing data -> THink of doing IBS at the same time 

#TODO="-doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1"
#
#FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd $MIN_N -setMinDepthInd 3 -SNP_pval 1e-6 -setMaxDepthInd $MAX_DEPTH"
#
#start=`date +%s`
#angsd -bam $INDIR/aspat_bam_clean.filelist.txt -out "${OUTDIR}GBR_aspatsamples_filt_05mis_uniq_${CONTIG}" -ref $REF_3 -r $CONTIG: -nThreads 5 $TODO $FILTERS
#end=`date +%s`
#echo Execution time was `expr $(( ($end - $start) / 60))` minutes.


#8.3 Run ANGSD analyses on NC + GBR samples jointly
#cat "/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/nc_aspat_noclones_amil.filelist.txt" > "/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/gbr_nc_aspat_noclones_amil.filelist.txt"
#cat "/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/gbr_aspat_noclones_amil.filelist.txt" >> "/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/gbr_nc_aspat_noclones_amil.filelist.txt"

readarray -t FILES <  "/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/gbr_nc_aspat_noclones_amil.filelist.txt"


INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/ANGSD_files/"
STATDIR="/nvme/disk0/lecellier_data/WGS_NC_data/BAM_statistics/"

#Get number of files
N_FILES="${#FILES[@]}"
MIN_N=$((95*N_FILES/100))

MAX_DEPTH=50


CHROMOSOME_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/chromosomes_header.txt"
CONTIG=`sed -n ${SLURM_ARRAY_TASK_ID}p $CHROMOSOME_FILE`

#Run a first time on all samples from this sequencing run 

TODO="-doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1"

FILTERS="-uniqueOnly 0 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd $MIN_N -setMinDepthInd 3 -SNP_pval 1e-6 -setMaxDepthInd $MAX_DEPTH"

INDIR="/nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_NC_data/ANGSD_files/"

start=`date +%s`
angsd -bam $INDIR/gbr_nc_aspat_noclones_amil.filelist.txt -out "${OUTDIR}GBR_NC_aspat_amil_samples_05mis_${CONTIG}" -ref $REF_3 -r $CONTIG: -nThreads 4 $TODO $FILTERS
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.
if else 

if [[ "$SLURM_ARRAY_TASK_ID" == 14 ]] 
then
  cat <(zcat /nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/GBR_NC_aspat_amil_samples_05mis_scaffold_1.beagle.gz | head -n 1) <(zcat /nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/GBR_NC_aspat_amil_samples_05mis_scaffold_*.beagle.gz | grep -v -w marker) | gzip > /nvme/disk0/lecellier_data/WGS_NC_data/Aligned_files/GBR_NC_aspat_amil_samples_05mis_all_chr.beagle.gz
  
else
  echo Not last file 
fi
