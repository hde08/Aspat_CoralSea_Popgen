#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Trimmomatic_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J Trimmomatic

### WALLTIME
#SBATCH -t 56:00:00

### MPI TASKS (cores)
#SBATCH -n 20

#Output and error directory
#SBATCH -o /home/hdenis/Slurm/outFile_%j.out

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

# This job's working directory :
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

ulimit -s unlimited

#### 2. Quality trimming and adapter removal using Trimmomatic v0.39
#mkdir /data1/WGS_Aspat_GBR/Trimmed_files
mkdir /nvme/disk0/lecellier_data/WGS_GBR_data/Postqfilt_quality_check/

cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/230404-A00151A_L001/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Trimmed_files/"

#List R1 files only 
FILES=($INDIR/*_1.fq.gz)

#Test scripts on a subset of 10 files 
#FILES=("${FILES[@]:160:20}")

#Save file IDs in file 
>/nvme/disk0/lecellier_data/WGS_GBR_data/Raw_data_processing/ids.txt
for FILE in ${FILES[@]}; do
	BASE=$(basename $FILE)
	BASE=${BASE%%_1*} 
	echo ${BASE} >> /nvme/disk0/lecellier_data/WGS_GBR_data/Raw_data_processing/ids.txt
done 


#Run Trimmomatic in parallel mode : assign each file to 1 CPU -> seems to work better than providing files one by one with -threads N argument
#Trim log removed as it takes too much storage -trimlog "${OUTDIR}{}_trim.log"
start=`date +%s`
cat /nvme/disk0/lecellier_data/WGS_GBR_data/Raw_data_processing/ids.txt | parallel --jobs $NPROCS "java -jar /home/hdenis/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 -summary "${OUTDIR}{}_sum.txt" "${INDIR}{}_1.fq.gz" "${INDIR}{}_2.fq.gz" "${OUTDIR}{}_R1_paired.fastq.gz" "${OUTDIR}{}_R1_unpaired.fastq.gz" "${OUTDIR}{}_R2_paired.fastq.gz" "${OUTDIR}{}_R2_unpaired.fastq.gz" ILLUMINACLIP:/home/hugo/PhD/Genomics/Raw_data_processing/Illumina_adapters_Iva_version.fa:2:30:10:4:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50"
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.


# Keep bases with phred-score quality > 20 in sliding window of 4 bp (average)
# Remove adapter sequences in user specified file
# Remove reads with length < 50bp following trimming'''

### Identify empty files after trimmomatic and delete them to avoid crashing multiqc
#Store their ID in a file to record files that have been eliminated
TRIMMED_FILES=(/nvme/disk0/lecellier_data/WGS_GBR_data/Trimmed_files/*_paired*.gz)
for FILE in ${TRIMMED_FILES[@]}; do
    NREADS=$(awk '{s++}END{print s/4}' $FILE)
    if [ "$NREADS" -eq 0 ]; then
        echo ${FILE} >> /nvme/disk0/lecellier_data/WGS_GBR_data/Raw_data_processing/trimming_empty_ids.txt
        rm $FILE
    fi   
done

#### 3. Check quality and adapter trimming
#Change threads number 
fastqc --noextract --outdir "Postqfilt_quality_check/" $(ls /nvme/disk0/lecellier_data/WGS_GBR_data/Trimmed_files/*_paired*.gz) --threads 5

#Generate one general html with multiqc
''' Warning : multiqc will crash if some fastqc reports are empty (0 sequences)'''

multiqc Postqfilt_quality_check/ -o Postqfilt_quality_check -f 
