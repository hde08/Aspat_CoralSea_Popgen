#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Bam_to_vcf/GATK_GatherVCFs_SelectVariants_slurm.sh

##SBATCH --clusters=dell_r740xd

### partition
#SBATCH --nodelist=R740xd
#SBATCH --partition=local

### JOB NAME
#SBATCH -J GatherVCFs

### WALLTIME
#SBATCH -t 3-00:00:00

#MEMORY
#SBATCH --mem=300G # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)


#Output and error directory
#SBATCH -o /home/hdenis/Slurm/gathvcf_%A_%a.o   #Standard output 
#SBATCH -e /home/hdenis/Slurm/gathvcf_%A_%a.e     # standard error

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

#11 Merge vcf files and select variants gatk 4.5.0.0

cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/GATK_files/Vcf_files/"

#A.millepora v3 reference genome 
REF_3="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amil_scaffolds_final_v3.fa"
REF_NAME="Amilleporav3"

#Associate slurm array index with contig (chr) name 
CHROMOSOME_FILE="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/chromosomes_header.txt"
CONTIG=`sed -n ${SLURM_ARRAY_TASK_ID}p $CHROMOSOME_FILE`

##Create list of files to be merged (corresponds to each intervals file names)
#N_VCF_FILES=14
#echo -n > "${INDIR}vcf_list.txt"
#for K in $(seq 1 $N_VCF_FILES);
#do
#  if (( $K < 5 )); then
#    printf "%s\t%s\n" "${INDIR}aspat_clean_scaffold_${K}a.vcf.gz" >> "${INDIR}vcf_list.txt"
#    printf "%s\t%s\n" "${INDIR}aspat_clean_scaffold_${K}b.vcf.gz" >> "${INDIR}vcf_list.txt"
#    printf "%s\t%s\n" "${INDIR}aspat_clean_scaffold_${K}c.vcf.gz" >> "${INDIR}vcf_list.txt"
#    printf "%s\t%s\n" "${INDIR}aspat_clean_scaffold_${K}d.vcf.gz" >> "${INDIR}vcf_list.txt"
#  else
#    printf "%s\t%s\n" "${INDIR}aspat_clean_scaffold_${K}a.vcf.gz" >> "${INDIR}vcf_list.txt"
#    printf "%s\t%s\n" "${INDIR}aspat_clean_scaffold_${K}b.vcf.gz" >> "${INDIR}vcf_list.txt"
#  fi
#
#done


#11.1 Combine multiple VCF (per chromosome) into single VCF
start=`date +%s`
echo start merging vcf files 

singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx200g" GatherVcfs --INPUT "${INDIR}vcf_list.txt" --OUTPUT "${INDIR}aspat_clean_all_chr.vcf.gz" --CREATE_INDEX true 

end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#11.2 Index combined file 
start=`date +%s`
echo start indexing merged file 

singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx200g" IndexFeatureFile --input "${INDIR}aspat_clean_all_chr.vcf.gz"

end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#11.3 Select variants (as we included monomorphic sites with --include-non-variant-sites)
#Output file with SNPs only
start=`date +%s`
echo start selecting SNPs 

singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx200g" SelectVariants --variant "${INDIR}aspat_clean_all_chr.vcf.gz" --output "${INDIR}aspat_clean_all_chr_SNP.vcf.gz" --select-type-to-include SNP

echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

#Output file with INDELs only 
start=`date +%s`
echo start selecting INDELs 

singularity exec --bind /nvme/disk0/lecellier_data:/nvme/disk0/lecellier_data /home/hdenis/gatk_latest.sif gatk --java-options "-Xmx200g" SelectVariants --variant "${INDIR}aspat_clean_all_chr.vcf.gz" --output "${INDIR}aspat_clean_all_chr_INDEL.vcf.gz" --select-type-to-include INDEL

end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

