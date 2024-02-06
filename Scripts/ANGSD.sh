#!/bin/bash
#chmod u+x /home/hdenis/Coral-Genomics/Scripts/ANGSD.sh

#### 8. Perform SNP and genotype calling in ANGSD v0.1.17

cd /nvme/disk0/lecellier_data/WGS_GBR_data/
INDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/Aligned_files/"
OUTDIR="/nvme/disk0/lecellier_data/WGS_GBR_data/ANGSD_files/"

#Create filelist with full path of each BAM file (one filepath per line) 
FILES=($INDIR*Amillepora_MARKED_DUP.bam)
printf "%s\n" "${FILES[@]}" > $INDIR/bam.filelist.txt
#Get number of files 
N_FILES="${#FILES[@]}"
MIN_N=$((95*N_FILES/100))

#Reference genomes -> Pick only one of the two
REF_1="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Amillepora_ncbi_dataset/data/GCA_013753865.1/GCA_013753865.1_Amil_v2.1_genomic.fna"
REF_2="/nvme/disk0/lecellier_data/WGS_GBR_data/Ref_genomes/Aspathulata_ncbi_dataset/data/GCA_031770025.1/GCA_031770025.1_AGI_CSIRO_Aspa_v1_genomic.fna"

#8.1 Re-index reference using samtools 
samtools faidx $REF_1

#8.2 Filtering of polymorphic sites 

#Parameters : 
#--minQ 30 base quality >30
#--minMapQ mapping quality >30
#--minInd MIN_N sites with data from at least 95% of individuals
#--setMinDepth 3 minimum of 3 reads 
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
#---doMajorMinor (infer from genotype likelihoods)
#--SNP_pval 1e-6 (variant likelihood ratio <0.0000001)

#All in one single command line 

#Edit the outdir path 
start=`date +%s`
angsd -bam $INDIR/bam.filelist.txt -out $OUTDIR -ref $REF_1 -uniqueOnly 0 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 30 -minInd MIN_N -setMinDepth 3  -doCounts 1 -GL 1 -doGlf 2 -doSNPstat -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -nThreads 20
end=`date +%s`
echo Execution time was `expr $(( ($end - $start) / 60))` minutes.

