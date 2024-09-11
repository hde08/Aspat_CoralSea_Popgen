
# Fastq_to_bam pipeline

### This pipeline uses codes modified from Iva Popovic custom scripts   [https://github.com/ivapops/fastq_toBam] to produce bam files with duplicates removed from raw WGS fastq files. 


### This notebook contains a summary of the programs versions, commands used at each step, along with parameters used. Detailed scripts are available in github repo https://github.com/hde08/Coral-Genomics/tree/master/Scripts/Fastq_to_bam
> Fastq_to_bam_pipeline_slurm.sh.

### At each step a summary of resulting statistics is provided for GBR and NC samples. Detailed statistics per file are also available in a separate csv file. 
> Pre_variant_per_sample_summary.csv 
---------------------------
## 1 Initial statistics before any filtration
**Fastqc v0.12.1 and MultiQC v1.21**
> Fastqc_multiqc.sh

    fastqc --noextract --outdir "Initial_quality_check/" "${INDIR}${BASE}_R1_paired.fastq.gz" --threads 1
    fastqc --noextract --outdir "Initial_quality_check/" "${INDIR}${BASE}_R2_paired.fastq.gz" --threads 1
    multiqc "${INDIR}" -o "${OUTDIR}" -f -d



|    | GBR | NC |
|----|----|----|
| N initial samples collected| 831 |303
| N initial samples sequenced| 830 |194
| N natural bleaching samples| 0 |20
| N additional technical duplicates| 2 |0
| **Total number of aspat samples sequenced**| **832** |**314**
| N other species samples as outgroups| 80 |7
| N controls from other sequencing facility | 0 |23
| **Total number of fastq files**| **912** |**344**
| N samples sucessfully sequenced | 912 |344
| Total number of reads prior filtering| 31.7*10^9 |14.2*10^9

## 2 Quality trimming and adapter removal
**Trimmomatic v0.39**
> Trimmomatic_slurm.sh

*Parameters:*    
*#adapters.fa:2:30:10:4:true -> clip adapters using 'palindrome mode'*     
*#LEADING:3 -> remove leading low quality bases below quality 3*   
*#TRAILING:3 -> remove trailing low quality bases below quality 3*    
*#SLIDINGWINDOW:4:20 -> 4 bp sliding window and min phred score quality 20*    
*#MINLEN:50 -> trim adapters*     

    java -jar /home/hdenis/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 -summary "${OUTDIR}${BASE}_sum.txt" "${INDIR}${BASE}_1.fastq.gz" "${INDIR}${BASE}_2.fastq.gz" "${OUTDIR}${BASE}_R1_paired.fastq.gz" "${OUTDIR}${BASE}_R1_unpaired.fastq.gz" "${OUTDIR}${BASE}_R2_paired.fastq.gz" "${OUTDIR}${BASE}_R2_unpaired.fastq.gz" ILLUMINACLIP:/nvme/disk0/lecellier_data/WGS_NC_data/Illumina_adapters_Iva_version.fa:2:30:10:4:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    fastqc --noextract --outdir "Postqfilt_quality_check/" "${OUTDIR}${BASE}_R1_paired.fastq.gz" --threads 1
    fastqc --noextract --outdir "Postqfilt_quality_check/" "${OUTDIR}${BASE}_R2_paired.fastq.gz" --threads 1 



|    | GBR | NC |
|----|----|----|
| Number of samples used as input | 912 |344
| Number of reads after reads quality filtering| 29.8*10^9 |14.20*10^9
| **Number of samples kept after >10M reads filtering**| **908** |**344**

## 3 Alignment to reference genome
**Bwa mem v0.7.17**
> Bwa_mem_slurm_v2.sh

*Parameters:*  
*Default parameters*  

    bwa mem -t 5 ${REF_3} -o "${OUTDIR}${BASE}_pe_aln_${REF_NAME}.sam" "${INDIR}${BASE}_R1_paired.fastq.gz" "${INDIR}${BASE}_R2_paired.fastq.gz" 
    samtools view -b --threads 5 "${OUTDIR}${BASE}_pe_aln_${REF_NAME}.sam" | samtools sort --threads 5  -O BAM -o "${OUTDIR}${BASE}_pe_aln_${REF_NAME}_UNDEDUP.bam"

  
## 4 Add read groups, mark and remove duplicates
**Samtools v1.20 and Picard v3.1.1**
   > picard_slurm.sh
   > Extract_bam_unmapped_slurm.sh
 
 *Parameters:*  
 *#RGID : unique read group identifier*  
  *#RGLB : library index*  
  *#RGPL : platform name*  
  *#RGPU : flowcell, lane, sample barcode*   
  *#RGSM : sample name (=RGID in this case)*  
  *#VALIDATION_STRINGENCY=LENIENT*  
 *#REMOVE_DUPLICATES=true*   

    #Add Read Group
    /home/hdenis/Programs/gatk-4.5.0.0/gatk AddOrReplaceReadGroups INPUT="${INDIR}${BASE}_UNDEDUP.bam" OUTPUT="${INDIR}${BASE}_UNDEDUP_RG.bam" RGID=$(echo ${BASE_REP} | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/") RGLB=$(echo ${BASE_REP} | cut -d"-" -f1,2,3) RGPL=ILLUMINA RGPU=$(echo ${BASE_REP} | cut -d"_" -f2) RGSM=$(echo ${BASE_REP} | cut -d"_" -f1 | sed "s/RRAP-.*-202.*-A/A/")
    #Mark and remove duplicate
     /home/hdenis/Programs/gatk-4.5.0.0/gatk MarkDuplicates INPUT="${INDIR}${BASE}_UNDEDUP_RG.bam" OUTPUT="${INDIR}${BASE}_MARKED_DUP.bam" METRICS_FILE="${OUTDIR}${BASE}_marked_dup_metrics.txt" REMOVE_DUPLICATES=true TMP_DIR="${INDIR}" VALIDATION_STRINGENCY=LENIENT
    
   

|    | GBR | NC |
|----|----|----|
| Number of samples used as input | 908 |344
| Number of reads after filtering for PCR and optical duplicates| 25.3*10^9 |10.5*10^9
| Number of reads that mapped| 24.3*10^9 |97.3*10^8
| Number of reads that primarly mapped| 23.1*10^9 |92.4*10^8
| Number of reads after filtering improperly paired reads|20.5*10^9  |80.7*10^8
| **Number of samples kept after >80% map and >10M mapped reads filtering**| **903** |**316**

## 5 Separate non coral reads
**Samtools v1.20 **
> BAM_statistics_slurm.sh

 *Parameters:*
 *#-f 12 extract only the reads when both paired reads are unmapped*    
 *#-F 256 discard secondary alignments*    
  

    # generate sample names
    echo ${BASE} >> ${OUTDIR}/sample_name.txt		
    # generate read counts
    samtools view -c --threads 2 "${INDIR}${BASE}_MARKED_DUP.bam" >> ${OUTDIR}allCounts.txt		
    # generate bam with unmapped reads = symbiont, microbes and other reads 
    samtools view -b -f 12 -F 256 --threads 2 "${INDIR}${BASE}_MARKED_DUP.bam" > "${INDIR}${BASE}.unmapped.bam"
    # generate unmapped read counts
    samtools view -c --threads 2 "${INDIR}${BASE}.unmapped.bam" >> "${OUTDIR}unmappedCounts.txt"

## 6 Generate BAM statistics
**Samtools v1.20 **
> Extract_bam_unmapped_slurm.sh



    #Generate index bai and statistics from BAM aligned reads 
    samtools index -@ 2 "${INDIR}${BASE}_MARKED_DUP.bam"
    #Output statistics 
    samtools idxstats --threads 2 "${INDIR}${BASE}_MARKED_DUP.bam" > "${OUTDIR}${BASE}-index_stats.txt"
    samtools coverage "${INDIR}${BASE}_MARKED_DUP.bam" > "${OUTDIR}${BASE}-coverage.txt"
    samtools flagstat --threads 2 -O tsv "${INDIR}${BASE}_MARKED_DUP.bam" > "${OUTDIR}${BASE}-flagstats.txt"	
    #Extract mapped reads lenght distribution
    samtools view -F12 "${INDIR}${BASE}_MARKED_DUP.bam" | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c > "${OUTDIR}${BASE}-readlenghts.txt" 

Average statistics for BAM files that passed the prior filters :


|    | GBR | NC |
|----|----|----|
| Number of samples used to compute statistics | 903 | 316
| Base quality all| 35.6 |39.4
| Base quality chr| 36.0 |39.5
| Coverage all (%)| 74.2 |75.0
| Coverage chr (%)| 87.1 |88.0
| Depth all (X)| 10.9 |12.3
| Depth chr (X)| 7.1 |7.9
| Mapping quality all| 27.6 |27.8
| Mapping quality chr| 41.0 |42.0
