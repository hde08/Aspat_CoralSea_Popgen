#!/bin/bash
#Execute file : sbatch /home/hdenis/Coral-Genomics/Scripts/Fastq_to_bam_pipeline_slurm.sh

#1. Quality check with fastqc/multiqc

#2. Quality trimming and adapter removal using Trimmomatic v0.39

#3. Post-trimming quality check

#4. Map to reference genome using Bwa mem (version 0.7.17-r1188) with default parameters
sbatch /home/hdenis/Coral-Genomics/Scripts/Bwa_mem_slurm_v2.sh

#5. Mark and delete PCR duplicates using picard version 2.27.5 
sbatch /home/hdenis/Coral-Genomics/Scripts/picard_slurm.sh

#6. Extract unmapped reads (non host reads) from bam files to separate files 
sbatch /home/hdenis/Coral-Genomics/Scripts/Extract_bam_unmapped_slurm.sh

#7. Generate BAM file statistics using samtools version 1.10 
sbatch /home/hdenis/Coral-Genomics/Scripts/BAM_statistics_slurm.sh