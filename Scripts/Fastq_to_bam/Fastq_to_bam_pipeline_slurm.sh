#!/bin/bash
#Execute file : sh /home/hdenis/Coral-Genomics/Scripts/Fastq_to_bam/Fastq_to_bam_pipeline_slurm.sh

#1. Quality check with fastqc/multiqc

#2. Quality trimming and adapter removal using Trimmomatic v0.39

#3. Post-trimming quality check

##4. Map to reference genome using Bwa mem (version 0.7.17-r1188) with default parameters
bwajob=$(sbatch /home/hdenis/Coral-Genomics/Scripts/Fastq_to_bam/Bwa_mem_slurm_v2.sh)
bwajob_id=${bwajob##* }

#5. Mark and delete PCR duplicates using picard version 2.27.5 
picardjob=$(sbatch --dependency=afterany:${bwajob_id} /home/hdenis/Coral-Genomics/Scripts/Fastq_to_bam/picard_slurm.sh)
picardjob_id=${picardjob##* }

##6. Extract unmapped reads (non host reads) from bam files to separate files 
extrjob=$(sbatch --dependency=afterany:${picardjob_id} /home/hdenis/Coral-Genomics/Scripts/Fastq_to_bam/Extract_bam_unmapped_slurm.sh)
extrjob_id=${extrjob##* }


##7. Generate BAM file statistics using samtools version 1.10 
bamstatjob=$(sbatch --dependency=afterany:${extrjob_id} /home/hdenis/Coral-Genomics/Scripts/Fastq_to_bam/BAM_statistics_slurm.sh)
bamstatjob_id=${bamstatjob##* }
