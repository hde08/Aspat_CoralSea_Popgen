# SNPs and Variants calling from Illumina WGS coral data.

## Scripts

### 1. Fastq_to_bam
### Modified from Iva Popovic custom scripts : https://github.com/ivapops/fastq_toBam
### Scripts format for execution using Slurm workload manager

1 Fastqc_multiqc.sh > Run initial quality check on raw sequencing data.  
2 Trimmomatic_slurm.sh & Trimmomatic_2nd_round_slurm.sh > Quality trimming and adapter removal.  
3 Bwa_mem_slurm_v2.sh.sh > Alignment to reference genome.  
4 picard_slurm.sh > Add read groups, mark and remove duplicates.  
5 Extract_bam_unmapped_slurm.sh > Extract unmapped reads in different bam file.  
6 BAM_statistics_slurm.sh > Generate statistics on raw sequence, quality, mapping, coverage etc. 
7 ANGSD.sh > Generate beagle file for population genetics analyses using genotype likelihoods 

### 1. Bam_to_vcf
### Personal scripts for SNPs and variant hard calling
### Following GATK 'Best practices' 
### https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels

1 GATK_HaplotypeCaller_slurm.sh > per sample variant calling (using complete reference chromosomes + scaffolds)  
2 GATK_GenomicsDBimport_slurm.sh > Create genomic database (chromosomes only, split by 4 or 2 intervals depending on size)  
3 GATK_GenotypeGVCFs_slurm.sh > joint genotyping across all samples (per intervals)  
4 VCF_filtration_slurm.sh > Combine multiple vcfs in single file (all chromosomes), select variable sites and produce filtered vcfs  

### 3. Other Analyses scripts.
### Various scripts for population demographics, connectivity, diversity inference and GEAs/GWAS





