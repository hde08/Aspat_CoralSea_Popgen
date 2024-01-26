## Scripts > Scripts for calling SNPs and Variants from Illumina WGS data .

### Modified from Iva Popovic custom scripts : https://github.com/ivapops/fastq_toBam

1 Fastqc_multiqc.sh > Run initial quality check on raw sequencing data.  
2 Trimmomatic.sh > Quality trimming and adapter removal.  
3 Bwa_mem.sh > Alignment to reference genome.  
4 picard.sh > Add read groups, mark and remove duplicates.  
5 Extract_bam_unmapped.sh > Extract unmapped reads in different bam file.  
6 BAM_statistics.sh > Generate statistics on raw sequence, quality, mapping, coverage etc. 

## Analyses scripts > Post fastq processing scripts .

1 BAM_statistics.rmd > Create diagnostic plots of mapped BAM files statistics.  



