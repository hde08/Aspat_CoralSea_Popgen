# Analysis of *Acropora spathulata* population genomics across the Western Pacific 

 ## DOI : 

**The authors would appreciate being notified if you intend to use these data or analyses in your own work**

### Title: 

Authors
Hugo Denis<sup>1,2,3</sup>, Katharine E. Prata<sup>4,5</sup>, Hisatake Ishida<sup>6</sup>, Iva Popovic<sup>4,5</sup>, Véronique J.L. Mocellin<sup>4</sup>, Magali Boussion<sup>1</sup>, Ilha Byrne<sup>5</sup>, Steven W. Purcell<sup>3</sup>, Line K. Bay<sup>4</sup>, Gaël Lecellier<sup>7</sup>, Cheong Xin Chan<sup>6</sup>, Cynthia Riginos<sup>4,5</sup>, Emily J. Howells<sup>3,8†</sup>, Véronique Berteaux-Lecellier<sup>1,†</sup>

1 Santé des Ecosystèmes tropicaux (SantEco) EMR9001, CNRS, IRD, Univ. de Nouvelle-Calédonie, Univ. la Reunion, UMR 250 ENTROPIE, Noumea 98800, New Caledonia, France.  
2 ED 129, SU Sorbonne Université, 4, Place Jussieu, 75252 Paris, France  
3 National Marine Science Centre, Faculty of Science and Engineering, Southern Cross University, Coffs Harbour, NSW, Australia  
4 Australian Institute of Marine Science, Townsville, QLD, Australia  
5 School of the Environment, University of Queensland, Brisbane, QLD 4072, Australia  
6 Australian Centre for Ecogenomics, School of Chemistry and Molecular Biosciences, The University of Queensland, Brisbane, QLD 4072, Australia  
7 Institut de Sciences Exactes et Appliquées (ISEA) EA7484, 145, Avenue James Cook – BP R4 98 851 Nouméa, New Caledonia  
8 Lizard Island Research Station, Australian Museum, Sydney, NSW 2010, Australia  
†: These authors contributed equally to this work  


Correspondence : Hugo DENIS, Institut de Recherche pour le Développement, Promenade Roger-Laroque, Noumea cedex, New Caledonia, France, Email : denis.hugo08@gmail.com ORCID : https://orcid.org/0000-0003-2909-6299


## WGS Data Processing & Analyses Scripts

### 1. Fastq_to_bam
### Modified from Iva Popovic custom scripts : https://github.com/ivapops/fastq_toBam
#### Scripts format for execution using Slurm workload manager

1 Fastqc_multiqc.sh > Run initial quality check on raw sequencing data.  
2 Trimmomatic_slurm.sh & Trimmomatic_2nd_round_slurm.sh > Quality trimming and adapter removal.  
3 Bwa_mem_slurm_v2.sh.sh > Alignment to reference genome.  
4 picard_slurm.sh > Add read groups, mark and remove duplicates.  
5 Extract_bam_unmapped_slurm.sh > Extract unmapped reads in different bam file.  
6 BAM_statistics_slurm.sh > Generate statistics on raw sequence, quality, mapping, coverage etc. 

### 2. Bam_to_vcf
### Script for variant 'hard calling'
#### Following GATK 'Best practices' 
#### https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels

1 GATK_HaplotypeCaller_slurm.sh > per sample variant calling (using complete reference chromosomes + scaffolds)    
2 GATK_GenomicsDBimport_slurm.sh > Create genomic database (chromosomes only, split by 4 or 2 intervals depending on size)    
3 GATK_GenotypeGVCFs_slurm.sh > joint genotyping across all samples (per intervals)   
4 GATK_GatherVCFs_SelectVariants_slurm.sh > Combine intervals vcf files and select variants     
5 VariantFiltration_slurm.sh > Filter final vcf by variants, genotypes and individuals   

### 3. Popgen_angsd
### Scripts to have a preliminary look at the data using a genotype likelihood framework
1 ANGSD_slurm.sh > Generate beagle file using ANGSD software  
2 PCangsd_slurm.sh > PCA and ADMIXTURE analyses using pcangsd  
3 IBS_slurm.sh > Compute Identity by state matrix in ANGSD  

### 4. Popgen_vcf
### Scripts to conduct final population genetics analyses presented in the manuscript using 'hard called' variants. 
1 Compare_sequencing_batch.sh > Confirm absence of strong batch effects 
2 ADMIXTURE_slurm.sh > Perform ADMIXTURE analyses on filtered vcf file   
3 Stats_pixy_slurm.sh > Compute Fst, Dxy, Pi in PIXY   
4 RAxML_slurm.sh > Compute ML-tree in RAxML    
5 Original scripts to perform dadi analyses can be found at https://github.com/kepra3/kp_dadi/tree/ahya-demo

### 5. Symbionts  
### Scripts to analyze symbiont communities composition and diversity from hologenome reads. 
#### Based on original scripts from Hisatake Ishida available at : https://github.com/hisatakeishida/Symb-SHIN
1 Preprocess_reads.sh > Separate Symbiodiniaceae reads from host reads  
2 Kmer_symbionts.sh > Generate k-mer counts per sample and compute D2S-derived distance    
3 GraftM_ITS2_slurm.sh > Recover ITS2 reads in hologenome data using sequence classifier tool GraftM    
4 psbAncr.sh > Map symbiont reads to psbAncr reference database and use uniquely mapped reads as relative abundance       
5 Symbionts_symportal_hologenome_comparison.Rmd > Compare Symbiodiniaceae diversity inferred from hologenome-based methods and ITS2 amplicon sequencing (Figure S8)  
6 Symbiont_NJ_tree_multi_markers_barplots.Rmd > Make NJ-tree based on D2S-derived distance and marker genes barplots (Figure S9a)    

### 6. Stats_Plots
### Scripts to conduct statistical analyses and create main plots from the publication
1 Host_structure_plot.Rmd > Make Figure 1  
2 IbD_Ne_plots.Rmd > Make Figure 2  
3 Host_symbiont_tanglegram.Rmd > Make Figure 3  
4 Symbiont_structure_plot.Rmd > Make Figure 4  
5 RDA_analyses.Rmd > Scripts to conduct RDA analyses and produce results from Table 2

