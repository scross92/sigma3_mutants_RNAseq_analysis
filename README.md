# sigma3_mutants_RNAseq_analysis
This repo has the scripts used to analyze the RNAseq data generated from Guo Y et al., PLoS Pathogens, 2021

The raw sequencing data can be found on NCBIs SRA under BioProject PRJNA699030

The following scripts perform:

-Creation of read count table using FeatureCounts after aligning with HISAT2 to the Mus musculus GRCm38 build

-RNAseq analysis using DESeq2

By running these scripts, the user should be able to reproduce the aspects of Figure 9 and Supplemental Table 3

