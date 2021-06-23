# sigma3_mutants_RNAseq_analysis
This repo has the scripts used to analyze the RNAseq data generated from Guo Y et al., PLoS Pathogens, 2021

The raw sequencing data can be found on NCBIs SRA under BioProject PRJNA699030

The following scripts perform:

-Initial quality filtering and adapter trimming (this is adapted from the Stenglein lab taxonomic assessment pipeline: https://github.com/stenglein-lab/taxonomy_pipeline)

-Alignment using HISAT2 and creation of read count table using FeatureCounts

-RNAseq analysis using DESeq2

