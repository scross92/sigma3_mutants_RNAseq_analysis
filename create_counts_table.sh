#!/bin/bash

#this script uses featureCounts to create a counts table from the HISAT2 output

gtf=./../references/Mus_musculus.GRCm38.95.gtf

featureCounts -a $gtf -o counts.txt -T 24 -p -t exon -g gene_id YG1_R1_f.fastq.YG1_R2_f.fastq.hs_mouse_genome.bam YG2_R1_f.fastq.YG2_R2_f.fastq.hs_mouse_genome.bam YG3_R1_f.fastq.YG3_R2_f.fastq.hs_mouse_genome.bam YG4_R1_f.fastq.YG4_R2_f.fastq.hs_mouse_genome.bam YG5_R1_f.fastq.YG5_R2_f.fastq.hs_mouse_genome.bam YG6_R1_f.fastq.YG6_R2_f.fastq.hs_mouse_genome.bam YG7_R1_f.fastq.YG7_R2_f.fastq.hs_mouse_genome.bam YG8_R1_f.fastq.YG8_R2_f.fastq.hs_mouse_genome.bam YG9_R1_f.fastq.YG9_R2_f.fastq.hs_mouse_genome.bam YG10_R1_f.fastq.YG10_R2_f.fastq.hs_mouse_genome.bam YG11_R1_f.fastq.YG11_R2_f.fastq.hs_mouse_genome.bam YG12_R1_f.fastq.YG12_R2_f.fastq.hs_mouse_genome.bam

