#Analyze RNA-seq data from various T1L reovirus mutants along with the WT virus

#Shaun Cross
#December 4, 2020

#save the workspace
save.image(file='RNAseq_yingying_env.RData')
#load saved environment
load('RNAseq_yingying_env.RData')

#import the featureCounts data
countdata <- read.table("counts.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Remove all characters after variable in filenames as it is additional info not needed
colnames(countdata) <- gsub("_.*", "", colnames(countdata))

library(tidyverse)

#Matrices of groups that I will compare and run the DESeq2 pipeline on
#Also, these have to be in the order of ***infected columns and then uninfected columns**** for DESeq2 analysis (see below)
head(countdata)
class(countdata)
as.matrix(countdata)
#group by experiment--3 groups per each condition (important for function below)
#just as a note, if you run the DESeq2 function below, it loads "MASS" library which masks select command from dplyr

mockvWT <- countdata %>% select(4:6, 1:3)
mockv287 <- countdata %>% select(7:9, 1:3)
mockv296 <- countdata %>% select(10:12, 1:3)
WTv287 <- countdata %>% select(7:9, 4:6)
WTv296 <- countdata %>% select(10:12, 4:6)
T1L287v296 <- countdata %>% select(7:12) ##This will show genes up or down reg in 287 mutants vs 296

#upload gene database for adding to DEGs
#this was taken from http://www.informatics.jax.org/downloads/reports/index.html#go MGI Gene Model Coordinates
#I then cut out just Ensembl IDs and gene names
gene_db <- read_tsv("ensembl_gene_reference.tsv")

#Turning DESeq2 Pipeline into a function. This way I can just run the function on the appropriate/desired comparisons
#DEPENDENCIES: In order to show log fold changes of infected as comapred to controls, matrices MUST be infected columns first, then uninfected

#This script is taken predominately from a DeSeq2 template by Stephen Turner found here:
#https://gist.github.com/stephenturner/f60c1934405c127f09a6
#additionally pulled from Mark Stenglein script found on the 104 server:
#/home/mdstengl/datasets/2020_Fabricio_RNA-seq/de_analysis_with_DeSeq_2.R

#make the matrix name a variable for downstream saving of plots
DESeq2_analysis <- function(countdata_group, condition_count) {
  var_name <- deparse(substitute(countdata_group))
  
  # Assign conditions
  #this is setting the reference level to the infected group
  #this means any log fold changes is how the infected flies appear compared to the control
  #e.g. gene x has a +2 log fold change, so that gene is expressed MORE in infected flies
  (condition <- factor(c(rep("Experiment", condition_count), rep("Control", condition_count))))
  
  # Analysis with DESeq2 ----------------------------------------------------
  library("DESeq2")
  # Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
  (coldata <- data.frame(row.names=colnames(countdata_group), condition))
  dds <- DESeqDataSetFromMatrix(countData=countdata_group, colData=coldata, design=~condition)
  dds
  
  # Run the DESeq pipeline
  dds <- DESeq(dds)
  
  #See the conditions. Should usually read "condition_Infect_vs_Uninfect"
  #return(resultsNames(dds))
  
  # Plot dispersions
  png(paste0("qc-dispersions_", var_name, ".png"), 1000, 1000, pointsize=20)
  plotDispEsts(dds, main="Dispersion plot")
  dev.off()
  
  # Regularized log transformation for clustering/heatmaps, etc
  rld <- rlogTransformation(dds)
  head(assay(rld))
  hist(assay(rld))
  
  # Colors for plots below
  library(RColorBrewer)
  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
  
  # Sample distance heatmap
  sampleDists <- as.matrix(dist(t(assay(rld))))
  library(gplots)
  pdf(paste0("qc-heatmap-samples_", var_name, ".pdf"), w=11, h=11, pointsize=20)
  heatmap.2(as.matrix(sampleDists), key=T, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[condition], RowSideColors=mycols[condition],
            density.info = "none",
            margin=c(10, 10), main="Sample Distance Matrix")
  dev.off()
  
  
  # Principal components analysis
  ## Could do with built-in DESeq2 function:
  
  p_PCA <- DESeq2::plotPCA(rld, intgroup="condition") +
    scale_color_manual(values = c("#56B4E9", "#D55E00")) +
    ylim(c(-20, 20)) +
    geom_text(aes(label=name),vjust=2)
  ggsave(paste0("qc-PCA-samples_", var_name, ".pdf"), plot = p_PCA, width = 5, height = 5, units = "in")
  
  dev.off()
  # Get differential expression results
  res <- results(dds)
  table(res$padj<0.05)
  ## Order by adjusted p-value
  res <- res[order(res$padj), ]
  
  ## Merge with normalized count data
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  head(resdata)
  
  ## Write results
  write.csv(resdata, file= paste0("diffexpr-results_", var_name, ".csv"))
  
  ### Run pheatmap
  ##look at heatmap clustering
  ### Set a color palette
  heat.colors <- brewer.pal(6, "RdBu")
  library(pheatmap)
  padj.cutoff <- 0.05
  lfc.cutoff <- .5849 # remember this is a log2 fold so this is a 1.5 fold log change
  resdata_threshold <- subset(resdata, abs(log2FoldChange) > lfc.cutoff & padj < padj.cutoff)
  row.names(resdata_threshold) <- resdata_threshold$Gene
  resdata_threshold <- resdata_threshold[,8:ncol(resdata_threshold)]
  
  png(paste0("threshold-heatmap-samples_", var_name, ".png"), w=1000, h=1000, pointsize=20)
  pheatmap(resdata_threshold, color = heat.colors, cluster_rows = T, show_rownames=F,
           border_color=NA, fontsize = 10, scale="row",
           fontsize_row = 10, height=20)
  dev.off()
  
  ## Examine plot of p-values
  hist(res$pvalue, breaks=50, col="grey")
  
  ## Examine independent filtering
  # attr(res, "filterThreshold")
  # plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")
  
  ## MA plot
  ## Could do with built-in DESeq2 function:
  #DESeq2::plotMA(dds, ylim=c(-1,1))
  ## This plot is from Stephen Turner script:
  library(calibrate)
  maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
    with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
    with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
    if (labelsig) {
      require(calibrate)
      with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
    }
  }
  png(paste0("diffexpr-maplot_", var_name, ".png"), 1500, 1000, pointsize=20)
  maplot(resdata, main="MA Plot")
  dev.off()
  
  #volcano plot using ggplot2 (these plots are easily modified)
  res_df <- as.data.frame(res)
  #not all padj values have a number...these are just highly insignificant
  #replace NAs with .99 for downstream plotting
  res_df$padj[is.na(res_df$padj)] <- 0.99
  #add a threshold column for proper coloring of dots
  #likely overly complex use of ifelse statements...
  res_df$threshold <- ifelse(res_df$padj<0.05, 
                             ifelse(abs(res_df$log2FoldChange)>1.5,"both", "pval"), 
                             ifelse(abs(res_df$log2FoldChange)>1.5,"lfc", "none"))
  # #now an issue arises if the padj is NA, then it doesn't work, so it then removes ~7K genes
  # #this next step addresses that issue
  # res_df$threshold <- ifelse(res_df$padj=="NA", 
  #                            ifelse(abs(res_df$log2FoldChange)>1,"lfc", "none"), 
  #                            res_df$threshold)
  
  
  
  p_volcano <- ggplot(res_df) +
    geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
    ggtitle(paste0("Volcano Plot ", var_name)) +
    xlab("Log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    scale_color_manual(values = c('#009E73', '#D55E00', 'black', '#E69F00')) +
    #scale_y_continuous(limits = c(0,50)) +
    theme_bw() +
    xlim(c(-5,5)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))  
  
  
  #p_volcano
  ggsave(paste0("diffexpr-volcanoplot_", var_name, ".pdf"), plot = p_volcano, width = 2.5, height = 2.5, units = "in")
  
  #Make files of upregulated and downregulated genes
  #combine with gene database so that gene names and descriptions are attached
  
  upreg_genes <- res[which(res$log2FoldChange > 1.5 & res$padj < .05),]
  downreg_genes <- res[which(res$log2FoldChange < -1.5 & res$padj < .05),]
  
  upreg_genes <-as.data.frame(upreg_genes)
  upreg_genes <- rownames_to_column(upreg_genes, "Ensembl gene id")
  
  downreg_genes <-as.data.frame(downreg_genes)
  downreg_genes <- rownames_to_column(downreg_genes, "Ensembl gene id")
  
  upreg_genes_info <- left_join(upreg_genes, gene_db)
  downreg_genes_info <- left_join(downreg_genes, gene_db)
  
  write.csv(upreg_genes_info, file= paste0("upreg_genes_", var_name, ".csv"))
  write.csv(downreg_genes_info, file= paste0("downreg_genes_", var_name, ".csv"))
  
}
dev.off()

#DESeq2 analysis for individual groups
DESeq2_analysis(mockvWT, condition_count = 3)
DESeq2_analysis(mockv287, condition_count = 3)
DESeq2_analysis(mockv296, condition_count = 3)
DESeq2_analysis(WTv287, condition_count = 3)
DESeq2_analysis(WTv296, condition_count = 3)
DESeq2_analysis(T1L287v296, condition_count = 3)


###VENN DIAGRAM OF SHARED AND UNIQUE GENES####

#the following code creates a venn diagram across the 3 different T1L viruses gene expression for upregulated and downregulated genes
#before this can be run, the DESeq2 function above needs to be run on the samples

library(VennDiagram)
#upregulated plot
upregWT <- read.csv('upreg_genes_mockvWT.csv')
upreg287 <- read.csv('upreg_genes_mockv287.csv')
upreg296 <- read.csv('upreg_genes_mockv296.csv')

pdf("venn_diagram_upregulated_genes.pdf")
venn_plot_upreg <- venn.diagram(list(upregWT$Ensembl.gene.id, upreg287$Ensembl.gene.id, upreg296$Ensembl.gene.id), 
                          NULL, fill=c("light grey", "#FFC20A", "#0C7BDC"), 
                          alpha=c(0.5,0.5, 0.5), cex = 2, 
                          cat.fontface=4, category.names=c("Wild-type", "287", "296"))
grid.draw(venn_plot_upreg)
dev.off()

#collect list of genes that are shared by all groups
shared_upreg_genes <- data.frame(intersect(intersect(upregWT$Ensembl.gene.id, upreg287$Ensembl.gene.id), upreg296$Ensembl.gene.id))
colnames(shared_upreg_genes) <- c("Ensembl gene id")
shared_upreg_genes <- left_join(shared_upreg_genes, gene_db)
write.csv(shared_upreg_genes, "shared_upreg_genes.csv")

#shared upreg genes between WT and 296
#this is done because 287 DOESN'T cause myocarditis, while 296 and WT do
shared_WT_296 <- data.frame(intersect(upregWT$Ensembl.gene.id, upreg296$Ensembl.gene.id))
colnames(shared_WT_296) <- c("Ensembl gene id")
shared_WT_296 <- subset(shared_WT_296, !(shared_WT_296$`Ensembl gene id` %in% upreg287$Ensembl.gene.id))
shared_WT_296 <- left_join(shared_WT_296, gene_db)
write.csv(shared_WT_296, "shared_upreg_genes_WT_296.csv")

#shared upreg genes between WT and 287
shared_WT_287 <- data.frame(intersect(upregWT$Ensembl.gene.id, upreg287$Ensembl.gene.id))
colnames(shared_WT_287) <- c("Ensembl gene id")
shared_WT_287 <- subset(shared_WT_287, !(shared_WT_287$`Ensembl gene id` %in% upreg296$Ensembl.gene.id))
shared_WT_287 <- left_join(shared_WT_287, gene_db)
write.csv(shared_WT_287, "shared_upreg_genes_WT_287.csv")

#shared upreg genes between 287 and 296
shared_287_296 <- data.frame(intersect(upreg296$Ensembl.gene.id, upreg287$Ensembl.gene.id))
colnames(shared_287_296) <- c("Ensembl gene id")
shared_287_296 <- subset(shared_287_296, !(shared_287_296$`Ensembl gene id` %in% upregWT$Ensembl.gene.id))
shared_287_296 <- left_join(shared_287_296, gene_db)
write.csv(shared_287_296, "shared_upreg_genes_287_296.csv")


#now collect list of genes that are unique for each group
#WT unique genes
unique_upreg_genes_WT <- subset(upregWT, !(Ensembl.gene.id %in% upreg287$Ensembl.gene.id))
unique_upreg_genes_WT <- subset(unique_upreg_genes_WT, !(Ensembl.gene.id %in% upreg296$Ensembl.gene.id))
write.csv(unique_upreg_genes_WT, "unique_upreg_genes_WT.csv")

#287 unique genes
unique_upreg_genes_287 <- subset(upreg287, !(Ensembl.gene.id %in% upregWT$Ensembl.gene.id))
unique_upreg_genes_287 <- subset(unique_upreg_genes_287, !(Ensembl.gene.id %in% upreg296$Ensembl.gene.id))
write.csv(unique_upreg_genes_287, "unique_upreg_genes_287.csv")

#296 unique genes
unique_upreg_genes_296 <- subset(upreg296, !(Ensembl.gene.id %in% upregWT$Ensembl.gene.id))
unique_upreg_genes_296 <- subset(unique_upreg_genes_296, !(Ensembl.gene.id %in% upreg287$Ensembl.gene.id))
write.csv(unique_upreg_genes_296, "unique_upreg_genes_296.csv")


#downregulated plot
downregWT <- read.csv('downreg_genes_mockvWT.csv')
downreg287 <- read.csv('downreg_genes_mockv287.csv')
downreg296 <- read.csv('downreg_genes_mockv296.csv')

pdf("venn_diagram_downregulated_genes.pdf")
venn_plot_downreg <- venn.diagram(list(downregWT$Ensembl.gene.id, downreg287$Ensembl.gene.id, downreg296$Ensembl.gene.id), 
                                NULL, fill=c("light grey", "#FFC20A", "#0C7BDC"), 
                                alpha=c(0.5,0.5, 0.5), cex = 2, 
                                cat.fontface=4, category.names=c("Wild-type", "287", "296"))
grid.draw(venn_plot_downreg)
dev.off()

#there aren't any shared downregulated genes across all three groups, so just pulling out unique genes
#WT unique genes
unique_downreg_genes_WT <- subset(downregWT, !(Ensembl.gene.id %in% downreg287$Ensembl.gene.id))
unique_downreg_genes_WT <- subset(unique_downreg_genes_WT, !(Ensembl.gene.id %in% downreg296$Ensembl.gene.id))
write.csv(unique_downreg_genes_WT, "unique_downreg_genes_WT.csv")

#287 unique genes
unique_downreg_genes_287 <- subset(downreg287, !(Ensembl.gene.id %in% downregWT$Ensembl.gene.id))
unique_downreg_genes_287 <- subset(unique_downreg_genes_287, !(Ensembl.gene.id %in% downreg296$Ensembl.gene.id))
write.csv(unique_downreg_genes_287, "unique_downreg_genes_287.csv")

#296 unique genes
unique_downreg_genes_296 <- subset(downreg296, !(Ensembl.gene.id %in% downregWT$Ensembl.gene.id))
unique_downreg_genes_296 <- subset(unique_downreg_genes_296, !(Ensembl.gene.id %in% downreg287$Ensembl.gene.id))
write.csv(unique_downreg_genes_296, "unique_downreg_genes_296.csv")


####PlOT LOG FOLD CHANGE COMPARED TO WT VIRUS#####
#Yingying is interested in seeing if logfold change differences of genes as compared to the WT virus may be responsible for myocarditis
#Here I am plotting the change in log fold difference of the mutants compared to the WT (Mutant_logfold - WT_logfold)
#plotted just as randomly assigned gene #s (1-212 for the shared genes)

#first pull out only the shared genes
shared_upreg_genes_287 <- left_join(shared_upreg_genes, upreg287, by = c("Ensembl gene id" = "Ensembl.gene.id"))
shared_upreg_genes_296 <- left_join(shared_upreg_genes, upreg296, by = c("Ensembl gene id" = "Ensembl.gene.id"))
shared_upreg_genes_WT <- left_join(shared_upreg_genes, upregWT, by = c("Ensembl gene id" = "Ensembl.gene.id"))

#then subtract the WT logfold change from the 287 and 296 and add extra column to specify difference for combining
difference_287 <- data.frame(shared_upreg_genes_287$log2FoldChange-shared_upreg_genes_WT$log2FoldChange)
colnames(difference_287) <- c("logfold_difference")
difference_287['mutant'] = "mut_287"
difference_287 <- tibble::rowid_to_column(difference_287, "gene_number")

difference_296 <- data.frame(shared_upreg_genes_296$log2FoldChange-shared_upreg_genes_WT$log2FoldChange)
colnames(difference_296) <- c("logfold_difference")
difference_296['mutant'] <- "mut_296"
difference_296 <- tibble::rowid_to_column(difference_296, "gene_number")

#then combine together in a new dataframe
combined_difference <- full_join(difference_287, difference_296)

#then plot
p_logfold_diff <- ggplot(combined_difference, aes(x=gene_number, y=logfold_difference, group = mutant)) +
  geom_point(aes(shape = mutant, color = mutant)) +
  scale_shape_manual(values=c(18, 17)) +
  theme_classic() +
  labs(title = "Logfold Differences of Mutants compared to WT",
       y = "Logfold difference (Mutant - WT)") +
  scale_color_manual(values = c("#FFC20A", "#0C7BDC")) +
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_logfold_diff

ggsave("logfold_difference_mutant_plot.pdf", height = 6, width = 10, units = "in")


#this logfold plot only gives a general idea of what sorts of changes are occurring.
#However, these are not statistically tested. Here I rerun the DESeq2 pipeline but:
#1) focus solely on WTvsMutants
#2) plot volcano plot of ONLY the 212 genes of interest

shared_genes_volcano_plot <- function(countdata_group, condition_count) {
  #make the matrix name a variable for downstream saving of plots
  var_name <- deparse(substitute(countdata_group))
  
  # Assign conditions
  #this is setting the reference level to the mutant group
  #this means any log fold changes is how the mutants appear compared to the wildtype virus
  #e.g. gene x has a +2 log fold change, so that gene is expressed MORE in the tested mutant
  (condition <- factor(c(rep("Mutant", condition_count), rep("WildType", condition_count))))
  
  # Analysis with DESeq2 ----------------------------------------------------
  library("DESeq2")
  # Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
  (coldata <- data.frame(row.names=colnames(countdata_group), condition))
  dds <- DESeqDataSetFromMatrix(countData=countdata_group, colData=coldata, design=~condition)
  dds
  
  # Run the DESeq pipeline
  dds <- DESeq(dds)
  
  # Regularized log transformation for clustering/heatmaps, etc
  rld <- rlogTransformation(dds)
  head(assay(rld))
  hist(assay(rld))
  
  # Colors for plots below
  library(RColorBrewer)
  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
  
  # Get differential expression results
  res <- results(dds)
  table(res$padj<0.05)
  ## Order by adjusted p-value
  res <- res[order(res$padj), ]
  
  ## Merge with normalized count data
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  head(resdata)
  
  ## Examine plot of p-values
  hist(res$pvalue, breaks=50, col="grey")
  
  #volcano plot using ggplot2 (these plots are easily modified)
  res_df <- as.data.frame(res)
  
  #make the row names its own column
  library(data.table)
  setDT(res_df, keep.rownames = "Gene")[]
  
  #not all padj values have a number...these are just highly insignificant
  #replace NAs with .99 for downstream plotting
  res_df$padj[is.na(res_df$padj)] <- 0.99
  
  #only pull out the 212 shared genes
  #these genes were calculated above
  res_df_subset <- semi_join(res_df, shared_upreg_genes, by = c("Gene" = 'Ensembl gene id'))
  
  #add a threshold column for proper coloring of dots
  #likely overly complex use of ifelse statements...
  res_df_subset$threshold <- ifelse(res_df_subset$padj<0.05, 
                             ifelse(abs(res_df_subset$log2FoldChange)>1.5,"both", "pval"), 
                             ifelse(abs(res_df_subset$log2FoldChange)>1.5,"lfc", "none"))
  View(res_df_subset)
  
  p_volcano <- ggplot(res_df_subset) +
    geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
    ggtitle(paste0("Volcano Plot ", var_name)) +
    xlab("Log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    scale_color_manual(values = c('#009E73', '#D55E00', 'black', '#E69F00')) +
    #scale_y_continuous(limits = c(0,50)) +
    theme_bw() +
    xlim(c(-5,5)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))  
  
  
  #p_volcano
  ggsave(paste0("shared_genes-volcanoplot_", var_name, ".pdf"), plot = p_volcano, width = 2.5, height = 2.5, units = "in")
}

shared_genes_volcano_plot(WTv287, condition_count = 3)
shared_genes_volcano_plot(WTv296, condition_count = 3)
