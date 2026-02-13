#-------------------------------------------------------DESeq2 (Differential Gene Expression Analysis)---------------------------------------------------

#Vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/ins t/doc/DESeq2.html#standard-workflow
#Data: Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications (GSE52778)

#Study Design
 #RNA-sequencing performed on 4 primary human airway smooth muscle cell lines treated with 1 micromolar dexamethasone for 18 hours.
 #For each of the 4 cell lines, study has a treated and untreated sample.
 #Goal: To understand the transcriptional changes occurring due to treatment with dexamethasone.

#Required R packages:
 #DESeq2 
 #tidyverse
 #airway

#SCRIPT

# script to perform differential gene expression analysis using DESeq2 package
setwd("D:/AMAN_PhD/Script/Analysis/DESeq2 ")

#Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("airway")
BiocManager::install("EnhancedVolcano")


# load libraries
library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)
library(EnhancedVolcano)

# script to get data from airway package

library(airway)

data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

# Step 1: preparing count data ----------------

# read in counts data
counts_data <- read.csv('D:/AMAN_PhD/Script/Analysis/counts_data.csv')
head(counts_data)


# read in sample info
colData <- read.csv('D:/AMAN_PhD/Script/Analysis/sample_info.csv')


# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level for referencen(by default use alphabetical order)
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res
#base mean= average of the normalize count
#log2FoldChange= change in expression of treated wrt untreated 
                 #positive value=upregulated gene (treated)
                 #negative value=downregulated gene(treated)
#lfcse= Standerd error estimate for log2FC 
#stat= wald test value 
#pvale= confidence of difference
#padj= corrected for multiple testing


# Explore Results ----------------

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

# e.g.: treated_4hrs, treated_8hrs, untreated

results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# MA(manhattan) plot 
plotMA(res)

#volcanoplot
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Treated vs Untreated",
  subtitle = "DESeq2",
  pointSize = 2,
  labSize = 3
)



#DATA SOURCES
#DESeq2 Vignette: 
  #https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow

#Link to Code:
  #Code to get data: https://github.com/kpatel427/YouTubeTutorials/blob/main/getData.R
  #DESeq2 workflow:  https://github.com/kpatel427/YouTubeTutorials/blob/main/runDESeq2.R

#Additional resources:
  #1. https://bioc.ism.ac.jp/packages/2.14/...
  #2. https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
  #3. https://lashlock.github.io/compbio/R_presentation.html
  #4. https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/




#-------------------------------------------------------DESeq2 (Differential Gene Expression Analysis)---------------------------------------------------
