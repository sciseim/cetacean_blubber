# reset R  
rm(list=ls())
library("DESeq2")
library(dplyr)
library("tibble")

# turn OFF scientific notation
options(scipen=999)
# revert back
# options(scipen=0)

setwd("~/Dropbox/Manuscripts/--IN_PREPARATION--/Tong Zhang - ADRB3 loss in whales/--DEC 2024/ADRB3_iWAT/cell_salmon/salmon-counts/")

# create output dir
system("mkdir -p ../output") # may not use ... just output to ../





# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Load the tximport package that we use to import Salmon counts
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(tximport)
library("tximportData")
# tximport is just for importing data, it does not contain statistical tests for comparing samples.
samples <- list.files(path = "./", full.names = T, pattern="salmon$")
files <- file.path(samples, "quant.sf")
## Since all quant files have the same name it is useful to have names for each element
library(stringr)
names(files) <- str_replace(samples, ".//", "") %>% 
  str_replace(".salmon", "")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load annotations
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(readr)
library("GenomicFeatures")
#  BiocManager::install("txdbmaker")
library("txdbmaker")
#@ TxDb <- makeTxDbFromGFF(file = "../gencode.vM36.annotation.gtf")
#@ k <- keys(TxDb, keytype = "TXNAME")
#@ tx2gene <- select(TxDb, k, "GENEID", "TXNAME")
#@ head(tx2gene)
# 1 Gb GTF ... let us save tx2gene object
saveRDS(tx2gene, "../tx2gene.rds") # 1.6 Mb only 
tx2gene <- readRDS("../tx2gene.rds")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Run tximport ... 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("TXNAME", "GENEID")], countsFromAbundance="lengthScaledTPM")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# let us check some genes before normalizing for fun
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Write the counts to an object
data <- txi$counts %>% 
  round() %>% 
  data.frame()

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# metadata
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(openxlsx)
# load sample info
curated.dataset <- read.xlsx("../sampleinfo.xlsx")
head(curated.dataset)

## Create a sampletable/metadata
sampletype <- factor(curated.dataset$sampletype)
meta <- data.frame(sampletype, row.names = colnames(txi$counts))
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Count normalization with DESeq2
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
dds <- DESeqDataSetFromTximport(txi, colData = curated.dataset, design = ~ sampletype)
# 3. Generate the normalized counts
dds <- estimateSizeFactors(dds)
# take a look at the normalization factor applied to each sample
sizeFactors(dds)
#
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="../normalized_counts.txt", sep="\t", quote=F, col.names=NA)
#@ head(normalized_counts)



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Create DESeq2Dataset object
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
dds <- DESeqDataSetFromTximport(txi, colData = curated.dataset, design = ~ sampletype)
#@ class(dds) # "DESeqDataSet"

dds$sampletype <- relevel(dds$sampletype, ref = "mADRB3") # set MOUSE ADRB3 AS THE REFERENCE CONTRAST!!!!
# THIS STEP IS VERY IMPORTANT SINCE DESeq apparently looks for a "control/Control column" and would output incorrect DEG here if we do not force a design!

## Run analysis
dds <- DESeq(dds)



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# PULL OUT COMPARISONS
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# resultsNames(dds) # 
# [1] "Intercept"                      "sampletype_control_vs_mADRB3"  
# [3] "sampletype_dADRB3_vs_mADRB3"    "sampletype_mADRB3.1d_vs_mADRB3"
# can see that we have here correctly compared all to mADRB3
# ... we need these later for shrinking log-values AND othe rcontrasts!!!

# set DEG cut-offs here
padj.cutoff <- 0.05
lfc.cutoff <- 0.58 # is 2^0.58 = 1.5-fold absolute Fc
#@ lfc.cutoff <- 2 # is 2^2 = 4-fold absolute Fc
#@ lfc.cutoff <- 1.5 # is 2^1.5 = 2.8-fold absolute Fc
#@ lfc.cutoff <- 1.0 # is 2^1 = 2-fold absolute Fc

# loop through comparisons
source("../sample_loop.R")


# save output
intersect(intersect(Ttru.ADRB3_vs_Mmus.ADRB3$mgi_symbol,NC_vs_Mmus.ADRB3$mgi_symbol),mADRB3.1d.ADRB3_vs_Mmus.ADRB3$mgi_symbol)
the.overlap <- intersect(intersect(Ttru.ADRB3_vs_Mmus.ADRB3$mgi_symbol,NC_vs_Mmus.ADRB3$mgi_symbol),mADRB3.1d.ADRB3_vs_Mmus.ADRB3$mgi_symbol)
write.table(the.overlap, paste("../",lfc.cutoff,"fc","-overlap.txt",sep=""), row.names=F, col.names=F, sep="\n", quote=F)
# pvalue: Wald test p-value: Indicates whether the gene analysed is likely to be differentially expressed in that comparison. The lower the more significant.
# padj: Bonferroni-Hochberg adjusted p-values (FDR): the lower the more significant. More robust that the regular p-value because it controls for the occurrence of false positives.
