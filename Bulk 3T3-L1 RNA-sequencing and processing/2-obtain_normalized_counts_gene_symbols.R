# reset R  
rm(list=ls())


setwd("~/Dropbox/Manuscripts/--IN_PREPARATION--/Tong Zhang - ADRB3 loss in whales/--DEC 2024/ADRB3_iWAT/cell_salmon/salmon-counts/")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load annotations
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
grcm38_gt <- readRDS("../grcm38_gt.rds")
# head(grcm38_gt) # ensembl_gene_id ensembl_transcript_id mgi_symbol   description external_gene_name hgnc_symbol
# tail(grcm38_gt


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load normalised counts
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
normalized_counts <- read.table(file="../normalized_counts.txt", sep="\t", header=TRUE)
#@ head(normalized_counts)

# add a gene column for easier joining below
# normalized_counts$gene <- row.names(normalized_counts) # for some reason set to X ...
colnames(normalized_counts)[1] ="gene"
#@ head(normalized_counts)


# remove ensembl tx suffix
# e.g. ENSMUSG00000064336.123 to ENSMUSG00000064336
normalized_counts$gene <- gsub("\\..*", "", normalized_counts$gene)

# look-up
library(dplyr)
# colnames(test.df) # "gene"           "baseMean"       "log2FoldChange" "lfcSE"          "pvalue"         "padj"  
normalized_counts.symbols <- normalized_counts %>% 
  dplyr::arrange(gene) %>% 
  dplyr::inner_join(grcm38_gt, by = c("gene" = "ensembl_gene_id")) 
# %>% 
#  dplyr::select(gene, baseMean,log2FoldChange, lfcSE, pvalue, padj, external_gene_name,mgi_symbol, description) 
# View(normalized_counts.symbols)
# 9999999999999999999999999999999999999999999999999999999999999999999999
#@ head(normalized_counts.symbols)

# not sure why some repeats ... 
normalized_counts.symbols <- normalized_counts.symbols[!duplicated(normalized_counts.symbols), ]
#@ head(normalized_counts.symbols)

# write output
write.table(normalized_counts.symbols, file="../normalized_counts_symbols.txt", sep="\t", quote=F, col.names=NA)