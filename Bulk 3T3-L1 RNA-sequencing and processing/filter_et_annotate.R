# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Summarize results
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
summary(res_temp_shrink, alpha = 0.05)
# check first rows of both results
head(res_temp_unshrunken)
head(res_temp_shrink)

library(tibble)
library(dplyr)
res_temp_shrink_tb <- res_temp_shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#@ padj.cutoff <- 0.05
# lfc.cutoff <- 0.58 # is 2^0.58 = 1.5-fold absolute Fc
#@ lfc.cutoff <- 2 # is 2^4 = 16-fold absolute Fc
sig.res_temp_shrink <- res_temp_shrink_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

#@ class(sig.res_temp_shrink) # [1] "tbl_df"     "tbl"        "data.frame"
#@ head(sig.res_temp_shrink)

# assign to new df
test.df <- as.data.frame(sig.res_temp_shrink)


# 9999999999999999999999999999999999999999999999999999999999999999999999
library(dplyr)
#@ grcm38_gt <- readRDS("../grcm38_gt.rds")
# head(grcm38_gt) # ensembl_gene_id ensembl_transcript_id mgi_symbol   description external_gene_name hgnc_symbol
# tail(grcm38_gt)

# remove ensembl tx suffix
# e.g. ENSMUSG00000064336.123 to ENSMUSG00000064336
test.df$gene <- gsub("\\..*", "", test.df$gene)

# colnames(test.df) # "gene"           "baseMean"       "log2FoldChange" "lfcSE"          "pvalue"         "padj"  
test2 <- test.df %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grcm38_gt, by = c("gene" = "ensembl_gene_id")) %>% 
  dplyr::select(gene, baseMean,log2FoldChange, lfcSE, pvalue, padj, external_gene_name,mgi_symbol, description) 
# View(test2)
# 9999999999999999999999999999999999999999999999999999999999999999999999
# write output
write.table(test2, paste("../",the.filename,"_",lfc.cutoff,"fc","-.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)