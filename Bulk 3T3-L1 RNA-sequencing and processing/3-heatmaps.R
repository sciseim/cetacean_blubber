# reset R  
rm(list=ls())

setwd("~/Dropbox/Manuscripts/--IN_PREPARATION--/Tong Zhang - ADRB3 loss in whales/--DEC 2024/ADRB3_iWAT/cell_salmon/salmon-counts/")


# load normalized counts
normalized_counts.symbols <- read.table(file="../normalized_counts_symbols.txt", sep="\t", header=TRUE,quote="")
# quote="", otherwise it will break if commas in a column!

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load the 3T3-L1 DEGs that overlap
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
the.overlap <- read.table(file="../0.58fc-overlap.txt", sep="\t", header=FALSE)
head(the.overlap)
the.overlap <- as.character(the.overlap$V1) # otherwise subsetting will fail!
head(the.overlap)
class(the.overlap)
them.counts <- normalized_counts.symbols[normalized_counts.symbols$mgi_symbol %in% the.overlap,]
#@ View(them.counts)

them.counts.heat <- them.counts
#@ row.names(them.counts.heat) <- them.counts.heat$mgi_symbol # FAILS!
# ... non-unique values when setting 'row.names': ‘Ddit3’, ‘Gcat’, ‘Pakap’ 
them.counts[which(them.counts$mgi_symbol=="Ddit3"),]
them.counts[which(them.counts$mgi_symbol=="Gcat"),]
them.counts[which(them.counts$mgi_symbol=="Pakap"),]
# hmmm, three overlapping genes have ~2-3 ENSEMBL gene IDs vs. symbols ... I guess just remove these from the heatmap. Not any of the ones of interest
# 28871	ENSMUSG00000025408	677.991186461835	620.83596361903	835.944885132939	2805.4673968513	1925.09937199302	3513.4670287027	1560.12217207744	2354.37505988043	2095.30272644778	3843.18638757527	2232.50282864609	2837.89195936205	Ddit3	DNA-damage inducible transcript 3	Ddit3
# 101680	ENSMUSG00000116429	0	0	0	0	0	0	0	0	0	0	0	0	Ddit3	DNA-damage inducible transcript 3	Ddit3
# let us remove ALL mgi_symbol duplicates ... 
length(them.counts.heat$mgi_symbol) # 380
them.counts.heat <- them.counts.heat[!duplicated(them.counts.heat$mgi_symbol), ]
length(them.counts.heat$mgi_symbol) # 376
# try again
row.names(them.counts.heat) <- them.counts.heat$mgi_symbol # success!
# 
head(them.counts.heat)
them.counts.heat <- them.counts.heat[order(them.counts.heat$mgi_symbol), ] #sort data frame column alphabetically
head(them.counts.heat)
# remove columns we do not want
colnames(them.counts.heat)
them.counts.heat <- subset(them.counts.heat, select = -c(gene, mgi_symbol, description,external_gene_name))



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Time to draw the heat map
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
them.counts.heat.matrix <- as.matrix(them.counts.heat)
#
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# install a nice pallete 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
#@ install.packages("ggsci")
#@ install.packages("gridExtra")
#@ library("ggplot2")
library(ArchR)
library(paletteer)





# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# sample info
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(openxlsx)
# load sample info
curated.dataset <- read.xlsx("../sampleinfo.xlsx")
head(curated.dataset)
sampletype <- factor(curated.dataset$tissuecolour) 
sampletype <- curated.dataset$tissuecolour
clab <- (cbind(sampletype))
clab <- as.matrix(clab)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 376 overlapping genes
# myclust=function(c) {hclust(c,method="average")}
#@ cols <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(256))
#@ cols<- colorRampPalette(c("red", "white", "blue"))(256)
# try to transpose instead
the.file.name <- paste("~/Downloads/test.pdf",sep="")
#@@@@ the.file.name <- paste("../","3T3_overlap","_DEG_heatmap.pdf",sep="")
mergedMATRIX <- them.counts.heat.matrix


# if continous palletes
# add a break so that colour is clearer
number.of.colors <- 50
cols <- paletteer_c("pals::kovesi.linear_bmy_10_95_c78",n = number.of.colors)
myBreaks <- seq(-2, 2, length.out=(number.of.colors+1))
#

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load heatmap.3
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
source("../heatmap.3.R")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# heat map of all overlapping DEGs
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# plot
pdf(the.file.name, height=10, width=10)
h <- heatmap.3((mergedMATRIX), col=cols, scale="row", trace="none",density.info="none", dendrogram="col",Colv=TRUE, key=TRUE,cexRow=1.5, hclustfun=myclust, distfun=mydist,labRow = FALSE,ColSideColors=clab,ColSideColorsSize=1,breaks=myBreaks )
dev.off()
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# heat map of themes
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# make them smaller
the.height <- 4
the.width <- 4

# thermogenesis
genes.of.interest <- c("Cmpk2","Trib1")
subset.title <- "thermogenesis"
the.file.name <- paste("../",subset.title,"3T3_overlap","_DEG_heatmap.pdf",sep="")
# for each
them.counts.heat.matrix.subset <- as.matrix ( them.counts.heat.matrix[which(row.names(them.counts.heat.matrix) %in% genes.of.interest),] )
class(them.counts.heat.matrix.subset) # matrix array
head(them.counts.heat.matrix.subset)
source("../SUB_npg-style.R") # npg style
#@ source("../SUB-fire.style.R") # fire-style 


