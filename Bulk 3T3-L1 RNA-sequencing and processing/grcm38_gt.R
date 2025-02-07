# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# let us obtain gene symbols
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(biomaRt)
# Mouse build 38 tx
remove(grcm38_gt)
grcm38_gt <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="mmusculus_gene_ensembl") %>% 
  getBM(mart=., attributes=c("ensembl_gene_id", "ensembl_transcript_id","mgi_symbol","description","external_gene_name"))
head(grcm38_gt)
tail(grcm38_gt)

# remove [Source:MGI Symbol
tail(grcm38_gt)
grcm38_gt$description <- gsub("\\[S.*","",grcm38_gt$description) # remove [Source:MGI
# remove trailing whitespace
grcm38_gt$description <- trimws(grcm38_gt$description)
#@ tail(grcm38_gt)
#@ View(grcm38_gt)

# get rid of transcripts and only keep one gene row ... not using ensembl_transcript_id
# get rid of ensembl_transcript_id 
grcm38_gt <- subset(grcm38_gt, select=-c(ensembl_transcript_id))
#
grcm38_gt <- subset(grcm38_gt, !duplicated(ensembl_gene_id))
#
saveRDS(grcm38_gt, "../grcm38_gt.rds") # ~3 Mb
grcm38_gt <- readRDS("../grcm38_gt.rds")



