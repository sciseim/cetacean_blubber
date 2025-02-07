# ..............................................................
# # "sampletype_dADRB3_vs_mADRB3"  
the.filename <- "Ttru-ADRB3_vs_Mmus-ADRB3"
# contrast: the column from the metadata that is used for the grouping of the samples 
contrast.temp<- c("sampletype", "dADRB3", "mADRB3")

res_temp <- results(dds, contrast=contrast.temp, alpha = 0.05)
#@ res_temp %>% data.frame() %>% View()
## Save the unshrunken results to compare
res_temp_unshrunken <- res_temp
## Summarize results
summary(res_temp_unshrunken, alpha = 0.05)
# Apply fold change shrinkage
library(apeglm)
res_temp_shrink <- lfcShrink(dds = dds,
                             coef="sampletype_dADRB3_vs_mADRB3",
                             type="apeglm")
# same for all 
source("../filter_et_annotate.R")
Ttru.ADRB3_vs_Mmus.ADRB3 <- test2
# ..............................................................


# ..............................................................
# "sampletype_control_vs_mADRB3" 
the.filename <- "NC_vs_Mmus-ADRB3"
# contrast: the column from the metadata that is used for the grouping of the samples 
contrast.temp<- c("sampletype", "control", "mADRB3")

res_temp <- results(dds, contrast=contrast.temp, alpha = 0.05)
#@ res_temp %>% data.frame() %>% View()
## Save the unshrunken results to compare
res_temp_unshrunken <- res_temp
## Summarize results
summary(res_temp_unshrunken, alpha = 0.05)
# Apply fold change shrinkage
library(apeglm)
res_temp_shrink <- lfcShrink(dds = dds,
                             coef="sampletype_control_vs_mADRB3",
                             type="apeglm")
# same for all 
source("../filter_et_annotate.R")
NC_vs_Mmus.ADRB3 <- test2
# ..............................................................

# ..............................................................
# "sampletype_mADRB3.1d_vs_mADRB3"
the.filename <- "mADRB3.1d-ADRB3_vs_Mmus-ADRB3"
# contrast: the column from the metadata that is used for the grouping of the samples 
contrast.temp<- c("sampletype", "mADRB3.1d", "mADRB3")

res_temp <- results(dds, contrast=contrast.temp, alpha = 0.05)
#@ res_temp %>% data.frame() %>% View()
## Save the unshrunken results to compare
res_temp_unshrunken <- res_temp
## Summarize results
summary(res_temp_unshrunken, alpha = 0.05)
# Apply fold change shrinkage
library(apeglm)
res_temp_shrink <- lfcShrink(dds = dds,
                             coef="sampletype_mADRB3.1d_vs_mADRB3",
                             type="apeglm")
# same for all 
source("../filter_et_annotate.R")
mADRB3.1d.ADRB3_vs_Mmus.ADRB3 <- test2
# ..............................................................