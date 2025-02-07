rm(list=ls())  # clear environment

setwd("~/Dropbox/Manuscripts/--IN_PREPARATION--/Tong Zhang - ADRB3 loss in whales/--JAN 2025/--Figures workings/Fig4-hypothesis/SST/")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load data from Nature MS 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# we have TEX86 Bayspar -derived SST data 
library(openxlsx)
# load sample info
the.dataset <- read.xlsx("41586_2022_5017_MOESM2_ESM.xlsx",sheet=1)
# get rid of first row
the.dataset <- the.dataset[-c(1),]


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Prepare the data set
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# fix names
colnames(the.dataset)[colnames(the.dataset) == 'Bayspar.5th.(°C)'] <- "Bayspar.5th"
colnames(the.dataset)[colnames(the.dataset) == "Bayspar.95th.(°C)"] <- "Bayspar.95th"
colnames(the.dataset)[colnames(the.dataset) == "SST.Bayspar.(°C)"] <- "SST.Bayspar"
colnames(the.dataset)[colnames(the.dataset) == "Age.(Ma)"] <- "Mya"
colnames(the.dataset)[colnames(the.dataset) == "SSTH.(°C)"] <- "SSTH"
colnames(the.dataset)[colnames(the.dataset) == "SSTL.(°C)"] <- "SSTL" 
# keep only the columns we want
the.dataset.subset <- subset(the.dataset, select = c(   "Mya","Latitude.category","SSTL","SSTH", "Bayspar.5th","SST.Bayspar","Bayspar.95th") )

# let us grab high latitude, middle latitude, and low latitude only!
library(tidyverse)
library(data.table)
the.dataset.subset <- the.dataset.subset[which(the.dataset.subset$Latitude.category %like% "latitude"),]

#@ class(the.dataset.subset$Latitude.category) # factor
the.dataset.subset$Latitude.category <- as.factor(the.dataset.subset$Latitude.category) # required to plot
#
the.dataset.subset$SSTL <- as.numeric(the.dataset.subset$SSTL) # required to plot
the.dataset.subset$SSTH <- as.numeric(the.dataset.subset$SSTH) # required to plot
the.dataset.subset$Bayspar.5th <- as.numeric(the.dataset.subset$Bayspar.5th) # required to plot
the.dataset.subset$SST.Bayspar <- as.numeric(the.dataset.subset$SST.Bayspar) # required to plot
the.dataset.subset$Bayspar.95th <- as.numeric(the.dataset.subset$Bayspar.95th) # required to plot

# OK, let us merge low and middle latitude
# from Nature MS: Unfortunately, low-latitude sites are scarce and do not provide a continuous SST record across the Cenozoic. Middle-latitude records are more abundant and show SST trends and absolute values that are very similar to those from low latitudes in the intervals in which both records overlap. Consequently, we combined the middle/low latitudes SST records to calculate the latitudinal SST gradient. 
the.dataset.subset$Latitude.category.NatureMS <- the.dataset.subset$Latitude.category
length(which(the.dataset.subset$Latitude.category %like% "mid")) # 1248 records
the.dataset.subset$Latitude.category.NatureMS <- gsub("mid", "low", the.dataset.subset$Latitude.category.NatureMS)

# CUT OFF THE DATA SET TO FIT WITH STEM CETACEAN EMERGENCE
# cut off at a sane time! mya
low.cut.off.time <- 1 # Mya
the.dataset.subset <- the.dataset.subset[which(the.dataset.subset$Mya >= low.cut.off.time ),] # 3,215 data points for 1M
high.cut.off.time <- 53 # Mya
the.dataset.subset <- the.dataset.subset[which(the.dataset.subset$Mya <= high.cut.off.time ),] # 3,215 data points for 1M
max(the.dataset.subset$Mya) # 55.9985

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# plot
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(ggplot2)

# Nature MS: LOESS smoothing (factor 0.02) and 90% confidence interval 
the.smoothing.factor <- 0.02
# A smaller span (e.g. span = 0.5) results in more local (flexible) smoothing, while a larger span (e.g. span = 1.5) produces more global (smooth) smoothing
# ... The degree of smoothing is controlled by the span parameter in the geom_smooth() layer.
the.confidence.interval <- 0.90


# let us flip the x-axis so we start with older
# ADDING DOTTED LINES HERE RESULTS IN A GRAPH THAT IS HARD TO READ ... better to manually indicate epochs!
# and add MINUS to put the dotted lines in here since we are reversing!
the.height <- 5
the.width <- 10
# save output
pdf("SST.pdf", height=the.height, width=the.width)
h <- ggplot(data = the.dataset.subset, aes(Mya, SST.Bayspar, color = Latitude.category.NatureMS))  + 
  geom_point(alpha = 1/8, show.legend = FALSE) + 
  geom_smooth(weight=5, method = "loess",  level=the.confidence.interval, span = the.smoothing.factor,linewidth=2 ) + 
  labs(x = "Mya", y = "SST (°C)") + 
  theme_classic () +  scale_y_continuous(breaks = seq(5, 45, by = 5)) + theme(ggh4x.axis.ticks.length.minor = rel(2)) + scale_color_manual(values=c("#3C5488B2", "#E64B35B2")) + coord_cartesian(ylim = c(10, 40)) + scale_x_reverse(guide=guide_axis(minor.ticks = TRUE),breaks=c(53,50,40,30,25,20,15,10,1))
print(h)
dev.off()
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



# Let us draw low-latitude only
# but, high-latitude data set records limited and highly variable for the time-periods we are interested in! Let us graph low-latitude only
# actually nicer with grey + dot transparent
low.only.df <- the.dataset.subset[ which(the.dataset.subset$Latitude.category.NatureMS %like% "low") ,]
#
pdf("low__SST.pdf", height=the.height, width=the.width)
h <- ggplot(data = low.only.df, aes(Mya, SST.Bayspar, color = Latitude.category.NatureMS))  + 
  geom_point(alpha = 1/8, show.legend = FALSE) + 
  geom_smooth(weight=5, method = "loess",  level=the.confidence.interval, span = the.smoothing.factor,linewidth=2 ) + 
  labs(x = "Mya", y = "SST (°C)") + 
  theme_classic () +  scale_y_continuous(breaks = seq(5, 45, by = 5)) + theme(ggh4x.axis.ticks.length.minor = rel(2)) + scale_color_manual(values=c("#E64B35B2", "#E64B35B2")) + coord_cartesian(ylim = c(10, 40)) + scale_x_reverse(guide=guide_axis(minor.ticks = TRUE),breaks=c(53,50,40,30,25,20,15,10,1))
print(h)
dev.off()