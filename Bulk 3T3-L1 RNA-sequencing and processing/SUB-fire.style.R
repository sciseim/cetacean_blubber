# fire style
cols <- paletteer_d("beyonce::X34", type="continuous")



pdf(the.file.name, height=the.height, width=the.width)
h <- heatmap.3((them.counts.heat.matrix.subset), col=cols, scale="row", trace="none",density.info="none", dendrogram="none",Colv=FALSE, key=TRUE,cexRow=1.5, hclustfun=myclust, distfun=mydist,labRow = rownames(them.counts.heat.matrix.subset),ColSideColors=clab,ColSideColorsSize=1.5 )
dev.off()


