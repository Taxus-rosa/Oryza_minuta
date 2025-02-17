setwd(".../.../...")
gdata <- read.csv("DEG_list.csv", header=T, row.names=1)
library(pheatmap)
bk <- c(seq(-12,0,by=0.01),seq(0.01,12,by=0.01))
pheatmap(gdata,color = c(colorRampPalette(colors = c("#00008B","White"))(length(bk)/2),
         colorRampPalette(colors = c("White","#FF4500"))(length(bk)/2)),
         legend_breaks=seq(-12,12,5),breaks=bk, cluster_cols=T,cluster_rows=T,cutree_rows = 5,
         treeheight_row = 100, treeheight_col = 50, cellwidth = 30,show_rownames = F,
		 filename="DEG_list_log2.cut6.pdf")

scale_rows <- function (x) {
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}

#gmat = switch('row', none = gdata, row = scale_rows(gdata), column = t(scale_rows(t(gdata)))) ###normalization

gd = dist(gdata, method = 'euclidean')
gtree = hclust(gd, method = 'complete')
cluster = cutree(gtree,6)
#ggaps = which((cluster[-1] - cluster[-length(gv)]) != 0)
gene.cluster <- as.data.frame(cluster)
write.csv(gene.cluster, file="DEG_list_log2.cut6.cluster.csv") 

### set colors
ancols.Cut <- c( "#E41A1C", "#4169E1", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33") 
names(ancols.Cut) <- unique(gene.cluster$cluster)
ann_colors <- list(cluster = ancols.Cut)
				   
pheatmap(gdata,color = c(colorRampPalette(colors = c("#00008B","White"))(length(bk)/2),
         colorRampPalette(colors = c("White","#FF4500"))(length(bk)/2)),  
         legend_breaks=seq(-15,15,5),breaks=bk, cluster_cols=T,cluster_rows=T,cutree_rows = 6,
         treeheight_row = 100, treeheight_col = 50, cellwidth = 30,show_rownames = F,
		 annotation_row = gene.cluster, angle_col = "45", annotation_colors=ann_colors,
		 filename="DEG_list_log2.cut6.pdf")
