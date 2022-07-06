######################################################################
###    script for heatmap of selected gene/TE expression pattern   ###
######################################################################

##### setup #####
library(ComplexHeatmap)

##### input count #####
### choose genes to plot
genelist_total = read.table("genelist.txt")

##### filt genes #####
gene.RPKM.df = read.csv("RPKM_onlygenewithReads.csv")

filt.gene.RPKM.df = gene.RPKM.df[match(genelist_total, gene.RPKM.df$SYMBOL),]
row.names(filt.gene.RPKM.df) = genelist_total
dim(filt.gene.RPKM.df)


##### plot by z-score #####
plot.df = filt.gene.RPKM.df[,8:19]
plot.df = plot.df[which(rowSums(plot.df)>0),]
zscore=function(x){
  (x-mean(x))/sd(x)
}
heatmapzscore.df=apply(plot.df,1,zscore)
heatmapzscore.df=as.data.frame(heatmapzscore.df)
heatmapzscore.df=t(heatmapzscore.df)
max(heatmapzscore.df)
min(heatmapzscore.df)

bk1 <- c(seq(floor(min(heatmapzscore.df)),0,by=0.01))
bk2 <- c(seq(0.01,ceiling(max(heatmapzscore.df)),by=0.01))
col1 = c(colorRampPalette(colors = c('#4393C3',"white"))(length(bk1)),
         colorRampPalette(colors = c("white","#D6604D"))(length(bk2)))

pdf("../heatmap/heatmap_Zscore.pdf",height=6,width=6)
Heatmap(heatmapzscore.df,
        col = col1,
        cluster_columns = T,
        cluster_rows = F,
        show_row_names = T
        )
dev.off()


##### plot by log2(RPKM+1) #####
plot.df = filt.gene.RPKM.df[,8:19]
plot.df = plot.df[which(rowSums(plot.df)>0),]
heatmaplog.df=log2(plot.df+1)

max(heatmaplog.df) 
min(heatmaplog.df) 


bk3 <- c(seq(0,1,by=0.01))
bk4 <- c(seq(1.01,ceiling(max(heatmapzscore.df)),by=0.01))
col2 = c(colorRampPalette(colors = c('#4393C3',"white"))(length(bk3)),
         colorRampPalette(colors = c("white","#D6604D"))(length(bk4)))                     

pdf("../heatmap/heatmap_log2RPKMplus1.pdf",height=6,width=6)
Heatmap(heatmaplog.df,
        col = col2,
        cluster_columns = T,
        cluster_rows = F,
        show_row_names = T
        )
dev.off()