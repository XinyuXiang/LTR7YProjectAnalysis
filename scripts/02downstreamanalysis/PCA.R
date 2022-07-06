######################################################################
###                     script for PCA plot                        ###
######################################################################

##### setup #####
library(DESeq2)
library(ggplot2)


####### PCA of genes #######
##### input count #####
configure_gene=read.delim("configure_gene.txt",header=FALSE)
# TesR_1.featureCounts_PCG.txt        rep1        WT_Primed

read.l=list()
count.l=list()
for (i in 1:nrow(configure_gene)){
  read.l[[i]]=read.table(as.character(configure_gene[i,1]), comment.char = "#",header=T)
  count.l[[i]]=read.l[[i]][,7] 
}
read.df=as.data.frame(do.call(cbind,count.l))
names(read.df)=paste0(configure_gene[,2],"_",configure_gene[,3])
rownames(read.df)=rownames(read.l[[1]])
read.df=data.frame(read.l[[1]][,1:6],read.df)
head(read.df)
dim(read.df)

count.df=read.df[which(rowSums(read.df[,-c(1:6)])>0),]
count.df = count.df[-(grep(pattern = "^GL|^KI|^MT", count.df$Chr)),] #delete scaffold and MT
summary(count.df$Chr)
dim(count.df)


##### calculate PCA #####
colData=data.frame(condition=configure_gene[,3],
                   type=configure_gene[,2],
                   row.names=names(read.df[,-c(1:6)]))


dds = DESeqDataSetFromMatrix(count.df[,-c(1:6)],colData = colData, design = ~ condition)
dds

rld <- rlogTransformation(dds)   # normalization
exprSet_new=assay(rld)   #extract normalized data 
head(exprSet_new)
#write.table(exprSet_new, file="gene_normalization.txt", sep="\t",quote=F)

pcaData <- plotPCA(rld, intgroup=c('condition'), returnData=TRUE)
#write.table(pcaData, "gene_pcaData.txt")
percentVar <- round(100 * attr(pcaData, "percentVar"))


##### plot PCA #####
pdf(paste0("PCA_DESeq2_gene.pdf"),height = 3,width = 4)
p=ggplot(data = pcaData,aes(x=PC1,y=PC2,colour=condition))
p+geom_point(alpha=0.4, size=1.5)+ #geom_text(hjust = 0.5, nudge_y = 1)+
  theme_bw()+
  scale_colour_manual(values=col1)+
  theme(panel.grid.major = element_line(colour =  NA))+
  theme(panel.grid.minor = element_line(colour =  NA))+
  theme(legend.title=element_blank())+
  labs(x = paste0("PC1 ",percentVar[1],"%"),
       y = paste0("PC2 ",percentVar[2],"%"),
       title="PCA")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()







####### PCA of TEs #######
##### input count #####
configure_TE=read.delim("configure_TE.txt",header=FALSE)
# TesR_1.featureCounts_exon_TEcy.txt        rep1        WT_Primed

read.l=list()
count.l=list()
for (i in 1:nrow(configure_gene)){
  read.l[[i]]=read.table(as.character(configure_gene[i,1]), comment.char = "#",header=T)
  count.l[[i]]=read.l[[i]][,7] 
}
read.df=as.data.frame(do.call(cbind,count.l))
names(read.df)=paste0(configure_gene[,2],"_",configure_gene[,3])
rownames(read.df)=rownames(read.l[[1]])
read.df=data.frame(read.l[[1]][,1:6],read.df)
head(read.df)
dim(read.df)

count.df=read.df[which(rowSums(read.df[,-c(1:6)])>0),]
count.df = count.df[-(grep(pattern = "^GL|^KI|^MT", count.df$Chr)),] #delete scaffold and MT
summary(count.df$Chr)
dim(count.df)


##### calculate PCA #####
colData=data.frame(condition=configure_gene[,3],
                   type=configure_gene[,2],
                   row.names=names(read.df[,-c(1:6)]))


dds = DESeqDataSetFromMatrix(count.df[,-c(1:6)],colData = colData, design = ~ condition)
dds

rld <- rlogTransformation(dds)   # normalization
exprSet_new=assay(rld)   #extract normalized data 
head(exprSet_new)
#write.table(exprSet_new, file="TE_normalization.txt", sep="\t",quote=F)

pcaData <- plotPCA(rld, intgroup=c('condition'), returnData=TRUE)
#write.table(pcaData, "TE_pcaData.txt")
percentVar <- round(100 * attr(pcaData, "percentVar"))


##### plot PCA #####
pdf(paste0("PCA_DESeq2_TE.pdf"),height = 3,width = 4)
p=ggplot(data = pcaData,aes(x=PC1,y=PC2,colour=condition))
p+geom_point(alpha=0.4, size=1.5)+ #geom_text(hjust = 0.5, nudge_y = 1)+
  theme_bw()+
  scale_colour_manual(values=col1)+
  theme(panel.grid.major = element_line(colour =  NA))+
  theme(panel.grid.minor = element_line(colour =  NA))+
  theme(legend.title=element_blank())+
  labs(x = paste0("PC1 ",percentVar[1],"%"),
       y = paste0("PC2 ",percentVar[2],"%"),
       title="PCA")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



