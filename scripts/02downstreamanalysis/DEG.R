######################################################################
###        script for differential expressed gene calling          ###
######################################################################

##### setup #####
library("clusterProfiler")
library(org.Hs.eg.db)
library(DESeq2)

##### input count #####
configure_gene=read.delim("configure.txt",header=FALSE)
# DKOa_gene.featureCounts_PCG.txt        rep1        double_KO

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
write.csv(count.df,"../DEG/count_onlygenewithReads.csv",quote=F,row.names = F) 


##### calculate RPKM #####
con=configure_gene[,3]
exonLen=count.df$Length
RPKM.df=round(count.df[,7:(7+nrow(configure_gene)-1)]/exonLen*10^3 ,3)
totalreads=apply(count.df[,7:(7+nrow(configure_gene)-1)],2,sum)

for (i in 1:ncol(RPKM.df)){
  RPKM.df[,i]=round(RPKM.df[,i]/totalreads[i]*10^6,3)}

RPKM.df=data.frame(RPKM.df,apply(RPKM.df[,c(1:3)],1,mean),
                   apply(RPKM.df[,c(4:6)],1,mean),
                   apply(RPKM.df[,c(7:9)],1,mean),
                   apply(RPKM.df[,c(10:12)],1,mean))
names(RPKM.df)=c(names(RPKM.df)[1:nrow(configure_gene)],paste0(unique(con),"_mean"))
RPKM.df=data.frame(count.df[,1:6],RPKM.df)

gene.df <- bitr(RPKM.df$Geneid, fromType = "ENSEMBL", 
                toType = c("SYMBOL"), 
                OrgDb = org.Hs.eg.db,
                drop = F)
geneid = gene.df[match(RPKM.df$Geneid,gene.df$ENSEMBL),]
RPKM.df = data.frame(SYMBOL = geneid$SYMBOL,RPKM.df)
head(RPKM.df)
dim(RPKM.df)
write.csv(RPKM.df,"../DEG/RPKM_onlygenewithReads.csv",quote=F,row.names = F)


##### extract basemean ##### 
colData=data.frame(condition=configure_gene[,3],
                   type=configure_gene[,2],
                   row.names=names(read.df[,-c(1:6)]))

dds = DESeqDataSetFromMatrix(count.df[,-c(1:6)],colData = colData, design = ~ condition)
dds
dds <- DESeq(dds)
res=results(dds)
summary(res)

base.l=list()
basemean.l=list()
for (i in 1:length(unique(con))){
  base.l[[i]]=counts(dds, normalized=TRUE)[,colData(dds)$condition == unique(con)[i]]
  if (is.vector(base.l[[i]])){
    basemean.l[[i]] = base.l[[i]]
  } else {
    basemean.l[[i]] = rowMeans(base.l[[i]])
  }
}
basemean.df=as.data.frame(do.call(cbind,basemean.l))
names(basemean.df)=unique(con)
head(basemean.df)
dim(basemean.df)

write.csv(cbind(RPKM.df[,1:7], basemean.df), "../DEG/BaseMean_onlygenewithReads.csv",quote=F,row.names = F)


######## DEG ######## 
# FC = A(former)/B(latter)
configure_gene=read.delim("configure.txt",header=FALSE)

gene.df <- bitr(count.df$Geneid, fromType = "ENSEMBL", 
                toType = c("SYMBOL"), 
                OrgDb = org.Hs.eg.db,
                drop = F)

idx=1:length(levels(configure_gene[,3]))
type=levels(con)
for (i in idx){
  for (j in idx[-i]){
    print(paste0("j=",j))
    print(paste0(levels(configure_gene[,3])[i],"_",levels(configure_gene[,3])[j]))
    deg = results(dds, contrast=c("condition",type[i],type[j]))
    degout=data.frame(count.df[,1:6],deg,
                      basemean.df[,which(names(basemean.df)==type[i])],
                      basemean.df[,which(names(basemean.df)==type[j])],
                      RPKM.df[,which(names(RPKM.df)==paste0(type[i],"_mean"))],
                      RPKM.df[,which(names(RPKM.df)==paste0(type[j],"_mean"))])
    names(degout)=c(names(count.df[,1:6]),names(deg),paste0(type[c(i,j)],"_baseMean"),paste0(type[c(i,j)],"_RPKM"))
    ### add SYMBOL ID
    geneid = gene.df[match(degout$Geneid,gene.df$ENSEMBL),]
    degout = data.frame(SYMBOL = geneid$SYMBOL,degout)
    
    write.csv(degout,paste0("../DEG/",type[i],"_vs_",type[j],"_desq2gene_all.csv"),quote=FALSE,row.names = FALSE)
    #DEG FC1.5
    deg15fc.dw=subset(degout,degout$log2FoldChange<=log2(1/1.5) & degout$padj<=0.05 & (degout[,ncol(degout)]>=1 | degout[,ncol(degout)-1]>=1))
    write.csv(deg15fc.dw,paste0("../DEG/desq2gene_",type[i],"_",type[j],"_downin",type[i],"_FC1.5_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    deg15fc.up=subset(degout,degout$log2FoldChange>=log2(1.5) & degout$padj<=0.05 & (degout[,ncol(degout)]>=1 | degout[,ncol(degout)-1]>=1))
    write.csv(deg15fc.up,paste0("../DEG/desq2gene_",type[i],"_",type[j],"_upin",type[i],"_FC1.5_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    #DEG FC 2
    deg2fc.dw=subset(degout,degout$log2FoldChange<=log2(1/2) & degout$padj<=0.05 & (degout[,ncol(degout)]>=1 | degout[,ncol(degout)-1]>=1))
    write.csv(deg2fc.dw,paste0("../DEG/desq2gene_",type[i],"_",type[j],"_downin",type[i],"_FC2_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    deg2fc.up=subset(degout,degout$log2FoldChange>=log2(2) & degout$padj<=0.05 & (degout[,ncol(degout)]>=1 | degout[,ncol(degout)-1]>=1))
    write.csv(deg2fc.up,paste0("../DEG/desq2gene_",type[i],"_",type[j],"_upin",type[i],"_FC2_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    #DEG FC 4
    deg4fc.dw=subset(degout,degout$log2FoldChange<=log2(1/4) & degout$padj<=0.05 & (degout[,ncol(degout)]>=1 | degout[,ncol(degout)-1]>=1))
    write.csv(deg4fc.dw,paste0("../DEG/desq2gene_",type[i],"_",type[j],"_downin",type[i],"_FC4_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    deg4fc.up=subset(degout,degout$log2FoldChange>=log2(4) & degout$padj<=0.05 & (degout[,ncol(degout)]>=1 | degout[,ncol(degout)-1]>=1))
    write.csv(deg4fc.up,paste0("../DEG/desq2gene_",type[i],"_",type[j],"_upin",type[i],"_FC4_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
  }
}
