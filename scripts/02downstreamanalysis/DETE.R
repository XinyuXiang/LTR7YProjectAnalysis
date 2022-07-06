######################################################################
###         script for differential expressed TE calling           ###
######################################################################

##### setup #####
library("clusterProfiler")
library(org.Hs.eg.db)
library(DESeq2)

##### input count #####
configure_TE=read.delim("configure.txt",header=FALSE)
read.l=list()
count.l=list()
for (i in 1:nrow(configure_TE)){
  read.l[[i]]=read.table(as.character(configure_TE[i,1]), comment.char = "#",header=T)
  count.l[[i]]=read.l[[i]][,7] 
}
read.df=as.data.frame(do.call(cbind,count.l))
names(read.df)=paste0(configure_TE[,2],"_",configure_TE[,3])
rownames(read.df)=rownames(read.l[[1]])

read.df=data.frame(Dupid = TE.df[,7], TEid = TE.df[,8], family = TE.df[,9], class = TE.df[,10], read.l[[1]][,2:6],read.df)
head(read.df)

count.df=read.df[which(rowSums(read.df[,-c(1:9)])>0),]
summary(count.df$Chr)
write.csv(count.df,"../DETE/count_onlyTEwithReads.csv",quote=F,row.names = F) 



##### calculate RPKM #####
con=configure_TE[,3]
exonLen=count.df$Length

RPKM.df =round(count.df[,10:(10+nrow(configure_TE)-1)]/exonLen*10^3 ,3)
totalreads=apply(count.df[,10:(10+nrow(configure_TE)-1)],2,sum)

for (i in 1:ncol(RPKM.df)){
  RPKM.df[,i]=round(RPKM.df[,i]/totalreads[i]*10^6,3)}

RPKM.df=data.frame(RPKM.df,apply(RPKM.df[,c(1:3)],1,mean),
                   apply(RPKM.df[,c(4:6)],1,mean),
                   apply(RPKM.df[,c(7:9)],1,mean),
                   apply(RPKM.df[,c(10:12)],1,mean))
names(RPKM.df)=c(names(RPKM.df)[1:nrow(configure_TE)],paste0(unique(con),"_mean"))
RPKM.df = cbind(count.df[,1:9], RPKM.df)
head(RPKM.df)
write.csv(RPKM.df,"../DETE/RPKM_onlyTEwithReads.csv",quote=F,row.names = F) 


##### extract basemean ##### 
colData=data.frame(condition=configure_TE[,3],
                   type=configure_TE[,2],
                   row.names=names(read.df[,-c(1:9)]))

dds = DESeqDataSetFromMatrix(count.df[,-c(1:9)],colData = colData, design = ~ condition)
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

write.csv(cbind(RPKM.df[,1:9], basemean.df), "../DETE/BaseMean_onlyTEwithReads.csv",quote=F,row.names = F)



######## DETE ######## 
# FC = A(former)/B(latter)
configure_TE=read.delim("configure.txt",header=FALSE)

idx=1:length(levels(configure_TE[,3]))
type=levels(con)
for (i in idx){
  for (j in idx[-i]){
    print(paste0("j=",j))
    print(paste0(levels(configure_TE[,3])[i],"_",levels(configure_TE[,3])[j]))
    DETE = results(dds, contrast=c("condition",type[i],type[j]))
    DETEout=data.frame(count.df[,1:9],DETE,
                       basemean.df[,which(names(basemean.df)==type[i])],
                       basemean.df[,which(names(basemean.df)==type[j])],
                       RPKM.df[,which(names(RPKM.df)==paste0(type[i],"_mean"))],
                       RPKM.df[,which(names(RPKM.df)==paste0(type[j],"_mean"))])
    names(DETEout)=c(names(count.df[,1:9]),names(DETE),paste0(type[c(i,j)],"_baseMean"),paste0(type[c(i,j)],"_RPKM"))
    
    write.csv(DETEout,paste0("../DETE/",type[i],"_vs_",type[j],"_desq2TE_all.csv"),quote=FALSE,row.names = FALSE)
    #DETE FC1.5
    DETE15fc.dw=subset(DETEout,DETEout$log2FoldChange<=log2(1/1.5) & DETEout$padj<=0.05 & (DETEout[,ncol(DETEout)]>=1 | DETEout[,ncol(DETEout)-1]>=1))
    write.csv(DETE15fc.dw,paste0("../DETE/desq2TE_",type[i],"_",type[j],"_downin",type[i],"_FC1.5_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    DETE15fc.up=subset(DETEout,DETEout$log2FoldChange>=log2(1.5) & DETEout$padj<=0.05 & (DETEout[,ncol(DETEout)]>=1 | DETEout[,ncol(DETEout)-1]>=1))
    write.csv(DETE15fc.up,paste0("../DETE/desq2TE_",type[i],"_",type[j],"_upin",type[i],"_FC1.5_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    #DETE FC 2
    DETE2fc.dw=subset(DETEout,DETEout$log2FoldChange<=log2(1/2) & DETEout$padj<=0.05 & (DETEout[,ncol(DETEout)]>=1 | DETEout[,ncol(DETEout)-1]>=1))
    write.csv(DETE2fc.dw,paste0("../DETE/desq2TE_",type[i],"_",type[j],"_downin",type[i],"_FC2_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    DETE2fc.up=subset(DETEout,DETEout$log2FoldChange>=log2(2) & DETEout$padj<=0.05 & (DETEout[,ncol(DETEout)]>=1 | DETEout[,ncol(DETEout)-1]>=1))
    write.csv(DETE2fc.up,paste0("../DETE/desq2TE_",type[i],"_",type[j],"_upin",type[i],"_FC2_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    #DETE FC 4
    DETE4fc.dw=subset(DETEout,DETEout$log2FoldChange<=log2(1/4) & DETEout$padj<=0.05 & (DETEout[,ncol(DETEout)]>=1 | DETEout[,ncol(DETEout)-1]>=1))
    write.csv(DETE4fc.dw,paste0("../DETE/desq2TE_",type[i],"_",type[j],"_downin",type[i],"_FC4_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
    DETE4fc.up=subset(DETEout,DETEout$log2FoldChange>=log2(4) & DETEout$padj<=0.05 & (DETEout[,ncol(DETEout)]>=1 | DETEout[,ncol(DETEout)-1]>=1))
    write.csv(DETE4fc.up,paste0("../DETE/desq2TE_",type[i],"_",type[j],"_upin",type[i],"_FC4_padj0.05_RPKMgt1.csv"),quote=FALSE,row.names = FALSE)
  }
}

