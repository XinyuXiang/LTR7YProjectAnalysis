######################################################################
###         script for bar plot of top DETE percentage             ###
######################################################################

##### setup #####
library(grid)
library(ggplot2)
library(gridExtra) 

subclass.total = read.table("GRCh38_allTE_CopyNumber.txt",header=F) # awk '{print $10}' GRCh38_rmsk_TE.gtf | sort | uniq -c > GRCh38_allTE_CopyNumber.txt
colnames(subclass.total) = c("Copynumber","TEid")

configure_DETE=read.table("configure_DETE.txt",header=FALSE)
# 20200507desq2TE_WT_double_KO_downinWT_FC2_padj0.05_RPKMgt1.csv FC2 WT_double_KO_downinWT_FC2_padj0.05 WT double_KO        down
configure_DETE_dw = configure_DETE[which(configure_DETE$V6=="down"),]
configure_DETE_up = configure_DETE[which(configure_DETE$V6=="up"),]

##### calculate DETE percentage and plot #####
# Filtering parameters for top DETE barplot: (1) DETE FC>=2 (2) total TE copy >=80 (3) up or down DETE copy >=8

for (i in 1:nrow(configure_DETE_dw)){
  # get up/dw TE id and freq
  DETEout_dw.df=read.csv(as.character(configure_DETE_dw[i,1]),header=T)
  DETEout_up.df=read.csv(as.character(configure_DETE_up[i,1]),header=T)
  
  DETEout_dw.subclass.df=as.data.frame(table(DETEout_dw.df$TEid))
  DETEout_up.subclass.df=as.data.frame(table(DETEout_up.df$TEid))
  union.subclass=union(DETEout_dw.subclass.df[,1],DETEout_up.subclass.df[,1])
  
  
  compare.subclass.df=data.frame(union.subclass,
                                 DETEout_dw.subclass.df[match(union.subclass,DETEout_dw.subclass.df$Var1),2],
                                 DETEout_up.subclass.df[match(union.subclass,DETEout_up.subclass.df$Var1),2],
                                 subclass.total[match(union.subclass,subclass.total$TEid),1])
  names(compare.subclass.df)=c("TEclass","DETEout_dw","DETEout_up","totalcopies")
  
  # calculate up_copy/total_copy% dw_copy/total_copy% 
  compare.subclass.df[is.na(compare.subclass.df)]=0
  compare.subclass.df=data.frame(compare.subclass.df,
                                 DETEout_dwperc = compare.subclass.df$DETEout_dw / compare.subclass.df$totalcopies,
                                 DETEout_upperc = compare.subclass.df$DETEout_up / compare.subclass.df$totalcopies)
  compare.subclass.df[,-1]=round(compare.subclass.df[,-1],3)
  compare.subclass.df=data.frame(compare.subclass.df,diff=compare.subclass.df$DETEout_upperc-compare.subclass.df$DETEout_dwperc)
  compare.subclass.df=compare.subclass.df[order(-compare.subclass.df$diff),]
  
  #only kept TE subclass with more than 80 copies
  #only kept TE subclass with more than 8 copies of DETEout_up or DETEout_dw
  compare.subclass.df=compare.subclass.df[which(compare.subclass.df$totalcopies>=80),]
  compare.subclass.df=compare.subclass.df[which(compare.subclass.df$DETEout_up>=8 | compare.subclass.df$DETEout_dw>=8),]
  
  # select first 10 and last 10 TE
  plot=rbind(compare.subclass.df[1:10,],compare.subclass.df[(nrow(compare.subclass.df)-9):nrow(compare.subclass.df),])
  plot$TEclass=factor(plot$TEclass,levels=c(rev(as.character(plot$TEclass))))
  write.csv(plot, paste0("barplot_top10varTE_",configure_DETE_up[i,3],".txt"))
  
  ymax = max(plot$DETEout_dwperc, plot$DETEout_upperc)
  
  # plot
  col = colorRampPalette(c("dodgerblue4", 'indianred4'))(20) 
  g.mid=ggplot(plot,aes(x=1,y=TEclass))+geom_text(aes(label=TEclass))+
    ggtitle("")+
    ylab(NULL)+
    theme(axis.title=element_blank(),
          panel.grid=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background=element_blank(),
          axis.text.x=element_text(color=NA),
          axis.ticks.x=element_line(color=NA),
          plot.margin = unit(c(1,-1,1,-1), "mm"))
  
  g1 = ggplot(data = plot, aes(x = TEclass, y = DETEout_upperc,fill=TEclass)) + 
    geom_bar(stat = "identity",show.legend = FALSE) + scale_fill_manual(values=col) + 
    ggtitle(paste0("up regulate TEs in ", configure_DETE_up[i,4], " (vs ", configure_DETE_up[i,5],") ", configure_DETE_up[i,2])) + 
    theme_bw()+
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_line(colour = NA),
          panel.grid.major = element_line(colour = NA),
          plot.margin = unit(c(1,-1,1,0), "mm")) + 
    scale_y_reverse(limits=c(ymax,0)) + coord_flip()
  
  g2 = ggplot(data = plot, aes(x = TEclass, y = DETEout_dwperc,fill=TEclass)) + 
    geom_bar(stat = "identity",show.legend = FALSE)+ scale_fill_manual(values=col) +
    ggtitle(paste0("down regulate TEs in ", configure_DETE_up[i,4], " (vs ", configure_DETE_up[i,5],") ", configure_DETE_up[i,2])) + theme_bw()+
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          plot.margin = unit(c(1,-1,1,0), "mm")) + 
    ylim(0,ymax) + coord_flip()
  
  gg1 = ggplot_gtable(ggplot_build(g1)) 
  gg2 = ggplot_gtable(ggplot_build(g2)) 
  gg.mid = ggplot_gtable(ggplot_build(g.mid)) 
  
  pdf(paste0("barplot_top10varTE_",configure_DETE_up[i,3],".pdf"),height=12,width=8)
  grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))
  dev.off()

}
