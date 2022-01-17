## Single AA mutations clustering


## required packages
require(ggplot2)
require(stringr)
require(reshape2)
require(dplyr)
require(grid)
require(ggdendro)


dir.create("_clustering")
path="_clustering"


##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")



#########################################################################################################################################

AB_wt<-"DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
ABseq<-unlist(strsplit(AB_wt, ""))
ABseq_pos<-paste0(ABseq, "\n", c(1:42))

all_aa=c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")



## NS color gradient 
min<-abs(min(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))
max<-abs(max(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))

cols <- c(colorRampPalette(c( "brown3", "grey95"))((min/(min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(min+max)*100)-0.5))







# clustering substitutions, insertions (considered before or after) and deletions at that position, by potision and by AA

#regions<-c("NT", "CT", "all")
regions<-c("all")
before_after<-c("before", "after")


for(bef_af in before_after){
  for(reg in regions){
  
  insertions_before_after<-bef_af
  
    
    clust_subs_df<-reshape(singles.df[,c("Pos", "Mut", "nscore_c")], idvar=c('Pos'), timevar='Mut', direction='wide')
    colnames(clust_subs_df)[2:length(clust_subs_df)]<- paste0("s", unique(singles.df$Mut))
    
    
    # insertions after
    clust_ins_df<-reshape(insertions_reps[,c("ins_pos", "ins_aa", "nscore_c")], idvar=c('ins_pos'), timevar='ins_aa', direction='wide')
    #insertions before
    if(insertions_before_after=="before"){ clust_ins_df$ins_pos<-as.numeric(clust_ins_df$ins_pos)+1}
    
    
    
    colnames(clust_ins_df)[2:length(clust_ins_df)]<- paste0("i", unique(insertions_reps$ins_aa))
    clust_ins_df$ins_pos<-as.numeric(clust_ins_df$ins_pos)
    
    merged_clust<-left_join(clust_subs_df, clust_ins_df, by=c("Pos"="ins_pos"))
    
    clust_deletions<-single_deletions_reps[,c("nscore_c", "del_pos")]
    clust_deletions$del_pos<-as.numeric(clust_deletions$del_pos)
    
    merged_clust<-left_join(merged_clust, clust_deletions[,c("nscore_c", "del_pos")], by=c("Pos"="del_pos"))
    colnames(merged_clust)[length(merged_clust)]<-"del"
    
    
    if(reg=="NT"){merged_clust<-merged_clust[merged_clust$Pos %in% c(1:28),]}
    if(reg=="CT"){merged_clust<-merged_clust[merged_clust$Pos %in% c(29:42),]}
    
    
    #make the positions the row numbers so I can see them in the dendogram
    rownames(merged_clust)<-merged_clust$Pos
    
    dendro_df <- as.dendrogram(hclust(d = dist(x = merged_clust)))
    dendro_plot <- ggdendrogram(data = dendro_df)
    
    print(dendro_plot)
    
    new_order<-order.dendrogram(dendro_df)
    levels=rownames(merged_clust)[new_order]
    
    
    
    # the same for the AA
    merged_clust_rev<-as.data.frame(t(merged_clust))
    merged_clust_rev<-merged_clust_rev[-1,]
    
    rev_dendro_df <- as.dendrogram(hclust(d = dist(x = merged_clust_rev)))
    dendro_plot_rev <- ggdendrogram(data = rev_dendro_df, rotate = T)
    
    print(dendro_plot_rev)
    
    new_order_aa<-order.dendrogram(rev_dendro_df)
    levels_aa=rownames(merged_clust_rev)[new_order_aa]
    
    
    
    melt_df<-melt(merged_clust, id="Pos")
    
    
    
    
    p_heatmap<-ggplot(melt_df)+
      geom_tile(aes(x=factor(Pos, levels=levels),y=factor(variable, levels=rev(levels_aa)),fill=value), size=0.1, color="white")+
      theme_minimal()+
      theme()+
      theme(axis.ticks.y=element_blank(),
            axis.ticks.x.top = element_line(),       
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position = "bottom")+
      labs(y="mutation", x="position", fill="Nucleation score", title = paste0("insertions ", insertions_before_after))+
      scale_fill_gradientn(colours=cols, limits=c(-min,max), na.value = "grey60") 
    p_heatmap
    
    
    
    pdf(paste0(path,"/clustering_ins_", insertions_before_after,"_", reg, ".pdf"), height = 12, width = 12)
    
    grid.newpage()
    print(p_heatmap, vp = viewport(x = 0.4, y = 0.4, width = 0.8, height = 0.8))
    print(dendro_plot, vp = viewport(x = 0.48, y = 0.9, width = 0.95, height = 0.2))
    print(dendro_plot_rev, vp=viewport(x=0.9, y=0.4, width = 0.2, height = 0.8))
    
    dev.off()
    
    
    
    
    
    # i will also keep versions with no clustering in positions
    pos_levels_no_clust<-c(1:42)
    
    if(reg=="NT"){pos_levels_no_clust<-c(1:28)}
    if(reg=="CT"){pos_levels_no_clust<-c(29:42)}
    
    p_heatmap<-ggplot(melt_df)+
      geom_tile(aes(x=factor(Pos, levels=pos_levels_no_clust),y=factor(variable, levels=levels_aa),fill=value), size=0.1, color="white")+
      labs(y="mutation", x="position", fill="Nucleation score", title = paste0("insertions ", insertions_before_after))+
      theme_minimal()+
      theme()+
      theme(axis.ticks.y=element_blank(),
            axis.ticks.x.top = element_line(),       
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position = "bottom")+
      scale_fill_gradientn(colours=cols, limits=c(-min,max), na.value = "grey60") 
    p_heatmap
    
    
    
    #pdf(paste0(path,"/clustering_onlyAA_ins_", insertions_before_after,"_", reg,".pdf"), height = 12, width = 12)
    
    #grid.newpage()
    #print(p_heatmap, vp = viewport(x = 0.4, y = 0.4, width = 0.8, height = 0.8))
    #print(dendro_plot_rev, vp=viewport(x=0.9, y=0.4, width = 0.2, height = 0.8))
    
    #dev.off()
    
    
    
    
    # and a version with no clustering in AA
    
    aa_levels_no_clust<-c(paste0("s", all_aa), paste0("i", all_aa), "del")
    
    
    p_heatmap<-ggplot(melt_df)+
      geom_tile(aes(x=factor(Pos, levels=levels),y=factor(variable, levels=rev(aa_levels_no_clust)),fill=value), size=0.1, color="white")+
      labs(y="mutation", x="position", fill="Nucleation score", title = paste0("insertions ", insertions_before_after))+
      theme_minimal()+
      theme()+
      theme(axis.ticks.y=element_blank(),
            axis.ticks.x.top = element_line(),       
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position = "bottom")+
      scale_fill_gradientn(colours=cols, limits=c(-min,max), na.value = "grey60") 
    p_heatmap
    
    
    #pdf(paste0(path,"/clustering_onlyPos_ins_", insertions_before_after,"_", reg,".pdf"), height = 12, width = 12)
    
    #grid.newpage()
    #print(p_heatmap, vp = viewport(x = 0.4, y = 0.4, width = 0.8, height = 0.8))
    #print(dendro_plot, vp = viewport(x = 0.48, y = 0.9, width = 0.95, height = 0.2))
    
    #dev.off()

  }}
