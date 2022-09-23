## NS distributions


## required packages
require(ggplot2)
require(stringr)
require(ggpubr)
require(dplyr)


dir.create("_NS distributions")
path="_NS distributions"


##  import required data 
load("INDEL.df.RData")



#########################################################################################################################################









# dataset distributions

levels_dataset=c("single", "insertion","single_deletion", "deletion", "truncation")
labels_dataset=c("1AA substitutions", "1AA insertions","1AA deletions", "Multi-AA deletions", "Truncations")


n_text<-as.data.frame(table(INDEL.df$dataset))
colnames(n_text)<-c("dataset", "n")

p_dist<-ggplot(INDEL.df[INDEL.df$dataset!="WT",], aes(x=nscore_c))+
  geom_histogram(bins=50, fill="grey70")+
  geom_vline(xintercept = 0,linetype="dashed", size=0.2)+
  facet_wrap(~factor(dataset, levels=levels_dataset, labels=labels_dataset),ncol=5,scales = "free_y")+
  theme_bw()+
  geom_text(data=n_text[n_text$dataset!="WT",],aes(label=paste0("n=",n), x=Inf, y=Inf),hjust=3,vjust=1.5, size=4, colour="black")+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=0.1),
        axis.ticks = element_line(size=0.1))+
  labs(x="Nucleation score", y="Counts")
p_dist

#ggsave(p_dist, file="p_NS_distributions.pdf", width = 9, height = 2, path=path)




# distributions for APRs and other regions


for(APR1 in c("17_21", "16_23")){
  dist_df<-INDEL.df
  dist_df$part<-""
  dist_df$part_fdr<-""
  
  if(APR1=="17_21"){
    reg1<-c(1:11)
    reg2<-c(12:16)
    reg3<-c(17:21)
    reg4<-c(22:28)
    reg5<-c(29:42)
  }
  
  
  if(APR1=="16_23"){
    reg1<-c(1:11)
    reg2<-c(12:15)
    reg3<-c(16:23)
    reg4<-c(24:28)
    reg5<-c(29:42)
  }
  
  
  for(i in 1:nrow(dist_df)){
    
    # singles
    if(dist_df[i,]$dataset=="single"){
      
      pos<-as.numeric(unlist(strsplit(dist_df[i,]$ID, "-"))[2])
      
      if(pos %in% reg1){dist_df[i,]$part<-paste0("reg_", reg1[1], "_", reg1[length(reg1)]) }
      if(pos %in% reg2){dist_df[i,]$part<-paste0("reg_", reg2[1], "_", reg2[length(reg2)])}
      if(pos %in% reg3){dist_df[i,]$part<-paste0("reg_", reg3[1], "_", reg3[length(reg3)])}
      if(pos %in% reg4){dist_df[i,]$part<-paste0("reg_", reg4[1], "_", reg4[length(reg4)])}
      if(pos %in% reg5){dist_df[i,]$part<-paste0("reg_", reg5[1], "_", reg5[length(reg5)])}
      
      
      if(pos<=28){dist_df[i,]$part_fdr<-"NT"}
      if(pos>=29){dist_df[i,]$part_fdr<-"CT"}
      
    }
    
    # single deletions
    if(dist_df[i,]$dataset=="single_deletion"){
      
      pos0<-unlist(strsplit(dist_df[i,]$ID, "_"))[3]
      pos<-as.numeric(unlist(strsplit(pos0, "-"))[1])
      
      
      if(pos %in% reg1){dist_df[i,]$part<-paste0("reg_", reg1[1], "_", reg1[length(reg1)]) }
      if(pos %in% reg2){dist_df[i,]$part<-paste0("reg_", reg2[1], "_", reg2[length(reg2)])}
      if(pos %in% reg3){dist_df[i,]$part<-paste0("reg_", reg3[1], "_", reg3[length(reg3)])}
      if(pos %in% reg4){dist_df[i,]$part<-paste0("reg_", reg4[1], "_", reg4[length(reg4)])}
      if(pos %in% reg5){dist_df[i,]$part<-paste0("reg_", reg5[1], "_", reg5[length(reg5)])}
      
      
      if(pos<=28){dist_df[i,]$part_fdr<-"NT"}
      if(pos>=29){dist_df[i,]$part_fdr<-"CT"}
    }
    
    
    # single insertions
    if(dist_df[i,]$dataset=="insertion"){
      
      pos<-as.numeric(unlist(strsplit(dist_df[i,]$ID, "_"))[3])
      
      if(pos %in% reg1){dist_df[i,]$part<-paste0("reg_", reg1[1], "_", reg1[length(reg1)]) }
      if(pos %in% reg2){dist_df[i,]$part<-paste0("reg_", reg2[1], "_", reg2[length(reg2)])}
      if(pos %in% reg3){dist_df[i,]$part<-paste0("reg_", reg3[1], "_", reg3[length(reg3)])}
      if(pos %in% reg4){dist_df[i,]$part<-paste0("reg_", reg4[1], "_", reg4[length(reg4)])}
      if(pos %in% reg5){dist_df[i,]$part<-paste0("reg_", reg5[1], "_", reg5[length(reg5)])}
      
      
      if(pos<=28){dist_df[i,]$part_fdr<-"NT"}
      if(pos>28){dist_df[i,]$part_fdr<-"CT"}
    }
    
    # deletions
    if(dist_df[i,]$dataset=="deletion"){
      
      pos0<-unlist(strsplit(dist_df[i,]$ID, ";"))[1]
      pos<-unlist(strsplit(pos0, "_"))[3]
      
      start<-as.numeric(unlist(strsplit(pos, "-"))[1])
      end<-as.numeric(unlist(strsplit(pos, "-"))[2])
      
      
      region=seq(start, end, 1)
      
      if(start<(reg4[1]) & end<(reg5[1]) & end>(reg3[1]-1)){dist_df[i,]$part<-"APR1"}
      if(start>(reg4[1]-1) & end>(reg5[1]-1)){dist_df[i,]$part<-"APR2"}
      if(start<(reg4[1]) & end>(reg5[1]-1)){dist_df[i,]$part<-"APR1_APR2"}
      
      if( sum(region %in% c(reg3, reg5))==0){dist_df[i,]$part<-"no APRs"}
      
      
      
      if(end<(reg5[1])){dist_df[i,]$part_fdr<-"NT"}
      if(start>(reg5[1]-1)){dist_df[i,]$part_fdr<-"CT"}
      if(start<(reg5[1]) & end>(reg5[1]-1)){dist_df[i,]$part_fdr<-"both"}
      
    }
    
    
    # truncations
    if( dist_df[i,]$dataset=="truncation"){
      
      pos00<-as.vector(unlist(strsplit(dist_df[i,]$ID, ";")))
      pos0<-pos00[!is.na(str_extract(pos00, "kmer"))]
      pos<-unlist(strsplit(pos0, "_"))[2]
      
      start<-as.numeric(unlist(strsplit(pos, "-"))[1])
      end<-as.numeric(unlist(strsplit(pos, "-"))[2])
      
      
      del_region=c(1:42)[-seq(start, end, 1)]
      
      if( sum(del_region %in% reg3)!=0){dist_df[i,]$part<-"APR1"}
      if( sum(del_region %in% reg5)!=0){dist_df[i,]$part<-"APR2"}
      if( sum(del_region %in% reg3)!=0 & sum(del_region %in% reg5)!=0  ){dist_df[i,]$part<-"APR1_APR2"}
      if( sum(del_region %in% c(reg3, reg5))==0){dist_df[i,]$part<-"no APRs"}
      
      
      dist_df[i,]$part_fdr<-"both"
      if(start<(reg5[1]) & end==42 ){dist_df[i,]$part_fdr<-"NT"}
      if(end>(reg5[1]-1) & start==1){dist_df[i,]$part_fdr<-"CT"}
      
    }
    
    
  }
  
  
  
  levels_part=c("APR1", "APR2", "APR1_APR2", "no APRs",
                paste0("reg_", reg1[1], "_", reg1[length(reg1)]),
                paste0("reg_", reg2[1], "_", reg2[length(reg2)]),
                paste0("reg_", reg3[1], "_", reg3[length(reg3)]),
                paste0("reg_", reg4[1], "_", reg4[length(reg4)]),
                paste0("reg_", reg5[1], "_", reg5[length(reg5)]) )
  
  
  p_dist_part<-ggplot(dist_df[dist_df$dataset!="WT" ,], aes(y=nscore_c, x=factor(part, levels=levels_part)))+
    geom_hline(yintercept = 0, size=0.1)+
    geom_violin(show.legend = F, size=0.2)+
    geom_jitter(size=0.1, width=0.1, color="grey80")+
    geom_boxplot(width=0.1,outlier.shape = NA, fill=NA, size=0.2)+
    facet_wrap(~factor(dataset, levels=levels_dataset,labels = labels_dataset),ncol=5, scales="free_x" )+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(size=0.2),
          axis.ticks = element_line(size=0.1),
          axis.text.x = element_text(angle = 45, hjust=1))+
    labs(y="Nucleation score", x="")
  p_dist_part
  #ggsave(p_dist_part, file=paste0("p_NS_distributions_NT_CT_",APR1, ".pdf"), width = 10, height = 3, path=path)
  
  
  # both distribution plots together
  p_both<-ggarrange(p_dist, p_dist_part, ncol=1, align = "v")
  p_both
  ggsave(p_both, file=paste0("p_NS_distributions_both_", APR1, ".pdf"), width = 8, height = 4, path=path)
  
  
  
}





# distributions for NT and CT separately. some mutations will affect both 
dist_df<-INDEL.df
dist_df$part<-""

for(i in 1:nrow(dist_df)){
  
  # singles
  if(dist_df[i,]$dataset=="single"){
    
    pos<-as.numeric(unlist(strsplit(dist_df[i,]$ID, "-"))[2])
    
    if(pos<=28){dist_df[i,]$part<-"NT"}
    if(pos>=29){dist_df[i,]$part<-"CT"}
  }
  
  # single deletions
  if(dist_df[i,]$dataset=="single_deletion"){
    
    pos0<-unlist(strsplit(dist_df[i,]$ID, "_"))[3]
    pos<-as.numeric(unlist(strsplit(pos0, "-"))[1])
    
    if(pos<=28){dist_df[i,]$part<-"NT"}
    if(pos>=29){dist_df[i,]$part<-"CT"}
  }
  
  
  # single insertions
  if(dist_df[i,]$dataset=="insertion"){
    
    pos<-as.numeric(unlist(strsplit(dist_df[i,]$ID, "_"))[3])
    
    if(pos<=28){dist_df[i,]$part<-"NT"}
    if(pos>28){dist_df[i,]$part<-"CT"}
  }
  
  # deletions
  if(dist_df[i,]$dataset=="deletion"){
    
    pos0<-unlist(strsplit(dist_df[i,]$ID, ";"))[1]
    pos<-unlist(strsplit(pos0, "_"))[3]
    
    start<-as.numeric(unlist(strsplit(pos, "-"))[1])
    end<-as.numeric(unlist(strsplit(pos, "-"))[2])
    
    if(end<29){dist_df[i,]$part<-"NT"}
    if(start>28){dist_df[i,]$part<-"CT"}
    if(start<29 & end>28){dist_df[i,]$part<-"both"}
  }
  
  
  # truncations
  if( dist_df[i,]$dataset=="truncation"){
    
    pos00<-as.vector(unlist(strsplit(dist_df[i,]$ID, ";")))
    pos0<-pos00[!is.na(str_extract(pos00, "kmer"))]
    pos<-unlist(strsplit(pos0, "_"))[2]
    
    start<-as.numeric(unlist(strsplit(pos, "-"))[1])
    end<-as.numeric(unlist(strsplit(pos, "-"))[2])
    
    dist_df[i,]$part<-"both"
    if(start<29 & end==42 ){dist_df[i,]$part<-"NT"}
    if(end>28 & start==1){dist_df[i,]$part<-"CT"}
  }
  
  #silent
  if( dist_df[i,]$dataset=="silent"){dist_df[i,]$part<-"all"}
  
}


# stacked barplot FDR10 NS categories for each dataset

fdr_categories<-dist_df[,c("aa_seq","ID", "nscore_c","sigma","zscore", "p.adjust", "category_fdr", "dataset", "part")]

fdr_categories$category<-"WT-like"
fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c<0),]$category<- "NS- 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c<0),]$category<- "NS- 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c<0),]$category<- "NS- 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c<0),]$category<- "NS- 1%"
fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c>0),]$category<- "NS+ 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c>0),]$category<- "NS+ 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c>0),]$category<- "NS+ 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c>0),]$category<- "NS+ 1%"



# Merge overall categories and splitted by region
categories <- as.data.frame(fdr_categories %>% group_by(dataset,category) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))
categories$part<-"all"

categories_regions <- as.data.frame(fdr_categories %>% group_by(part,dataset,category) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))

categories_allmut <- as.data.frame(fdr_categories %>% group_by(part,category) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))
categories_allmut$dataset<-"all"

categories_all_all <- as.data.frame(fdr_categories %>% group_by(category) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))
categories_all_all$part<-"all"
categories_all_all$dataset<-"all"


all_categories_regions<-rbind(categories, categories_regions, categories_allmut, categories_all_all)



# plot
colors_fdr <- c(colorRampPalette(c( "brown3", "grey95", "darkblue"))(9))

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")

levels_part=c("all", "NT", "CT", "both")
labels_part=c("All", "N-terminus", "C-terminus", "N and C-terminus")

levels_plot=c(levels_dataset, "all")
labels_plot=c(labels_dataset, "All")

p_categories_region<-ggplot(all_categories_regions[all_categories_regions$dataset %in% levels_plot & all_categories_regions$part %in% levels_part,], 
                            aes(fill=factor(category, levels=rev(levels)), 
                                y=factor(dataset, levels=rev(levels_plot), labels = rev(labels_plot)), 
                                x=freq)) + 
  facet_wrap(~factor(part, levels=levels_part, labels=labels_part), ncol=4)+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(color='black', size=0.25),
        strip.background = element_rect(fill=NA, color=NA))+
  scale_fill_manual("FDR category", values=rev(colors_fdr))+
  labs( x= "Frequency",y="")
p_categories_region

ggsave(p_categories_region, file="p_fdr_datasets_regions.pdf", width = 10, height = 2, path=path)






# count number of variants
print(as.data.frame(dist_df %>% group_by(category_fdr) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100)))
print(as.data.frame(dist_df %>% group_by(category_fdr, dataset) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100)))
print(as.data.frame(dist_df %>% group_by(category_fdr, part) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100)))
print(as.data.frame(dist_df %>% group_by(dataset, category_fdr) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100)))
print(as.data.frame(dist_df %>% group_by(dataset,part, category_fdr) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100)))
print(as.data.frame(dist_df %>% group_by(dataset, category_fdr, part) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100)))






# barplot FDR categories
plot_df<-as.data.frame(dist_df[dist_df$dataset!="WT",] %>% group_by(category_fdr, part) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100))

levels_fdr = c("NS_dec", "WT-like", "NS_inc")
labels_fdr=c( "NS-", "WT-like", "NS+")

colors=c("grey55","grey68", "grey80")
levels_part=c("NT", "CT", "both")
labels_part=c("N-terminus", "C-terminus", "N and C-terminus")


p_barplot<-ggplot(plot_df, aes(x=factor(category_fdr, levels=levels_fdr, labels=labels_fdr), 
                         fill=factor(part, levels=levels_part, labels=labels_part), 
                         y=freq)) + 
  geom_bar(position="stack", stat="identity", width = 0.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.y = element_line(color='black', size=0.25),
        axis.ticks.x = element_blank())+
  scale_fill_manual("", values=colors)+
  geom_text(aes(label = n),size=7,  position = position_stack(vjust = .5)) +
  labs(x="", y="Frequency")
p_barplot

ggsave(p_barplot,file="p_barplot_fdr_freq_regions.pdf", width =6, height = 4, path=path )




# pie charts NS+ mutations

inc_all_df<-as.data.frame(dist_df[dist_df$category_fdr=="NS_inc",] %>% group_by(dataset) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100))
inc_all_df$part<-"all"
inc_regions_df<-as.data.frame(dist_df[dist_df$category_fdr=="NS_inc",] %>% group_by(part, dataset) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100))

inc_plot_df<-rbind(inc_all_df, inc_regions_df)


p_pie_all<-ggplot(inc_plot_df, aes(x="", fill=factor(dataset, levels=levels_dataset, labels = labels_dataset), y=freq)) + 
  geom_bar(stat="identity")+
  facet_wrap(~factor(part, levels=c("all", levels_part), labels=c("All", labels_part)), ncol = 4)+
  coord_polar("y", start=0) +
  theme_void()+
  
  scale_fill_brewer("Dataset", palette = "YlGn", direction = -1)+
  geom_text(aes(label = n),  position = position_stack(vjust = .5)) +
  labs(x="", y="Frequency")
p_pie_all

ggsave(p_pie_all, file="p_pie_inc.pdf", width = 8, height = 2, path = path)



# pie charts NS- mutations

dec_all_df<-as.data.frame(dist_df[dist_df$category_fdr=="NS_dec",] %>% group_by(dataset) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100))
dec_all_df$part<-"all"
dec_regions_df<-as.data.frame(dist_df[dist_df$category_fdr=="NS_dec",] %>% group_by(part, dataset) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100))

dec_plot_df<-rbind(dec_all_df, dec_regions_df)


p_pie_all<-ggplot(dec_plot_df, aes(x="", fill=factor(dataset, levels=levels_dataset, labels = labels_dataset), y=freq)) + 
  geom_bar(stat="identity")+
  facet_wrap(~factor(part, levels=c("all", levels_part), labels=c("All", labels_part)), ncol = 4)+
  coord_polar("y", start=0) +
  theme_void()+
  scale_fill_brewer("Dataset", palette = "YlGn", direction = -1)+
  geom_text(aes(label = n),  position = position_stack(vjust = .5)) +
  labs(x="", y="Frequency")
p_pie_all

ggsave(p_pie_all, file="p_pie_dec.pdf", width = 8, height = 2, path = path)













