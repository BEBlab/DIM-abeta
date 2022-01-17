## Single AA mutations : substitutions, insertions and deletions


## required packages
require(ggplot2)
require(ggpubr)
require(stringr)
require(dplyr)


dir.create("_Single aa mutations")
path="_Single aa mutations"


##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")
load("required data/non_nucleating_variants.RData")



#########################################################################################################################################


AB_wt<-"DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
ABseq<-unlist(strsplit(AB_wt, ""))
ABseq_pos<-paste0(ABseq, "\n", c(1:42))

all_aa=c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

disease_mutations<-INDEL.df[INDEL.df$fAD=="fAD_d" & INDEL.df$dataset=="single",]$ID



## NS color gradient 
min<-abs(min(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))
max<-abs(max(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))

cols <- c(colorRampPalette(c( "brown3", "grey95"))((min/(min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(min+max)*100)-0.5))










### single AA substitutions

#all possible singles
single_subst=c()
for (i in c(1:42)){
  wt_AA<-unlist(strsplit(AB_wt, ""))[i]
  
  for (j in all_aa){
    new_mutation<-paste0(wt_AA,"-", i, "-", j)
    single_subst<-c(single_subst, new_mutation)
    
  }}




singles<-singles.df[,c("aa_seq", "ID", "nscore_c", "sigma")]

missing_singles<-data.frame("aa_seq"=NA,
                            "ID"=single_subst[!single_subst %in% singles$ID],
                            "nscore_c"=NA,
                            "sigma"=0)

heatmap_singles<-rbind(singles, missing_singles)


heatmap_singles$WT_AA<-""
heatmap_singles$Pos<-""
heatmap_singles$Mut<-""

heatmap_singles$ID<-as.character(heatmap_singles$ID)
for(i in 1:nrow(heatmap_singles)){
  heatmap_singles[i,]$WT_AA<-unlist(strsplit(heatmap_singles[i,]$ID, "-"))[1]
  heatmap_singles[i,]$Pos<-unlist(strsplit(heatmap_singles[i,]$ID, "-"))[2]
  heatmap_singles[i,]$Mut<-unlist(strsplit(heatmap_singles[i,]$ID, "-"))[3]
}
heatmap_singles$Pos<-as.numeric(heatmap_singles$Pos)




## heatmap

#add info fAD
heatmap_singles$box<-"VUS"
heatmap_singles[heatmap_singles$ID %in% disease_mutations,]$box<-"Dominant"


#add info syn and info on non-nucleating variants
heatmap_singles$label<-""
heatmap_singles[heatmap_singles$WT_AA == heatmap_singles$Mut,]$label<-"*"
heatmap_singles[heatmap_singles$WT_AA == heatmap_singles$Mut,]$nscore_c<-0

heatmap_singles[heatmap_singles$ID %in% non_nuc_df$ID,]$label<-"-"




p_heatmap<-ggplot(heatmap_singles)+
  geom_tile(aes(x=Pos,y=factor(Mut, levels=rev(all_aa)),fill=nscore_c), size=0.1, color="white")+
  geom_tile(data=heatmap_singles[heatmap_singles$box=="Dominant",], 
            aes(x=Pos,y=factor(Mut, levels=rev(all_aa))),color="grey20", fill=NA, size=1)+
  theme_minimal()+
  theme()+
  scale_x_continuous(breaks=seq(1:42), labels = ABseq_pos, expand = c(0,0))+
  labs(x="AB(1-42) WT amino acid and position", y="Mutant amino acid", fill="Nucleation\nscore")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x.top = element_line(),       
        plot.title =element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15), 
        axis.text = element_text( size=12),
        axis.title = element_text(size = 20))+
  geom_text(aes(Pos,factor(Mut, levels=rev(all_aa)),label=label), size=7)+
  scale_fill_gradientn(colours=cols, limits=c(-min,max), na.value = "grey60") 
p_heatmap

ggsave(p_heatmap, file="p_heatmap_single_subst.pdf",width=14, height=8, path=path)







#### single AA insertions
# use insertion_reps dataframe to visualise duplicated coding sequences

insertions<-insertions_reps[,c("aa_seq", "ins_pos", "ins_aa", "nscore_c", "sigma")]

my_ins<-c(paste0(insertions$ins_pos,"-", insertions$ins_aa))

#add missing insertions to heatmap
missing_ins=c()
for (i in c(1:41)){
  for (j in all_aa){
    this_ins<-paste0(i,"-",j)
    if(! this_ins %in% my_ins){  missing_ins=c(missing_ins,i,j)}
  }}


my_df<-as.data.frame(matrix(missing_ins, ncol=2, byrow = T))

missing_insertions<-data.frame("aa_seq"=NA,
                                "ins_pos"=my_df$V1,
                               "ins_aa"=my_df$V2,
                               "nscore_c"=NA,
                               "sigma"=NA)

heatmap_insertions<-rbind(insertions, missing_insertions)

heatmap_insertions$ID<-paste0("Ins_k1_", heatmap_insertions$ins_pos, "_", heatmap_insertions$ins_aa)

#add info on non-nucleating variants
heatmap_insertions$label<-""
heatmap_insertions[heatmap_insertions$ID %in% non_nuc_df$ID,]$label<-"-"


p_heatmap<-ggplot(heatmap_insertions)+
  geom_tile(aes(factor(ins_pos, levels=c(1:41), labels=c(2:42)),
                factor(ins_aa, levels=rev(all_aa)),fill=nscore_c), color="white", size=0.1)+
  theme_minimal()+
  theme()+
  labs(x="Position of inserted amino acid", y="Inserted amino acid", fill="Nucleation\nscore")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x.top = element_line(),       
        plot.title =element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15), 
        axis.title = element_text(size = 16))+
  geom_text(aes(factor(ins_pos, levels=c(1:41), labels=c(2:42)),
                factor(ins_aa, levels=rev(all_aa)),label=label), size=7)+
  scale_fill_gradientn(colours=cols, limits=c(-min,max), na.value = "grey60") 
p_heatmap

ggsave(p_heatmap, file="p_heatmap_insertions.pdf",width=10, height=6, path=path)






# violin plots

subs_median_df<-as.data.frame(singles.df %>% group_by(Pos) %>% dplyr::summarise(median=median(nscore_c)))
subs_plot_df<-left_join(singles.df, subs_median_df)

p_subs_violin<-ggplot(subs_plot_df, aes(x=factor(Pos, levels=c(1:42), labels=ABseq_pos), 
                            y=nscore_c, group=factor(Pos, levels=c(1:42), labels=ABseq_pos)))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_violin(scale = "width", aes(fill=median), size=0.2)+
  geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
  theme_bw()+
  labs(x="AB(1-42) WT amino acid and position", y="Nucleation score", fill="Median NS")+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line())+
  scale_fill_gradientn(colours=cols, limits=c(-min,max), na.value = "grey60") 
p_subs_violin

ggsave(p_subs_violin, file="p_subs_violin.pdf", width = 12, height = 2.5, path=path)






ins_median_df<-as.data.frame(insertions_reps %>% group_by(ins_pos) %>% dplyr::summarise(median=median(nscore_c)))
ins_plot_df<-left_join(insertions_reps, ins_median_df)

p_ins_violin<-ggplot(ins_plot_df, aes(x=factor(ins_pos, levels=c(1:41), labels=c(2:42)), 
                                           y=nscore_c, group=factor(ins_pos, levels=c(2:42), labels=c(2:42))))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_violin(scale = "width", aes(fill=median), size=0.2)+
  geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
  theme_bw()+
  labs(x="Position of inserted amino acid", y="Nucleation score", fill="Median NS")+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line())+
  scale_fill_gradientn(colours=cols, limits=c(-min,max), na.value = "grey60") 
p_ins_violin

ggsave(p_ins_violin, file="p_ins_violin.pdf", width = 12, height = 2.5, path=path)








#### single AA deletions
# use single_deletions_reps dataframe to visualise duplicated coding sequences

p_single_deletion<-ggplot(single_deletions_reps, aes(x=factor(del_pos, levels = c(1:42)), y=nscore_c ))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_errorbar(aes(ymin=nscore_c-1.96*sigma, ymax=nscore_c+1.96*sigma), width=0, size=0.1)+
  
  geom_point(data=single_deletions_reps[single_deletions_reps$category_fdr %in% c("NS_inc", "NS_dec"),],
             aes(fill=nscore_c),size=4, shape=21, stroke=1.2)+
  geom_point(data=single_deletions_reps[!single_deletions_reps$category_fdr %in% c("NS_inc", "NS_dec"),],
             aes(fill=nscore_c),size=4, shape=21, stroke=0.2)+
  
  labs(y="Nucleation score",x="Deleted amino acid and position", fill="Nucleation score")+
  theme_bw()+
  theme(axis.text.x=element_text(size=8),
        panel.border = element_blank(),
        axis.line.y = element_line(color="black"),
        panel.grid.major = element_line(size=0.2), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom")+
  scale_x_discrete(breaks=c(1:42), labels=ABseq_pos)+
  scale_fill_gradientn(colors=cols, limits=c(-min,max))

p_single_deletion

ggsave(p_single_deletion, file="p_single_deletions.pdf", width = 8, height = 3, path=path)











# stacked barplots for each AA


levels_fdr = c("NS_dec", "WT-like", "NS_inc")
labels_fdr=c("NS-", "WT-like", "NS+")
colors<-c( "brown3", "grey90", "darkblue")

levels_part<-c("all", "NT", "CT")
labels_part=c("All", "N-terminus", "C-terminus")


# substitutions

aa_counts_subs<-as.data.frame(singles.df %>% group_by(Mut, category_fdr) %>% dplyr::summarise(n=n()) %>%  dplyr::mutate(freq=n / sum(n)))
aa_counts_subs$part<-"all"
aa_counts_subs_nt<-as.data.frame(singles.df[singles.df$Pos %in% c(1:28),] %>% group_by(Mut, category_fdr) %>% dplyr::summarise(n=n()) %>%  dplyr::mutate(freq=n / sum(n)))
aa_counts_subs_nt$part<-"NT"
aa_counts_subs_ct<-as.data.frame(singles.df[singles.df$Pos %in% c(29:42),] %>% group_by(Mut, category_fdr) %>% dplyr::summarise(n=n()) %>%  dplyr::mutate(freq=n / sum(n)))
aa_counts_subs_ct$part<-"CT"
aa_counts_plot<-rbind(aa_counts_subs, aa_counts_subs_nt, aa_counts_subs_ct)



p_counts_aa_subs<-ggplot(aa_counts_plot, aes(fill=factor(category_fdr, levels=levels_fdr, labels=labels_fdr), 
                                y=freq, x=factor(Mut, levels=all_aa))) + 
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  facet_wrap(~factor(part, levels=levels_part, labels=labels_part), ncol = 1)+
  
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(color='black', size=0.25),
        strip.background = element_rect(fill=NA, color=NA))+
  scale_fill_manual("FDR category", values=colors)+
  labs(y="Frequency",x="Mutant amino acid")
p_counts_aa_subs

ggsave(p_counts_aa_subs, file="p_counts_aa_subs.pdf",  width = 8, height = 4, path=path)



# substituted (WT) AA

aa_counts_subs<-as.data.frame(singles.df %>% group_by(WT_AA, category_fdr) %>% dplyr::summarise(n=n()) %>%  dplyr::mutate(freq=n / sum(n)))
aa_counts_subs$part<-"all"
aa_counts_subs_nt<-as.data.frame(singles.df[singles.df$Pos %in% c(1:28),] %>% group_by(WT_AA, category_fdr) %>% dplyr::summarise(n=n()) %>%  dplyr::mutate(freq=n / sum(n)))
aa_counts_subs_nt$part<-"NT"
aa_counts_subs_ct<-as.data.frame(singles.df[singles.df$Pos %in% c(29:42),] %>% group_by(WT_AA, category_fdr) %>% dplyr::summarise(n=n()) %>%  dplyr::mutate(freq=n / sum(n)))
aa_counts_subs_ct$part<-"CT"
aa_counts_plot<-rbind(aa_counts_subs, aa_counts_subs_nt, aa_counts_subs_ct)



p_counts_aa_subs<-ggplot(aa_counts_plot, aes(fill=factor(category_fdr, levels=levels_fdr, labels=labels_fdr), 
                                             y=freq, x=factor(WT_AA, levels=all_aa))) + 
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  facet_wrap(~factor(part, levels=levels_part, labels=labels_part), ncol = 1)+
  
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(color='black', size=0.25),
        strip.background = element_rect(fill=NA, color=NA))+
  scale_fill_manual("FDR category", values=colors)+
  labs(y="Frequency",x="Substituted amino acid")
p_counts_aa_subs

ggsave(p_counts_aa_subs, file="p_counts_aa_subs_removed_aa.pdf",  width = 8, height = 4, path=path)





# insertions

aa_counts_ins<-as.data.frame(insertions_reps %>% group_by(ins_aa, category_fdr) %>% dplyr::summarise(n=n()) %>%  dplyr::mutate(freq=n / sum(n)))
aa_counts_ins$part<-"all"
aa_counts_ins_nt<-as.data.frame(insertions_reps[insertions_reps$ins_pos %in% c(1:28),] %>% group_by(ins_aa, category_fdr) %>% dplyr::summarise(n=n()) %>%  dplyr::mutate(freq=n / sum(n)))
aa_counts_ins_nt$part<-"NT"
aa_counts_ins_ct<-as.data.frame(insertions_reps[insertions_reps$ins_pos %in% c(29:42),] %>% group_by(ins_aa, category_fdr) %>% dplyr::summarise(n=n()) %>%  dplyr::mutate(freq=n / sum(n)))
aa_counts_ins_ct$part<-"CT"
aa_counts_plot<-rbind(aa_counts_ins, aa_counts_ins_nt, aa_counts_ins_ct)



p_counts_aa_ins<-ggplot(aa_counts_plot, aes(fill=factor(category_fdr, levels=levels_fdr, labels=labels_fdr), 
                                             y=freq, x=factor(ins_aa, levels=all_aa))) + 
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  facet_wrap(~factor(part, levels=levels_part, labels=labels_part), ncol = 1)+
  
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line.x = element_line(color='black', size=0.25),
        strip.background = element_rect(fill=NA, color=NA))+
  scale_fill_manual("FDR category", values=colors)+
  labs(y="Frequency",x="Inserted amino acid")
p_counts_aa_ins

ggsave(p_counts_aa_ins, file="p_counts_aa_ins.pdf",  width = 8, height = 4, path=path)











### supplementary dotplots




levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")

#colors <- c(colorRampPalette(c( "brown3", "grey95", "darkblue"))(9))
colors<-c( "NS- 1%"="#CD3333", 
           "NS- 5%"="#D66262", 
           "NS- 10%"="#DF9292", 
           "NS- 25%"="#E8C2C2",
           "WT-like"="#F2F2F2", 
           "NS+ 25%"="#B5B5D8",
           "NS+ 10%"="#7979BE",
           "NS+ 5%"="#3C3CA4",
           "NS+ 1%"="#00008B")



fdr_categories<-singles.df

fdr_categories$category<-"WT-like"
fdr_categories[(!is.na(fdr_categories$nscore_c) & fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c<0),]$category<- "NS- 25%"
fdr_categories[(!is.na(fdr_categories$nscore_c) & fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c<0),]$category<- "NS- 10%"
fdr_categories[(!is.na(fdr_categories$nscore_c) & fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c<0),]$category<- "NS- 5%"
fdr_categories[(!is.na(fdr_categories$nscore_c) & fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c<0),]$category<- "NS- 1%"
fdr_categories[(!is.na(fdr_categories$nscore_c) & fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c>0),]$category<- "NS+ 25%"
fdr_categories[(!is.na(fdr_categories$nscore_c) & fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c>0),]$category<- "NS+ 10%"
fdr_categories[(!is.na(fdr_categories$nscore_c) & fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c>0),]$category<- "NS+ 5%"
fdr_categories[(!is.na(fdr_categories$nscore_c) & fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c>0),]$category<- "NS+ 1%"





### master boxplot for each AA

plot_list<-list()
for (AA in unique(fdr_categories$Mut)){
  
  p_AA<-ggplot(fdr_categories, aes(x = factor(Pos, levels=c(1:42), labels=ABseq_pos), y = nscore_c)) +
    
    geom_hline(aes(yintercept = 0), size=0.15, colour="black")+
    geom_boxplot( size=0.2, outlier.shape = NA, color="grey90", fill="white") +
    
    geom_jitter(data=fdr_categories[fdr_categories$Mut!=AA,],  color="grey90", fill="grey90")+
    geom_jitter(data=fdr_categories[fdr_categories$Mut==AA,], aes(color=factor(category, levels=levels)))+
    
    scale_color_manual("FDR",values=colors)+
    
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title=element_text(size=25),
          axis.line.y = element_line(color="black", size=0.15))+
    #labs(title=paste0(AA), y="Nucleation score",x="AB(1-42) WT amino acid and position")
  labs(title=paste0(AA), y="",x="")
  
  p_AA
  
  plot_list[[AA]]<-p_AA
  
  #ggsave(p_AA, path=path, file=paste0(AA, "_AA_boxplot.pdf"), width = 12, height = 4)
  
  
}

aa_plot<-unique(fdr_categories$Mut)
p_all<-ggarrange(plotlist = plot_list, ncol=2, nrow = 10, common.legend = T)


ggsave(p_all, file="p_all_dotplots_subst.pdf", path=path, width = 18, height = 25)





# same for insertions

fdr_ins<-insertions_reps

fdr_ins$category<-"WT-like"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.25 & fdr_ins$nscore_c<0),]$category<- "NS- 25%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.1 & fdr_ins$nscore_c<0),]$category<- "NS- 10%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.05 & fdr_ins$nscore_c<0),]$category<- "NS- 5%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.01 & fdr_ins$nscore_c<0),]$category<- "NS- 1%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.25 & fdr_ins$nscore_c>0),]$category<- "NS+ 25%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.1 & fdr_ins$nscore_c>0),]$category<- "NS+ 10%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.05 & fdr_ins$nscore_c>0),]$category<- "NS+ 5%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.01 & fdr_ins$nscore_c>0),]$category<- "NS+ 1%"




plot_list_ins<-list()
for (AA in unique(fdr_ins$ins_aa)){
  
  p_AA<-ggplot(fdr_ins, aes(x = factor((as.numeric(ins_pos)+1), levels=c(2:42)), y = nscore_c)) +
    
    geom_hline(aes(yintercept = 0), size=0.15, colour="black")+
    geom_boxplot( size=0.2, outlier.shape = NA, color="grey90", fill="white") +
    
    geom_jitter(data=fdr_ins[fdr_ins$ins_aa!=AA,],  color="grey90", fill="grey90")+
    geom_jitter(data=fdr_ins[fdr_ins$ins_aa==AA,], aes(color=factor(category, levels=levels)))+
    
    scale_color_manual("FDR",values=colors)+
    
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title=element_text(size=25),
          axis.line.y = element_line(color="black", size=0.15))+
    #labs(title=paste0(AA), y="Nucleation score",x="Position of inserted AA")
    labs(title=paste0(AA), y="",x="")
  p_AA
  
  plot_list_ins[[AA]]<-p_AA
  
  #ggsave(p_AA, path=path, file=paste0(AA, "_AA_boxplot.pdf"), width = 12, height = 4)
  
  
}


aa_plot<-unique(fdr_ins$ins_aa)
p_all_ins<-ggarrange(plotlist = plot_list_ins, ncol=2, nrow = 10, common.legend = T)


ggsave(p_all_ins, file="p_all_dotplots_ins.pdf", path=path, width = 18, height = 25)









