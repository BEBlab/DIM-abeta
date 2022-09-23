## Comparison of single AA mutations effects : substitutions, insertions and deletions


## required packages
require(stringr)
require(ggplot2)
require(ggpubr)
require(dplyr)
require(ggrepel)



dir.create("_Mutation types comparison")
path="_Mutation types comparison"


##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")



#########################################################################################################################################


AB_wt<-"DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
ABseq<-unlist(strsplit(AB_wt, ""))
ABseq_pos<-paste0(ABseq, "\n", c(1:42))

all_aa<-c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

AA_type<-data.frame("AA"= all_aa,
                    "name_AA"=c("Glycine", "Alanine","Valine","Leucine","Methionine","Isoleucine","Phenylalanine",
                                "Tyrosine","Tryptophan","Lysine","Arginine","Aspartic acid","Glutamic acid","Serine","Threonine",
                                "Cysteine","Asparagine","Glutamine","Histidine","Proline"),
                    "type"=c("glycine",rep("aliphatic",5),rep("aromatic",3),rep("positive",2),rep("negative",2),rep("polar",6),"proline"))

color_AAtype<-c("aliphatic"="darkgrey",
               "aromatic"="#9A703EFF",
               "negative"="#EE0011FF",
               "positive"="#0C5BB0FF", 
               "polar"="#15983DFF", 
               "glycine"="grey30", 
               "proline"="#FEC10BFF")


## NS color gradient 
min<-abs(min(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))
max<-abs(max(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))

cols <- c(colorRampPalette(c( "brown3", "grey95"))((min/(min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(min+max)*100)-0.5))










# merge data substitutions, insertions and deletions

insertions_reps$new_pos<-as.numeric(insertions_reps$ins_pos)+1

singles_insertions_dels<-as.data.frame(singles.df[,c("nscore_c", "sigma", "WT_AA", "Pos", "Mut")])
colnames(singles_insertions_dels)[c(1:2)]<-c("nscore_single", "sigma_single")


singles_insertions_dels$nscore_ins_before<-""
singles_insertions_dels$sigma_ins_before<-""
singles_insertions_dels$nscore_ins_after<-""
singles_insertions_dels$sigma_ins_after<-""

singles_insertions_dels$nscore_del<-""
singles_insertions_dels$sigma_del<-""

for(i in 1:nrow(singles_insertions_dels)){
  
  which_pos<-singles_insertions_dels[i,]$Pos
  which_mut<-singles_insertions_dels[i,]$Mut
  
  tryCatch({
    singles_insertions_dels[i,]$nscore_del<-single_deletions_reps[single_deletions_reps$del_pos==which_pos,]$nscore_c
    singles_insertions_dels[i,]$sigma_del<-single_deletions_reps[single_deletions_reps$del_pos==which_pos,]$sigma
  }, error=function(e){})
  
 
  tryCatch({
    singles_insertions_dels[i,]$nscore_ins_after<-insertions_reps[insertions_reps$ins_pos==which_pos & insertions_reps$ins_aa==which_mut ,]$nscore_c
    singles_insertions_dels[i,]$sigma_ins_after<-insertions_reps[insertions_reps$ins_pos==which_pos &  insertions_reps$ins_aa==which_mut ,]$sigma
  }, error=function(e){})
  
  
  tryCatch({
    singles_insertions_dels[i,]$nscore_ins_before<-insertions_reps[insertions_reps$new_pos==which_pos & insertions_reps$ins_aa==which_mut ,]$nscore_c
    singles_insertions_dels[i,]$sigma_ins_before<-insertions_reps[insertions_reps$new_pos==which_pos &  insertions_reps$ins_aa==which_mut ,]$sigma
  }, error=function(e){})
  
}


for(i in c(6:11)){singles_insertions_dels[[i]]<-as.numeric(singles_insertions_dels[[i]])}






#### substitutions and insertions

singles_insertions_dels$part<-""
singles_insertions_dels[singles_insertions_dels$Pos %in% c(1:28),]$part<-"NT"
singles_insertions_dels[singles_insertions_dels$Pos %in% c(29:42),]$part<-"CT"



corr_before=cor.test(singles_insertions_dels$nscore_single, singles_insertions_dels$nscore_ins_before, use="complete.obs")$estimate
p_value_before=cor.test(singles_insertions_dels$nscore_single, singles_insertions_dels$nscore_ins_before, use="complete.obs")$p.value

corr_before_nt=cor.test(singles_insertions_dels[singles_insertions_dels$part=="NT",]$nscore_single, 
                        singles_insertions_dels[singles_insertions_dels$part=="NT",]$nscore_ins_before, use="complete.obs")$estimate
p_value_before_nt=cor.test(singles_insertions_dels[singles_insertions_dels$part=="NT",]$nscore_single, 
                           singles_insertions_dels[singles_insertions_dels$part=="NT",]$nscore_ins_before, use="complete.obs")$p.value
corr_before_ct=cor.test(singles_insertions_dels[singles_insertions_dels$part=="CT",]$nscore_single, 
                        singles_insertions_dels[singles_insertions_dels$part=="CT",]$nscore_ins_before, use="complete.obs")$estimate
p_value_before_ct=cor.test(singles_insertions_dels[singles_insertions_dels$part=="CT",]$nscore_single, 
                           singles_insertions_dels[singles_insertions_dels$part=="CT",]$nscore_ins_before, use="complete.obs")$p.value



p_before<-ggplot(singles_insertions_dels, aes(x=nscore_single, y=nscore_ins_before))+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(xintercept = 0,linetype="dashed", color="black")+
  geom_point(aes(color=factor(part, levels=c("NT", "CT"))))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  scale_color_manual("", values=c("grey40","grey70"))+
 
  annotate("text", label=paste0("R= ", round(corr_before,2)), x=-4, y=2, size=6)+
  annotate("text", label=paste0("p= ", format(p_value_before, digits = 2, scientific = T)), x=-4, y=1.5, size=6)+
  
  annotate("text", label=paste0("[Nt] R= ", round(corr_before_nt,2)), x=-4, y=1, size=6, color="grey40")+
  annotate("text", label=paste0("[Nt] p= ", format(p_value_before_nt, digits = 2, scientific = T)), x=-4, y=0.5, size=6,color="grey40")+
  
  annotate("text", label=paste0("[Ct] R= ", round(corr_before_ct,2)), x=-4, y=0, size=6,color="grey70")+
  annotate("text", label=paste0("[Ct] p= ", format(p_value_before_ct, digits = 2, scientific = T)), x=-4, y=-0.5, size=6,color="grey70")+
  
  labs(x="Substitutions", y="Insertion before position i", color="position")
p_before
#ggsave(p_before, file="p_subst_ins_before.pdf", width = 5, height = 4, path=path)




corr_after=cor.test(singles_insertions_dels$nscore_single, singles_insertions_dels$nscore_ins_after, use="complete.obs")$estimate
p_value_after=cor.test(singles_insertions_dels$nscore_single, singles_insertions_dels$nscore_ins_after, use="complete.obs")$p.value

corr_after_nt=cor.test(singles_insertions_dels[singles_insertions_dels$part=="NT",]$nscore_single, 
                       singles_insertions_dels[singles_insertions_dels$part=="NT",]$nscore_ins_after, use="complete.obs")$estimate
p_value_after_nt=cor.test(singles_insertions_dels[singles_insertions_dels$part=="NT",]$nscore_single, 
                          singles_insertions_dels[singles_insertions_dels$part=="NT",]$nscore_ins_after, use="complete.obs")$p.value
corr_after_ct=cor.test(singles_insertions_dels[singles_insertions_dels$part=="CT",]$nscore_single, 
                       singles_insertions_dels[singles_insertions_dels$part=="CT",]$nscore_ins_after, use="complete.obs")$estimate
p_value_after_ct=cor.test(singles_insertions_dels[singles_insertions_dels$part=="CT",]$nscore_single, 
                          singles_insertions_dels[singles_insertions_dels$part=="CT",]$nscore_ins_after, use="complete.obs")$p.value


p_after<-ggplot(singles_insertions_dels, aes(x=nscore_single, y=nscore_ins_after))+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(xintercept = 0,linetype="dashed", color="black")+
  geom_point(aes(color=factor(part, levels=c("NT", "CT"))))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  scale_color_manual("", values=c("grey40","grey70"))+
  
  annotate("text", label=paste0("R= ", round(corr_after,2)), x=-4, y=2, size=6)+
  annotate("text", label=paste0("p= ", format(p_value_after, digits = 2, scientific = T)), x=-4, y=1.5, size=6)+
  
  annotate("text", label=paste0("[Nt] R= ", round(corr_after_nt,2)), x=-4, y=1, size=6, color="grey40")+
  annotate("text", label=paste0("[Nt] p= ", format(p_value_after_nt, digits = 2, scientific = T)), x=-4, y=0.5, size=6,color="grey40")+
  
  annotate("text", label=paste0("[Ct] R= ", round(corr_after_ct,2)), x=-4, y=0, size=6,color="grey70")+
  annotate("text", label=paste0("[Ct] p= ", format(p_value_after_ct, digits = 2, scientific = T)), x=-4, y=-0.5, size=6,color="grey70")+
  
  labs(x="Substitutions", y="Insertion after position i", color="Position")
p_after
#ggsave(p_after, file="p_subst_ins_after.pdf", width = 5, height = 4, path=path)


p_both<-ggarrange(p_before, p_after, ncol=1, common.legend = T)
p_both
ggsave(p_both, file="p_subst_ins_before_after.pdf", width = 3, height = 6, path=path)





#### substitutions and insertions
# by aa mut

AA_loop<- all_aa[all_aa %in% unique(singles_insertions_dels$Mut)]


corr_vector=c()
for(i in AA_loop){
  corr<-cor.test(singles_insertions_dels[singles_insertions_dels$Mut==i,]$nscore_single,
                 singles_insertions_dels[singles_insertions_dels$Mut==i,]$nscore_ins_before, use="complete.obs")$estimate
  p_value<-cor.test(singles_insertions_dels[singles_insertions_dels$Mut==i,]$nscore_single,
                    singles_insertions_dels[singles_insertions_dels$Mut==i,]$nscore_ins_before, use="complete.obs")$p.value
  corr_vector=c(corr_vector,i,corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Mut", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}


p_before<-ggplot(singles_insertions_dels, aes(x=nscore_single, y=nscore_ins_before))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0,linetype="dashed", color="darkgrey")+
  geom_point(aes(color=as.numeric(Pos)))+
  theme_bw()+
  facet_wrap(~factor(Mut, levels=AA_loop))+
  theme(panel.grid = element_blank())+
  scale_color_viridis_c()+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.15,vjust=1.5, size=4, colour="black")+
  geom_text(data=corr_text, aes(label=paste0("p= ", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.9,vjust=1.5, size=4, colour="black")+
  labs(x="Substitutions", y="Insertion before position i", color="position")
p_before
ggsave(p_before, file="p_subst_ins_before_by mut.pdf", width = 10, height = 8, path=path)



corr_vector=c()
for(i in AA_loop){
  corr<-cor.test(singles_insertions_dels[singles_insertions_dels$Mut==i,]$nscore_single,
                 singles_insertions_dels[singles_insertions_dels$Mut==i,]$nscore_ins_after, use="complete.obs")$estimate
  p_value<-cor.test(singles_insertions_dels[singles_insertions_dels$Mut==i,]$nscore_single,
                    singles_insertions_dels[singles_insertions_dels$Mut==i,]$nscore_ins_after, use="complete.obs")$p.value
  corr_vector=c(corr_vector,i,corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Mut", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}


p_after<-ggplot(singles_insertions_dels, aes(x=nscore_single, y=nscore_ins_after))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0,linetype="dashed", color="darkgrey")+
  geom_point(aes(color=as.numeric(Pos)))+
  theme_bw()+
  facet_wrap(~factor(Mut, levels=AA_loop))+
  theme(panel.grid = element_blank())+
  scale_color_viridis_c()+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.15,vjust=1.5, size=4, colour="black")+
  geom_text(data=corr_text, aes(label=paste0("p= ", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.9,vjust=1.5, size=4, colour="black")+
  labs(x="Substitutions", y="Insertion after position i", color="position")
p_after
ggsave(p_after, file="p_subst_ins_after_by mut.pdf", width = 10, height = 8, path=path)






#### substitutions and insertions
# by pos

corr_vector=c()
for(i in c(2:42)){
  corr<-cor.test(singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_single,
                 singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_ins_before, use="complete.obs")$estimate
  p_value<-cor.test(singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_single,
                    singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_ins_before, use="complete.obs")$p.value
  corr_vector=c(corr_vector,i,corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(1:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}


p_before<-ggplot(singles_insertions_dels, aes(x=nscore_single, y=nscore_ins_before))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0,linetype="dashed", color="darkgrey")+
  geom_point()+
  theme_bw()+
  facet_wrap(~as.numeric(Pos))+
  theme(panel.grid = element_blank())+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.15,vjust=1.5, size=4, colour="black")+
  geom_text(data=corr_text, aes(label=paste0("p= ", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.9,vjust=1.5, size=4, colour="black")+
  labs(x="Substitutions", y="Insertion before position i")
p_before
ggsave(p_before, file="p_subs_ins_before_by pos.pdf", width = 12, height = 10, path=path)



corr_vector=c()
for(i in c(1:41)){
  corr<-cor.test(singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_single,
                 singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_ins_after, use="complete.obs")$estimate
  p_value<-cor.test(singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_single,
                    singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_ins_after, use="complete.obs")$p.value
  corr_vector=c(corr_vector,i,corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(1:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}


p_after<-ggplot(singles_insertions_dels, aes(x=nscore_single, y=nscore_ins_after))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0,linetype="dashed", color="darkgrey")+
  geom_point()+
  theme_bw()+
  facet_wrap(~as.numeric(Pos))+
  theme(panel.grid = element_blank())+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.15,vjust=1.5, size=4, colour="black")+
  geom_text(data=corr_text, aes(label=paste0("p= ", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.9,vjust=1.5, size=4, colour="black")+
  labs(x="Substitutions", y="Insertion after position i")
p_after
ggsave(p_after, file="p_subs_ins_after_by pos.pdf", width = 12, height = 10, path=path)












# compare average effects - by mutant AA- and then by NT/CT as well

my_vector=c()
for(i in all_aa){
  
  # this AA as a mutant in the substitutions
  
  mean_subs<-mean(singles.df[singles.df$Mut==i,]$nscore_c, na.rm = T)
  mean_subs_nt<-mean(singles.df[singles.df$Mut==i & singles.df$Pos %in% c(1:28),]$nscore_c, na.rm = T)
  mean_subs_ct<-mean(singles.df[singles.df$Mut==i & singles.df$Pos %in% c(29:42),]$nscore_c, na.rm = T)
  
  # this AA as a deletion
  
  this_pos<-as.vector(which(ABseq ==i))
  mean_del<-mean(single_deletions_reps[single_deletions_reps$del_pos %in% this_pos,]$nscore_c, na.rm=T)
  
  this_pos_nt<- this_pos[which(this_pos<=28)]
  mean_del_nt<-mean(single_deletions_reps[single_deletions_reps$del_pos %in% this_pos_nt,]$nscore_c, na.rm=T)
  
  this_pos_ct<- this_pos[which(this_pos>28)]
  mean_del_ct<-mean(single_deletions_reps[single_deletions_reps$del_pos %in% this_pos_ct,]$nscore_c, na.rm=T)
  
  # this AA as an insertion
  
  mean_ins<-mean(insertions_reps[insertions_reps$ins_aa==i,]$nscore_c, na.rm = T)
  mean_ins_nt<-mean(insertions_reps[insertions_reps$ins_aa==i & insertions_reps$ins_pos %in% c(1:28) ,]$nscore_c, na.rm = T)
  mean_ins_ct<-mean(insertions_reps[insertions_reps$ins_aa==i& insertions_reps$ins_pos %in% c(29:42),]$nscore_c, na.rm = T)
  
  
  my_vector=c(my_vector,i, mean_subs, mean_subs_nt, mean_subs_ct,mean_ins, mean_ins_nt, mean_ins_ct, mean_del, mean_del_nt, mean_del_ct)
  
}


means_df<-as.data.frame(matrix(my_vector, ncol=10, byrow = T))
colnames(means_df)<-c("AA", "mean_subs", "mean_subs_nt", "mean_subs_ct","mean_ins", "mean_ins_nt", "mean_ins_ct", "mean_del", "mean_del_nt", "mean_del_ct")
for(i in 2:length(means_df)){means_df[,i]<-as.numeric(as.character(means_df[,i]))}

means_df<-left_join(means_df, AA_type[,c("AA", "type")], by="AA")



list<-list( c("mean_subs","mean_del"),
                 c("mean_subs_nt","mean_del_nt"),
                 c("mean_subs_ct","mean_del_ct"),
                 
                 c("mean_subs","mean_ins"),
                 c("mean_subs_nt","mean_ins_nt"),
                 c("mean_subs_ct","mean_ins_ct"),
                 
                 c("mean_del","mean_ins"),
                 c("mean_del_nt","mean_ins_nt"),
                 c("mean_del_ct","mean_ins_ct"),
                 
                 c("mean_subs_nt","mean_subs_ct"),
                 c("mean_ins_nt","mean_ins_ct"),
                 c("mean_del_nt","mean_del_ct"))



plot_list<-list()

for(i in 1:length(list)){
  
  corr=cor.test(means_df[[list[[i]][1]]], means_df[[list[[i]][2]]], use="complete.obs")$estimate
  pval=cor.test(means_df[[list[[i]][1]]], means_df[[list[[i]][2]]], use="complete.obs")$p.value
  
  
  center_x<-(min(na.omit(means_df[[list[[i]][1]]]))+max(na.omit(means_df[[list[[i]][1]]])))/2
  #center_y<-(min(na.omit(means_df[[list[[i]][2]]]))+max(na.omit(means_df[[list[[i]][2]]])))/2
  limits_x=c((center_x-3), (center_x+3))
  #limits_y=c((center_y-3), (center_y+3))
  #limits_x=c(-5.3,2)
  limits_y=c(-5.3,2)
  
  
  p_corr_aa<-ggplot(means_df)+
    geom_hline(yintercept = 0, linetype="dashed", color="black")+
    geom_vline(xintercept = 0,linetype="dashed", color="black")+
    geom_point(aes( x=.data[[list[[i]][1]]], y=.data[[list[[i]][2]]], color=type))+
    geom_text_repel(aes( x=.data[[list[[i]][1]]], y=.data[[list[[i]][2]]], label=AA, color=type), show.legend = F)+
    theme_bw()+
    scale_color_manual(values=color_AAtype)+
    
    scale_x_continuous(limits = limits_x)+
    scale_y_continuous(limits = limits_y)+
    
        theme(panel.grid= element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())+
    labs(x=list[[i]][1], x=list[[i]][2],color="AA type")+
    annotate("text", label=paste0("R= ", round(corr,2)), x=Inf, y=-Inf,hjust=1.2, vjust=-1.8, size=5)+
    annotate("text", label=paste0("p= ", format(pval, digits = 2, scientific = T)), x=Inf, y=-Inf,hjust=1.2, vjust=-0.5, size=5)
  p_corr_aa
  
  plot_list[[i]]<-p_corr_aa
  
  
}


p_all<-ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]], 
                 plot_list[[4]],plot_list[[5]],plot_list[[6]],
                 plot_list[[7]],plot_list[[8]],plot_list[[9]], 
                 plot_list[[10]],plot_list[[11]],plot_list[[12]], ncol=3, nrow = 4, common.legend = T)
p_all

ggsave(p_all, file="p_subs_ins_dels_corr_by_aa_nt_ct.pdf", width=8, height = 11, path=path)





#adjusting X-axis limits for plots #5 and #6 x main figure

plot_list<-list()

for(i in c(5,6)){
  
  corr=cor.test(means_df[[list[[i]][1]]], means_df[[list[[i]][2]]], use="complete.obs")$estimate
  pval=cor.test(means_df[[list[[i]][1]]], means_df[[list[[i]][2]]], use="complete.obs")$p.value
  
  
  if(i==5){limits_x=c(-1.65,0.7)}
  if(i==6){limits_x=c(-4,0.2) } 
  limits_y=c(-5.3,2)
  
  
  p_corr_aa<-ggplot(means_df)+
    geom_hline(yintercept = 0, linetype="dashed", color="black")+
    geom_vline(xintercept = 0,linetype="dashed", color="black")+
    geom_point(aes( x=.data[[list[[i]][1]]], y=.data[[list[[i]][2]]], color=type))+
    geom_text_repel(aes( x=.data[[list[[i]][1]]], y=.data[[list[[i]][2]]], label=AA, color=type), show.legend = F)+
    theme_bw()+
    scale_color_manual(values=color_AAtype)+
    
    scale_x_continuous(limits = limits_x)+
    scale_y_continuous(limits = limits_y)+
    
    theme(panel.grid= element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())+
    labs(x=list[[i]][1], x=list[[i]][2],color="AA type")+
    annotate("text", label=paste0("R= ", round(corr,2)), x=Inf, y=-Inf,hjust=1.2, vjust=-1.8, size=5)+
    annotate("text", label=paste0("p= ", format(pval, digits = 2, scientific = T)), x=Inf, y=-Inf,hjust=1.2, vjust=-0.5, size=5)
  p_corr_aa
  
  plot_list[[i]]<-p_corr_aa
  
  
}

p_all<-ggarrange(plot_list[[5]],plot_list[[6]],
               ncol=1,common.legend = T)
p_all

ggsave(p_all, file="p_subs_ins_dels_corr_by_aa_nt_ct_main.pdf", width=3, height = 6, path=path)







# compare average effects - by position

positions<-as.vector(unique(singles_insertions_dels$Pos))
my_vector=c()
for(i in positions){
  
  mean_singles<-mean(singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_single, na.rm = T)
  mean_ins_before<-mean(singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_ins_before, na.rm = T)
  mean_ins_after<-mean(singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_ins_after, na.rm = T)
  nscore_del<-mean(singles_insertions_dels[singles_insertions_dels$Pos==i,]$nscore_del, na.rm = T)
  
  my_vector=c(my_vector,i, mean_singles, mean_ins_before, mean_ins_after, nscore_del)
  
}



means_df<-as.data.frame(matrix(my_vector, ncol=5, byrow = T))
colnames(means_df)<-c("Pos", "mean_single", "mean_ins_before", "mean_ins_after", "nscore_del")
for(i in 1:length(means_df)){means_df[,i]<-as.numeric(as.character(means_df[,i]))}


means_df$part<-""
means_df[means_df$Pos %in% c(1:28),]$part<-"NT"
means_df[means_df$Pos %in% c(29:42),]$part<-"CT"

means_df_nt<-means_df[means_df$part=="NT",]
means_df_ct<-means_df[means_df$part=="CT",]




list<-list(c("mean_single","mean_ins_before"),
           c("mean_single","mean_ins_after"),
           c("mean_single","nscore_del"),
           c("nscore_del","mean_ins_before"),
           c("nscore_del","mean_ins_after")
)


plot_list<-list()

for(i in 1:length(list)){
  
  corr=cor.test(means_df[[list[[i]][1]]], means_df[[list[[i]][2]]], use="complete.obs")$estimate
  corr_nt=cor.test(means_df_nt[[list[[i]][1]]], means_df_nt[[list[[i]][2]]], use="complete.obs")$estimate
  corr_ct=cor.test(means_df_ct[[list[[i]][1]]], means_df_ct[[list[[i]][2]]], use="complete.obs")$estimate
  
  pval=cor.test(means_df[[list[[i]][1]]], means_df[[list[[i]][2]]], use="complete.obs")$p.value
  pval_nt=cor.test(means_df_nt[[list[[i]][1]]], means_df_nt[[list[[i]][2]]], use="complete.obs")$p.value
  pval_ct=cor.test(means_df_ct[[list[[i]][1]]], means_df_ct[[list[[i]][2]]], use="complete.obs")$p.value
  
  p_corr_aa<-ggplot(means_df)+
    geom_hline(yintercept = 0, linetype="dashed", color="black")+
    geom_vline(xintercept = 0,linetype="dashed", color="black")+
    geom_point(aes( x=.data[[list[[i]][1]]], y=.data[[list[[i]][2]]], color=part))+
    theme_bw()+
    scale_color_manual("", values=c("grey70","grey40"))+
    
    scale_x_continuous(limits = c(-5.3,2))+
    scale_y_continuous(limits = c(-5.3,2))+
    
    theme(panel.grid= element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())+
    labs(x=list[[i]][1], x=list[[i]][2], color="")+
    annotate("text", label=paste0("R= ", round(corr,2)), x=Inf, y=-Inf,hjust=1.2, vjust=-10, size=5, color="black")+
    annotate("text", label=paste0("p= ", format(pval, digits = 2, scientific = T)), x=Inf, y=-Inf,hjust=1.2, vjust=-8.5, size=5, color="black")+
    annotate("text", label=paste0("R= ", round(corr_nt,2)), x=Inf, y=-Inf,hjust=1.2, vjust=-7, size=5, color="grey70")+
    annotate("text", label=paste0("p= ", format(pval_nt, digits = 2, scientific = T)), x=Inf, y=-Inf,hjust=1.2, vjust=-5.5, size=5, color="grey70")+
    annotate("text", label=paste0("R= ", round(corr_ct,2)), x=Inf, y=-Inf,hjust=1.2, vjust=-4, size=5, color="grey40")+
    annotate("text", label=paste0("p= ", format(pval_ct, digits = 2, scientific = T)), x=Inf, y=-Inf,hjust=1.2, vjust=-2.5, size=5, color="grey40")
  p_corr_aa
  
  plot_list[[i]]<-p_corr_aa
  
  
}

p_all<-ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]], 
                 plot_list[[4]],plot_list[[5]], common.legend = T)
p_all


ggsave(p_all, file="p_subs_ins_dels_corr_by_pos.pdf", width=8.5, height = 6, path=path)




