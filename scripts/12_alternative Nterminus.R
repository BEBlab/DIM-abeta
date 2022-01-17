## new N-terminus sequences and residues flanking the amyloid core


## required packages and functions
require(ggplot2)
require(stringr)
require(reshape2)
require(dplyr)
require(DescTools)


lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 



dir.create("_alternative nt")
path="_alternative nt"




##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")
load("deletions_map.RData")
load("required data/core_charge.RData")





#########################################################################################################################################


AB_wt<-"DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
ABseq<-unlist(strsplit(AB_wt, ""))
ABseq_pos<-paste0(ABseq, "\n", c(1:42))
ABseq_pos_v<-paste0(ABseq, c(1:42))
ABseq_pos_and_dist<-paste0(c(28:1),"\n", ABseq[1:28], "\n", c(1:28))


all_aa<-c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "-")

AA_type<-data.frame("AA"= all_aa,
                    "name_AA"=c("Glycine", "Alanine","Valine","Leucine","Methionine","Isoleucine","Phenylalanine",
                                "Tyrosine","Tryptophan","Lysine","Arginine","Aspartic acid","Glutamic acid","Serine","Threonine",
                                "Cysteine","Asparagine","Glutamine","Histidine","Proline", "deletion"),
                    "type"=c("glycine",rep("aliphatic",5),rep("aromatic",3),rep("positive",2),rep("negative",2),rep("polar",6),"proline", "deletion"))

color_AAtype<-c("aliphatic"="darkgrey",
                "aromatic"="#9A703EFF",
                "negative"="#EE0011FF",
                "positive"="#0C5BB0FF", 
                "polar"="#15983DFF", 
                "glycine"="grey30", 
                "proline"="#FEC10BFF",
                "deletion"="darkgrey")

## NS color gradient 
min<-abs(min(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))
max<-abs(max(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))

cols <- c(colorRampPalette(c( "brown3", "grey95"))((min/(min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(min+max)*100)-0.5))














# median NS at each position for NT deletions - new NT tails

nt_dels<-deletions_map[deletions_map$del_end<=28,]


#I will create a matrix with all positions
all_seqs=c()
for(i in 1:nrow(nt_dels)){
  this_seq<-unlist(strsplit(nt_dels[i,]$aa_seq, ""))
  to_add<-42-length(this_seq)
  this_seq<-c(rep("del", to_add), this_seq)
  
  all_seqs=c(all_seqs, this_seq)
}

matrix_seqs<-as.data.frame(matrix(all_seqs, ncol=42, byrow = T))
colnames(matrix_seqs)<-paste0("pos_", c(1:42))

nt_dels_matrix<-cbind(nt_dels, matrix_seqs[1:28])
melt_nt_dels_matrix<-melt(nt_dels_matrix, id=c(colnames(nt_dels)))








all_pos<-as.vector(unique(melt_nt_dels_matrix$variable))
all_mut<-as.vector(unique(melt_nt_dels_matrix$value))
all_mut<-all_mut[all_mut!="del"]

wt_df<-data.frame("pos"=all_pos,"wt_AA"=unlist(strsplit(AB_wt, ""))[1:28] )

my_vector=c()
for(aa in all_mut){
  for(pos in all_pos){
    
    this_median<-median(melt_nt_dels_matrix[melt_nt_dels_matrix$variable==pos & melt_nt_dels_matrix$value==aa,]$nscore_c )
    
    wt<-""
    if(aa==wt_df[wt_df$pos==pos,]$wt_AA){wt<-"wt"}  
    
    my_vector=c(my_vector, this_median, aa, pos, wt)
  }}


medians_df<-as.data.frame(matrix(my_vector, ncol=4, byrow = T))
colnames(medians_df)<-c("median", "aa", "pos", "wt")
medians_df$median<-as.numeric(as.character(medians_df$median))
medians_df<-medians_df[!is.na(medians_df$median),]

#p_median_heatmap<-ggplot(medians_df, aes(x=factor(pos, levels=all_pos, labels=ABseq_pos_and_dist), 
#                                         y=factor(aa, levels=rev(all_aa))))+
#  geom_tile(aes(fill=median), color="black")+
#  geom_text(aes(label=factor(wt, labels = c("", "*"))))+
#  scale_fill_gradientn(colours=cols, limits=c(-min,max))+
#  theme_bw()+
#  theme(panel.border = element_blank())+
#  scale_y_discrete(position="right")+
#  labs(x="Distance from the Ct\nWT aa and position", y="AA", fill="Median\nNucleation\nscore")
#p_median_heatmap

#ggsave(p_median_heatmap, file="p_median_heatmap.pdf", height = 4, width = 8, path=path)







# violin grid

melt_nt_dels_matrix$is_wt<-F
for(i in 1:nrow(melt_nt_dels_matrix)){
  this_pos<-ABseq[as.numeric(unlist(strsplit(as.character(melt_nt_dels_matrix[i,]$variable), "_"))[2])]
  if(this_pos==as.character(melt_nt_dels_matrix[i,]$value)){melt_nt_dels_matrix[i,]$is_wt<-T}
}

melt_nt_dels_matrix$ID_new<-paste0(melt_nt_dels_matrix$variable, "_", melt_nt_dels_matrix$value)
medians_df$ID_new<-paste0(medians_df$pos, "_", medians_df$aa)

melt_nt_dels_matrix<- left_join(melt_nt_dels_matrix, medians_df[,c("ID_new", "median")], by="ID_new")

p<-ggplot(melt_nt_dels_matrix[melt_nt_dels_matrix$value!="del",], aes(x=factor(variable, levels=all_pos, labels=ABseq_pos_and_dist), y=nscore_c))+
  facet_grid(factor(value, levels = rev(all_aa))~factor(variable, levels=all_pos, labels=ABseq_pos_and_dist), scales = "free_x")+
  geom_hline(yintercept = 0, size=0.2, color="black")+
  geom_violin(aes(fill=median))+
  geom_point(size=0.2)+
  labs(x="Distance from the Ct\nWT amino acid and position", y="Nucleation score")+
  theme_bw()+
  scale_fill_gradientn(colours=cols, limits=c(-min,max))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank())
#p

ggsave(p, file="p_violin_matrix.pdf", width = 12, height = 10, path=path)












# print seqs deleting K28

seqs_del_28<-deletions_map[deletions_map$del_end==28 & deletions_map$nscore_c>0,]
length(unique(seqs_del_28$aa_seq))

seqs_del_28<-seqs_del_28[order(seqs_del_28$nscore_c, decreasing = T),]


my_vector=c()
for(i in 1:nrow(seqs_del_28)){
  
  this_seq<-unlist(strsplit(seqs_del_28[i,]$aa_seq,""))
  start<-as.numeric(seqs_del_28[i,]$del_start)
  end<-as.numeric(seqs_del_28[i,]$del_end)
  
  if(start!=1){  new_seq<-c(this_seq[1:(start-1)], rep("-",(end-start+1)), this_seq[start:length(this_seq)] )}
  if(start==1){new_seq<-c(rep("-",(end-start+1)), this_seq[start:length(this_seq)] )}
  
  my_vector=c(my_vector,seqs_del_28[i,]$aa_seq,seqs_del_28[i,]$ID, new_seq)
}


seqs_df<-as.data.frame(matrix(my_vector, ncol=44, byrow = T))
colnames(seqs_df)<-c("aa_seq","ID", c(1:42))
melt_seqs_df<-melt(seqs_df, id=c("aa_seq", "ID")) 
melt_seqs_df<-left_join(melt_seqs_df, AA_type[,c("AA", "type")], by=c("value"="AA"))


levels=as.vector(unique(seqs_del_28[order(seqs_del_28$nscore_c),]$ID))

p_written_seqs<-ggplot(melt_seqs_df, aes(x=variable, y=factor(ID, levels=levels)))+
  geom_text(aes(label=value, color=type))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left")+
  labs(x="", y="", color="AA type")+
  scale_color_manual(values=color_AAtype)
p_written_seqs


p_scores<-ggplot(seqs_del_28)+
  geom_tile(aes(x=1, y=factor(ID, levels=levels), fill=nscore_c), color="white", size=0.5)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank())+
  scale_fill_gradientn(colors=cols, limits=c(-min,max))+
  labs(x="", y="", fill="Nucleation\nscore")
p_scores


ggsave(lay_out(list(p_written_seqs, 1, 1:5),list(p_scores, 1, 6)), 
       file="p_written_seqs_delete28.pdf", width = 10, height = 2, path=path)









# individual testing core-charge charge-core constructs

core_charge_melt<-melt(core_charge, id=c("ID"))


levels_core=c("Core","SupN","Abeta","K-Core","R-Core","D-Core","E-Core","Core-K","Core-R","Core-D","Core-E")
duntest<-DunnettTest(core_charge_melt$value, factor(core_charge_melt$ID, levels=levels_core), na.rm=T)

levels_supn=c("SupN","Core","Abeta","K-Core","R-Core","D-Core","E-Core","Core-K","Core-R","Core-D","Core-E")
duntest_supn<-DunnettTest(core_charge_melt$value, factor(core_charge_melt$ID, levels=levels_supn), na.rm=T)

levels_abeta=c("Abeta","Core","SupN","K-Core","R-Core","D-Core","E-Core","Core-K","Core-R","Core-D","Core-E")
duntest_abeta<-DunnettTest(core_charge_melt$value, factor(core_charge_melt$ID, levels=levels_abeta), na.rm=T)



duntest<-as.data.frame(duntest$Core)
duntest$ID<-str_remove(row.names(duntest), "-Core")

#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
duntest$sig<-""
duntest[duntest$pval<0.05,]$sig<-"*"
duntest[duntest$pval<0.01,]$sig<-"**"
duntest[duntest$pval<0.001,]$sig<-"***"


levels_plot=c("SupN","Core","Abeta", "K-Core","R-Core","D-Core","E-Core","Core-K","Core-R","Core-D","Core-E")
labels_plot=c("Sup35N","core","AB42", "K-core","R-core","D-core","E-core","core-K","core-R","core-D","core-E")

p_core_charge<-ggplot(core_charge_melt, aes(x=factor(ID, levels=levels_plot, labels=labels_plot), y=value))+
    theme_bw()+
    geom_boxplot(aes(color=factor(ID, levels=levels_plot, labels=labels_plot)), show.legend = F, width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=factor(ID, levels=levels_plot, labels=labels_plot)),width=0.2, size=1.5, show.legend = F)+
    #scale_color_manual(values=c("grey",rep("#ef4d59",2), rep("darkblue", 4), rep("cornflowerblue", 4)))+
  scale_color_manual(values=c("grey",rep("#ffd93d",2), rep("#8736aa", 8)))+
 
    
  theme(axis.title.x=element_blank(),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x=element_text(angle=45, vjust = 1, hjust=1))+
  geom_text(data=duntest, aes(x=factor(ID, levels=levels_plot, labels=labels_plot), y=55, label=sig), size=6)+
  labs(y="% -Adenine growth")
p_core_charge

ggsave(p_core_charge,file="p_core_charge.pdf", width=5, height=2.5, path=path)
  
  
