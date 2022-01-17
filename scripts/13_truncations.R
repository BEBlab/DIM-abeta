## Truncations


## required packages and functions
require(ggplot2)
require(stringr)
require(reshape2)
require(dplyr)
require(ggpubr)

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




dir.create("_Truncations")
path="_Truncations"


##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")




#########################################################################################################################################



AB_wt<-"DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
ABseq<-unlist(strsplit(AB_wt, ""))
ABseq_pos<-paste0(ABseq, "\n", c(1:42))
ABseq_pos_v<-paste0(ABseq, c(1:42))

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














# I add single AA deletions of positions 1 and 42 to truncations x visualisation 

to_add<-single_deletions_reps[single_deletions_reps$del_pos %in% c("1", "42"),]

to_add<-select(to_add,!c("del_pos"))
to_add$k_length<-"41"

to_add$k_start<-"1"
to_add[to_add$ID=="Del_k1_1-1",]$k_start<-"2"

to_add$k_end<-"41"
to_add[to_add$ID=="Del_k1_1-1",]$k_end<-"42"

to_add$side<-"Cterminal"
to_add[to_add$ID=="Del_k1_1-1",]$side<-"Nterminal"

truncations_map<-rbind(truncations.df, to_add)




### all truncations

p_truncations<-ggplot(truncations_map, 
                      aes(x=factor(k_start, levels=c(1:42)), y=factor(k_end, levels=c(1:42)), fill=nscore_c))+
  theme_bw()+
  geom_tile(color="black")+
  scale_fill_gradientn(colors=cols, limits=c(-min,max), na.value = "grey60")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_discrete(labels=ABseq_pos, breaks=c(1:42))+
  scale_y_discrete(labels=ABseq_pos_v, breaks=c(1:42))+
  labs(x="First position of the peptide", y="Last position of the peptide", fill="Nucleation\nscore")
p_truncations

ggsave(p_truncations, file="p_truncations_both_sides.pdf",width=8, height=6.5, path=path)




### Nterminal truncations
truncations_nt<-truncations_map[truncations_map$side=="Nterminal",]

p_trunc_nt<-ggplot(truncations_nt, aes(x=factor(k_start, levels=c(1:42)), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.1, color="black")+
  geom_errorbar(aes(ymin=nscore_c-1.96*sigma, ymax=nscore_c+1.96*sigma),width=0, size=0.1)+
 
  geom_point(data=truncations_nt[truncations_nt$category_fdr %in% c("NS_inc", "NS_dec"),],
             aes(fill=nscore_c),size=4, shape=21, stroke=1.2)+
  
  geom_point(data=truncations_nt[!truncations_nt$category_fdr %in% c("NS_inc", "NS_dec"),],
             aes(fill=nscore_c),size=4, shape=21, stroke=0.2)+
  
  theme_bw()+
  theme(panel.border = element_blank(),
    panel.grid.major = element_line(size=0.2), 
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color='black'))+
  labs(x="First position of the peptide", y="Nucleation score", fill="Nucleation score")+
  scale_y_continuous(limits = c(-5.4,3.6), breaks = c(-4,-2,0,2))+
  scale_fill_gradientn(colors=cols,limits=c(-min, max))
p_trunc_nt




### Cterminal truncations
truncations_ct<-truncations_map[truncations_map$side=="Cterminal",]


p_trunc_ct<-ggplot(truncations_ct, aes(x=factor(k_end, levels=c(1:42)), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.1, color="black")+
  geom_errorbar(aes(ymin=nscore_c-1.96*sigma, ymax=nscore_c+1.96*sigma), width=0, size=0.1)+
  
  geom_point(data=truncations_ct[truncations_ct$category_fdr %in% c("NS_inc", "NS_dec"),],
             aes(fill=nscore_c),size=4, shape=21, stroke=1.2)+
  
  geom_point(data=truncations_ct[!truncations_ct$category_fdr %in% c("NS_inc", "NS_dec"),],
             aes(fill=nscore_c),size=4, shape=21, stroke=0.2)+

    theme_bw()+
  theme(panel.border = element_blank(),
    panel.grid.major = element_line(size=0.2), 
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color='black'))+
  labs(x="Last position of the peptide", y="Nucleation score",fill="Nucleation score")+
  scale_y_continuous(limits = c(-5.4,3.6), breaks = c(-4, -2, 0, 2))+
  scale_fill_gradientn(colors=cols, limits=c(-min,max))
p_trunc_ct


p_trunc<-ggarrange(p_trunc_nt, p_trunc_ct, nrow = 2, align="hv", common.legend = T, legend = "bottom")
p_trunc

ggsave(p_trunc, file="p_truncations_one_side.pdf", width = 8, height = 5, path=path)





# print NS+ truncations

inc_trunc<-truncations_map[truncations_map$category_fdr=="NS_inc",]

my_vector=c()
for(i in 1:nrow(inc_trunc)){
  this_seq<-unlist(strsplit(inc_trunc[i,]$aa_seq,""))
  to_add<-42-length(this_seq)
  new_seq<-c(rep("", to_add), this_seq)
  my_vector=c(my_vector,inc_trunc[i,]$aa_seq, new_seq)
}

seqs_df<-as.data.frame(matrix(my_vector, ncol=43, byrow = T))
colnames(seqs_df)<-c("aa_seq", c(1:42))
melt_seqs_df<-melt(seqs_df, id="aa_seq")
melt_seqs_df<-left_join(melt_seqs_df, AA_type[,c("AA", "type")], by=c("value"="AA"))

levels=as.vector(inc_trunc[order(inc_trunc$nscore_c),]$aa_seq)

p_written_seqs<-ggplot(melt_seqs_df, aes(x=variable, y=factor(aa_seq, levels=levels)))+
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


p_scores<-ggplot(inc_trunc)+
  geom_tile(aes(x=1, y=factor(aa_seq, levels=levels), fill=nscore_c), color="white", size=0.5)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank())+
  scale_fill_gradientn(colors=cols, limits=c(-min,max))+
  labs(x="", y="", fill="Nucleation\nscore")
p_scores



ggsave(lay_out(list(p_written_seqs, 1, 1:5),list(p_scores, 1, 6)), 
       file="p_written_truncations_inc_NS.pdf", width = 10, height = 4, path=path)


