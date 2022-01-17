## alternative amyloid cores


## required packages and functions
require(ggplot2)
require(stringr)
require(reshape2)
require(dplyr)


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


dir.create("_alternative core")
path="_alternative core"




##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")
load("deletions_map.RData")





#########################################################################################################################################


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










# deletions in both NT and CT
entering_ct<-deletions_map[as.numeric(deletions_map$del_start)<29 &
                             as.numeric(deletions_map$del_start)!=1 &
                             as.numeric(deletions_map$del_end)>28 &
                             as.numeric(deletions_map$del_end)!=42,]

print(paste0("deletions in both N and C-terminus = ", length(unique(entering_ct$aa_seq))))
print(paste0("NS+ deletions in both N and C-terminus = ", length(unique(entering_ct[entering_ct$category_fdr=="NS_inc",]$aa_seq))))





## print sequences in both Nt and Ct and FDR NS+

entering_ct_inc<-entering_ct[entering_ct$category_fdr=="NS_inc",]
entering_ct_inc<-entering_ct_inc[order(entering_ct_inc$nscore_c, decreasing = T),]

my_vector=c()
for(i in 1:nrow(entering_ct_inc)){
  
  this_seq<-unlist(strsplit(entering_ct_inc[i,]$aa_seq,""))
  start<-as.numeric(entering_ct_inc[i,]$del_start)
  end<-as.numeric(entering_ct_inc[i,]$del_end)
  
  if(start!=1){new_seq<-c(this_seq[1:(start-1)], rep("-",(end-start+1)), this_seq[start:length(this_seq)] )}
  if(start==1){new_seq<-c(rep("-",(end-start+1)), this_seq[start:length(this_seq)] )}
  
  my_vector=c(my_vector,entering_ct_inc[i,]$aa_seq,entering_ct_inc[i,]$ID, new_seq)
}

seqs_df<-as.data.frame(matrix(my_vector, ncol=44, byrow = T))
colnames(seqs_df)<-c("aa_seq","ID", c(1:42))
melt_seqs_df<-melt(seqs_df, id=c("aa_seq", "ID")) 
melt_seqs_df<-left_join(melt_seqs_df, AA_type[,c("AA", "type")], by=c("value"="AA"))


levels=as.vector(unique(entering_ct_inc[order(entering_ct_inc$nscore_c),]$ID))

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


p_scores<-ggplot(entering_ct_inc)+
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
       file="p_written_seqs_both_nt_ct.pdf", width = 10, height = 2.5, path=path)





# take the part that is not the remaining core and check whats in there
# the first position of the remaining core will be the position where the deletion started in the original sequence

entering_ct$nt_tail<-""
entering_ct$ct_remaining<-""

entering_ct$nt_remaining_length<-""
entering_ct$ct_remaining_length<-""

for(i in 1:nrow(entering_ct)){
  
  unlist_seq<-(unlist(strsplit(entering_ct[i,]$aa_seq, "")))
  
  entering_ct[i,]$nt_tail<-paste0(unlist_seq[1:(as.numeric(entering_ct[i,]$del_start)-1)], collapse = "")
  entering_ct[i,]$ct_remaining<-paste0(unlist_seq[(as.numeric(entering_ct[i,]$del_start)):length(unlist_seq)], collapse = "")
  
  entering_ct[i,]$nt_remaining_length<-length(unlist_seq[1:(as.numeric(entering_ct[i,]$del_start)-1)])
  entering_ct[i,]$ct_remaining_length<-length(unlist_seq[(as.numeric(entering_ct[i,]$del_start)):length(unlist_seq)])
}

entering_ct$nt_remaining_length<-as.numeric(entering_ct$nt_remaining_length)
entering_ct$ct_remaining_length<-as.numeric(entering_ct$ct_remaining_length)






## length alternative cores

# take those ones that could be alternative cores, they keep K16 with some aliphatics from the LVFFA
# this means that deletions started at 18,19,20,21,22 

possible_cores<-entering_ct[entering_ct$del_start %in% c(18:22) ,]
possible_cores$length_after_K<-(possible_cores$nt_remaining_length) -16 + (possible_cores$ct_remaining_length)


p_possible_cores<-ggplot(possible_cores, 
                         aes(x=length_after_K, y=nscore_c, 
                              color=factor(category_fdr, levels=c("NS_inc", "WT-like", "NS_dec"), labels=c("NS+", "WT-like", "NS-"))))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_vline(xintercept = 14, size=0.1)+
  geom_point(size=2)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line())+
  scale_x_continuous(breaks = c(1:20))+
  scale_color_manual("FDR=0.1", values=c("NS+"="darkblue", "WT-like"="grey", "NS-"="brown3"))+
  labs(x="Length of the alternative core", y="Nucleation score")
p_possible_cores

ggsave(p_possible_cores, file="p_length_possible_cores.pdf", height = 3, width = 6, path=path)





## print possible cores that cannot nucleate and FDR NS-

possible_cores_dec<-possible_cores[possible_cores$category_fdr=="NS_dec",]
possible_cores_dec<-possible_cores_dec[order(possible_cores_dec$nscore_c, decreasing = T),]

my_vector=c()
for(i in 1:nrow(possible_cores_dec)){
  
  this_seq<-unlist(strsplit(possible_cores_dec[i,]$aa_seq,""))
  start<-as.numeric(possible_cores_dec[i,]$del_start)
  end<-as.numeric(possible_cores_dec[i,]$del_end)
  
  if(start!=1){new_seq<-c(this_seq[1:(start-1)], rep("-",(end-start+1)), this_seq[start:length(this_seq)] )}
  if(start==1){new_seq<-c(rep("-",(end-start+1)), this_seq[start:length(this_seq)] )}
  
  my_vector=c(my_vector,possible_cores_dec[i,]$aa_seq,possible_cores_dec[i,]$ID, new_seq)
}

seqs_df<-as.data.frame(matrix(my_vector, ncol=44, byrow = T))
colnames(seqs_df)<-c("aa_seq","ID", c(1:42))
melt_seqs_df<-melt(seqs_df, id=c("aa_seq", "ID")) 
melt_seqs_df<-left_join(melt_seqs_df, AA_type[,c("AA", "type")], by=c("value"="AA"))


levels=as.vector(unique(possible_cores_dec[order(possible_cores_dec$nscore_c),]$ID))

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


p_scores<-ggplot(possible_cores_dec)+
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
       file="p_written_seqs_possible_cores_dec.pdf", width = 10, height = 12, path=path)








