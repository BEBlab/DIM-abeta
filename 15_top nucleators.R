## Top 1% nucleators in the dataset


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




dir.create("_top nucleators")
path="_top nucleators"


##  import required data 
load("INDEL.df.RData")





#########################################################################################################################################


all_aa<-c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "-")

AA_type<-data.frame("AA"= all_aa,
                    "name_AA"=c("Glycine", "Alanine","Valine","Leucine","Methionine","Isoleucine","Phenylalanine",
                                "Tyrosine","Tryptophan","Lysine","Arginine","Aspartic acid","Glutamic acid","Serine","Threonine",
                                "Cysteine","Asparagine","Glutamine","Histidine","Proline", "deletion"),
                    "type"=c("glycine",rep("aliphatic",5),rep("aromatic",3),rep("positive",2),rep("negative",2),rep("polar",6),"proline", "deletion"))

color_AAtype=c("aliphatic"="darkgrey",
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












# take the 1% top scorers

ranked<-INDEL.df[order(INDEL.df$nscore_c, decreasing = T),]

how_many<-round(nrow(INDEL.df)*0.01)
ranked_top<-ranked[c(1:how_many),]


# print the sequences

ranked_top$new_ID<-""

my_vector=c()

for(i in 1:nrow(ranked_top)){
  this_seq<-unlist(strsplit(ranked_top[i,]$aa_seq,""))

  if(ranked_top[i,]$dataset=="truncation"){
  
    to_add<-42-length(this_seq)
    new_seq<-c(rep("-", to_add), this_seq, "")
    
    new_ID<-paste0("truncation ",to_add+1, "-", 42 )
    ranked_top[i,]$new_ID<-new_ID
    
    my_vector=c(my_vector,ranked_top[i,]$aa_seq,new_ID, new_seq)

  }
  
  
  if(ranked_top[i,]$dataset=="deletion"){
  
    ID<-ranked_top[i,]$ID
    
    if(!is.na(str_extract(ID, ";"))){  ID<-unlist(str_split(ID, ";"))[1]}
  
    start<-as.numeric(unlist(strsplit(unlist(strsplit(ID, "_"))[3], "-"))[1])
    end<-as.numeric(unlist(strsplit(unlist(strsplit(ID, "_"))[3], "-"))[2])
    new_seq<-c(this_seq[1:(start-1)], rep("-",(end-start+1)), this_seq[start:length(this_seq)] , "")
    
    new_ID<-paste0("deletion ", start, "-", end)
    ranked_top[i,]$new_ID<-new_ID
    
    my_vector=c(my_vector,ranked_top[i,]$aa_seq,new_ID, new_seq)
    
  }
  
  if(ranked_top[i,]$dataset=="single_deletion"){
    
    ID<-ranked_top[i,]$ID
    
    del_pos<-as.numeric(unlist(strsplit(unlist(strsplit(ID, "_"))[3], "-"))[1])
    new_seq<-c(this_seq[1:(del_pos-1)], "-", this_seq[(del_pos):length(this_seq)] , "")
    
    new_ID<-paste0("single deletion ", del_pos)
    ranked_top[i,]$new_ID<-new_ID
    
    my_vector=c(my_vector,ranked_top[i,]$aa_seq,new_ID, new_seq)
    
  }
  
  
  if(ranked_top[i,]$dataset=="single"){new_seq<-c(this_seq, "")
                                        new_ID<-paste0("substitution ",ranked_top[i,]$ID)
                                        ranked_top[i,]$new_ID<-new_ID
                                        
                                        my_vector=c(my_vector,ranked_top[i,]$aa_seq, new_ID, new_seq)
  }
 
   if(ranked_top[i,]$dataset=="insertion"){new_seq<-c(this_seq)
                                          new_ID<-paste0("insertion ",unlist(strsplit(ranked_top[i,]$ID, "_"))[3], unlist(strsplit(ranked_top[i,]$ID, "_"))[4])
                                          ranked_top[i,]$new_ID<-new_ID
                                          
                                          my_vector=c(my_vector,ranked_top[i,]$aa_seq,new_ID, new_seq)
  }
      
}
  
  


seqs_df<-as.data.frame(matrix(my_vector, ncol=45, byrow = T))
colnames(seqs_df)<-c("aa_seq","ID", c(1:43))
melt_seqs_df<-melt(seqs_df, id=c("aa_seq", "ID"))
melt_seqs_df<-left_join(melt_seqs_df, AA_type[,c("AA", "type")], by=c("value"="AA"))

labels=as.vector(unique(ranked_top[order(ranked_top$nscore_c),]$new_ID))
levels=as.vector(unique(ranked_top[order(ranked_top$nscore_c),]$aa_seq))


p_written_seqs<-ggplot(melt_seqs_df, aes(x=variable, y=factor(aa_seq, levels=levels, labels=labels)))+
  geom_text(aes(label=value, color=type))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
      axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left")+
  labs(x="", y="", color="AA type")+
  scale_color_manual(values=color_AAtype)
p_written_seqs


p_scores<-ggplot(ranked_top)+
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
       file="p_written_top_seqs.pdf", width = 10, height = 8, path=path)



