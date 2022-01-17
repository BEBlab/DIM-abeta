## Deletions


## required packages and functions
require(ggplot2)
require(stringr)
require(ggpubr)
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



dir.create("_Internal deletions")
path="_Internal deletions"


##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")




#########################################################################################################################################


AB_wt<-"DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
ABseq<-unlist(strsplit(AB_wt, ""))
ABseq_pos<-paste0(ABseq, "\n", c(1:42))
ABseq_pos_v<-paste0(ABseq, c(1:42))


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











# I add single AA deletions and truncations to internal deletions x visualisation 

deletions<-deletions_reps[,c("aa_seq","nscore_c", "sigma","p.adjust","category_fdr", "ID", "del_length", "del_start", "del_end", "dataset")]


truncations_to_add<-truncations.df[truncations.df$k_end==42 | truncations.df$k_start==1,
                                c("aa_seq","nscore_c", "sigma","p.adjust","category_fdr", "ID", "k_length", "k_start", "k_end", "dataset")]

truncations_to_add$del_length<-42-(as.numeric(truncations_to_add$k_length))

truncations_to_add$del_start<-0
truncations_to_add[truncations_to_add$k_end==42,]$del_start<-1
truncations_to_add[truncations_to_add$k_start==1,]$del_start<-(truncations_to_add[truncations_to_add$k_start==1,]$k_end)+1

truncations_to_add$del_end<-0
truncations_to_add[truncations_to_add$k_start==1,]$del_end<-42
truncations_to_add[truncations_to_add$k_end==42,]$del_end<-(truncations_to_add[truncations_to_add$k_end==42,]$k_start)-1


single_dels_to_add<-single_deletions_reps

single_dels_to_add$del_length<-1
single_dels_to_add$del_start<-single_dels_to_add$del_pos
single_dels_to_add$del_end<-single_dels_to_add$del_pos


deletions_map<-rbind(deletions, 
                     truncations_to_add[,c("aa_seq","nscore_c", "sigma","p.adjust","category_fdr", "ID", "del_length", "del_start", "del_end", "dataset")],
                     single_dels_to_add[,c("aa_seq","nscore_c", "sigma", "p.adjust", "category_fdr","ID", "del_length", "del_start", "del_end", "dataset")])

save(deletions_map, file="deletions_map.RData")


#load("deletions_map.RData")












## big matrix deletions

# add missing deletions

all_dels<-expand.grid(c(1:42), c(1:42))
colnames(all_dels)<-c("del_start", "del_end")
all_dels<-all_dels[all_dels$del_start<=all_dels$del_end,]
all_dels$del_length<-all_dels$del_end-all_dels$del_start+1
all_dels<-all_dels[all_dels$del_length<40,]
all_dels$ID<-paste0("Del_k", all_dels$del_length, "_", all_dels$del_start, "-", all_dels$del_end)

all_dels[all_dels$del_start==1 & all_dels$del_end!=1,]$ID<-paste0("kmer", 42-all_dels[all_dels$del_start==1& all_dels$del_end!=1,]$del_end, "_",
                                            all_dels[all_dels$del_start==1& all_dels$del_end!=1,]$del_end+1, "-42")

all_dels[all_dels$del_end==42 & all_dels$del_start!=42,]$ID<-paste0("kmer", 42-all_dels[all_dels$del_end==42& all_dels$del_start!=42,]$del_length, 
                                           "_1-",42-all_dels[all_dels$del_end==42& all_dels$del_start!=42,]$del_length)

missing_dels<-all_dels[!all_dels$ID %in% deletions_map$ID,]
missing_dels[,c("aa_seq", "nscore_c", "sigma", "p.adjust", "category_fdr")]<-NA
missing_dels$dataset<-"deletion"


heatmap_deletions_map<-rbind(deletions_map, missing_dels)



#add info syn and info on non-nucleating variants
heatmap_deletions_map$box<-"VUS"
heatmap_deletions_map[heatmap_deletions_map$ID %in% c("Del_k1_22-22", "Del_k6_19-24"),]$box<-"fAD"

heatmap_deletions_map$label<-""
heatmap_deletions_map[heatmap_deletions_map$ID %in% non_nuc_df$ID,]$label<-"-"





p_deletions<-ggplot(heatmap_deletions_map, aes(x=factor(del_start, levels=c(1:42)), 
                                               y=factor(del_end, levels = c(1:42)), fill=nscore_c))+
  theme_bw()+
  geom_tile(color="white")+
  geom_tile(data=heatmap_deletions_map[heatmap_deletions_map$box=="fAD",],color="black", fill=NA, size=1)+
  scale_fill_gradientn(colors=cols, limits=c(-min,max), na.value = "grey60")+
  labs(x="First deleted position", y="Last deleted position", fill="Nucleation\nscore")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x=factor(del_start, levels=c(1:42)), 
                y=factor(del_end, levels = c(1:42)),label=label), size=7)+
  scale_x_discrete(labels=ABseq_pos)+
  scale_y_discrete(labels=ABseq_pos_v, position = "left")

p_deletions

ggsave(p_deletions, file="p_deletions.pdf",width=9, height  = 7.5, path=path)





## brick plot for short deletions (up to 6AA deletion)

dels_map_max6<-deletions_map[as.numeric(deletions_map$del_length)<=6,]
dels_map_max6$row<-""

for(i in 1:nrow(dels_map_max6)){
  this_length<-as.numeric(dels_map_max6[i,]$del_length)
  
  # check which row this deletion belong to
  this_start<-as.numeric(dels_map_max6[i,]$del_start)
  
  while(this_start>this_length){this_start<-this_start-this_length}
  
  dels_map_max6[i,]$row<-paste(this_length, this_start, sep="_")
  
}


my_levels=c()
for(i in 1:6){
  for(j in 1:6){
    my_levels=c(my_levels,paste0(i,"_",j))      }}

my_labels=c()
for(i in 1:6){ my_labels=c(my_labels, rep(i, i))    }


p_short_dels<-ggplot(dels_map_max6,
      aes(xmin=as.numeric(del_start)-0.4,
            xmax=as.numeric(del_end)+0.4,
            ymin=factor(row, levels=my_levels),
            ymax=factor(row, levels=my_levels), 
            color=nscore_c))+
      theme_bw()+
      geom_rect( size=3)+
      scale_color_gradientn(colors=cols, limits=c(-min,max))+
      labs( y="Deletion length", color="Nucleation\nscore")+
      theme(panel.border = element_blank(),
            panel.grid.minor = element_blank())+
      scale_x_continuous(breaks=c(1:42), labels = ABseq_pos)+
      scale_y_discrete(labels=my_labels)
p_short_dels

ggsave(p_short_dels, file="p_short_deletions.pdf", height=3, width=7, path=path)





## FDR NS+ and NS- stacked barplot for each starting,ending and deleted position

#regions=c("all", "NT")
regions=c( "NT")

for(reg in regions){

  subset_df<-deletions_map
  
  if(reg=="NT"){
    
    subset_df<-deletions_map[as.numeric(deletions_map$del_end)<29,]
  }
  
    counts_start<-as.data.frame(subset_df %>% group_by(del_start, category_fdr) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)))
    counts_end<-as.data.frame(subset_df %>% group_by(del_end, category_fdr) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)))
    
    
    # deletions reps
    
    if(reg=="NT"){positions=c(1:28)}
    if(reg=="all"){positions=c(1:42)}
    
    my_vector<-c()
    for(i in positions){
      
      total<-nrow(subset_df[subset_df$del_start<=i & subset_df$del_end>=i,])
      
      for(cat in c("NS_inc", "NS_dec", "WT-like")){
      
      this_n<-nrow(subset_df[subset_df$del_start<=i & subset_df$del_end>=i & subset_df$category_fdr==cat,])
      
      my_vector<-c(my_vector, i,cat, this_n, this_n/total)
      
      
      }
    }
    
    df<-as.data.frame(matrix(my_vector, ncol=4, byrow = T))
    colnames(df)<-c("pos",  "category_fdr","n", "freq")
    
    df$freq<-as.numeric(as.character(df$freq))
    
    
    
    levels = c("NS_inc", "WT-like", "NS_dec")
    labels=c("NS+", "WT-like", "NS-")
    colors<-c( "darkblue", "grey90", "brown3")
    
    p_start<-ggplot(counts_start,  aes(x=factor(del_start, levels=positions), fill=factor(category_fdr, levels=levels, labels=labels),y=freq)) +
      geom_bar(position=position_fill(reverse = T), stat="identity", alpha=0.8, width=0.8)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            axis.line.y = element_line(color='black', size=0.25),
            axis.ticks.x = element_blank(),
            strip.background = element_rect(fill=NA, color=NA))+
      scale_fill_manual(values=colors)+
      scale_y_continuous(position = "right")+
      labs(x="First deleted position", y="Frequency", fill="FDR=0.1")
    p_start
    
    p_end<-ggplot(counts_end,  aes(x=factor(del_end, levels=positions),fill=factor(category_fdr, levels=levels , labels=labels),y=freq)) + 
      geom_bar(position=position_fill(reverse = T), stat="identity", alpha=0.8, width=0.8)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            axis.line.y = element_line(color='black', size=0.25),
            axis.ticks.x = element_blank(),
            strip.background = element_rect(fill=NA, color=NA))+
      scale_fill_manual(values=colors)+
      scale_y_continuous(position = "right")+
      labs(x="Last deleted position", y="Frequency", fill="FDR=0.1")
    p_end
    
    
    p_containing<-ggplot(df,  aes(x=factor(pos, levels=positions), fill=factor(category_fdr, levels=levels, labels=labels),y=freq)) +
      geom_bar(position=position_fill(reverse = T), stat="identity", alpha=0.8, width=0.8)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            axis.line.y = element_line(color='black', size=0.25),
            axis.ticks.x = element_blank(),
            strip.background = element_rect(fill=NA, color=NA))+
     scale_fill_manual(values=colors)+
      scale_y_continuous(position = "right")+
      labs(x="Containing deletion", y="Frequency", fill="FDR=0.1")
    p_containing
    
    
    
    p_all<-ggarrange(p_start, p_end,p_containing, ncol=1, common.legend = T, legend = "right")
    p_all
    ggsave(p_all, file=paste0("p_hotspot_freq_",reg, ".pdf"), height = 5, width = 8, path=path)
    




# define hotspot region for internal deletions
# hotspot contains those consecutive residues where [freq NS+ position i >  (max freq NS+ any position)/2]  at first and last deleted positions 

max_freq_start<-max(counts_start[counts_start$category_fdr=="NS_inc" ,]$freq)/2
top_pos_start<-unique(as.numeric(counts_start[counts_start$category_fdr=="NS_inc" & counts_start$freq>max_freq_start,]$del_start))
top_pos_start<-top_pos_start[order(top_pos_start, decreasing = F)]

consecutive_start<-split(top_pos_start, cumsum(c(1, diff(top_pos_start) != 1)))
hotspot_start<-as.vector(consecutive_start[[which.max(as.vector(unlist(lapply(consecutive_start, function(x) length(x)))))]])


max_freq_end<-max(counts_end[counts_end$category_fdr=="NS_inc" ,]$freq)/2 
top_pos_end<-unique(as.numeric(counts_end[counts_end$category_fdr=="NS_inc" & counts_end$freq>max_freq_end,]$del_end))
top_pos_end<-top_pos_end[order(top_pos_end, decreasing = F)]

consecutive_end<-split(top_pos_end, cumsum(c(1, diff(top_pos_end) != 1)))
hotspot_end<-as.vector(consecutive_end[[which.max(as.vector(unlist(lapply(consecutive_end, function(x) length(x)))))]])


print(paste0("region: ",reg, "_hotspot start: ", hotspot_start[1], "-", hotspot_start[length(hotspot_start)]))
print(paste0("region: ",reg, "_hotspot end: ", hotspot_end[1], "-", hotspot_end[length(hotspot_end)]))


}






## print sequences in the hotspot (NS+ and NS-)

hotspot_seqs_inc<-deletions_map[deletions_map$category_fdr=="NS_inc" & 
                                  deletions_map$del_start %in% hotspot_start & 
                                  deletions_map$del_end %in% hotspot_end,]
print(paste0("NS+ seqs in hotspot = ", length(unique(hotspot_seqs_inc$aa_seq))))




hotspot_seqs<-hotspot_seqs_inc


hotspot_seqs<-hotspot_seqs[order(hotspot_seqs$nscore_c, decreasing = T),]

my_vector=c()
for(i in 1:nrow(hotspot_seqs)){
 
  this_seq<-unlist(strsplit(hotspot_seqs[i,]$aa_seq,""))
  start<-as.numeric(hotspot_seqs[i,]$del_start)
  end<-as.numeric(hotspot_seqs[i,]$del_end)
  
  if(start!=1){new_seq<-c(this_seq[1:(start-1)], rep("-",(end-start+1)), this_seq[start:length(this_seq)] )  }
  if(start==1){new_seq<-c(rep("-",(end-start+1)), this_seq[start:length(this_seq)] )   }
  
  my_vector=c(my_vector,hotspot_seqs[i,]$aa_seq,hotspot_seqs[i,]$ID, new_seq)
}


seqs_df<-as.data.frame(matrix(my_vector, ncol=44, byrow = T))
colnames(seqs_df)<-c("aa_seq","ID", c(1:42))
melt_seqs_df<-melt(seqs_df, id=c("aa_seq", "ID")) 
melt_seqs_df<-left_join(melt_seqs_df, AA_type[,c("AA", "type")], by=c("value"="AA"))


levels=as.vector(unique(hotspot_seqs[order(hotspot_seqs$nscore_c),]$ID))

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


p_scores<-ggplot(hotspot_seqs)+
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
       file="p_written_seqs_hotspot.pdf", width = 10, height = 8, path=path)









# take deletions that remove ED

subset_df<-deletions_map[as.numeric(deletions_map$del_start)<24 & 
                             as.numeric(deletions_map$del_end)>21 & 
                             as.numeric(deletions_map$del_end)<29,]


subset_df$residue<-"no_charge"
subset_df$distance<-"no_charge"

for(i in 1:nrow(subset_df)){
  
  this_seq<-unlist(strsplit(subset_df[i,]$aa_seq, ""))
  no_core<-this_seq[1:(length(this_seq)-14)]
  
  if(length(which(no_core %in% c("D", "E"))==T)!=0){
    
    pos<-as.vector(which(no_core %in% c("D", "E")))
    pos<-pos[length(pos)]
    
    subset_df[i,]$residue<- no_core[pos]
    subset_df[i,]$distance<-(length(no_core)-pos+1)
    
  }
  
}


p_remove_ed<-ggplot(subset_df, aes(x=factor(distance, levels=c(1:42, "no_charge")), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_point(aes(shape=residue, 
                 color=factor(category_fdr, levels=c("NS_inc", "WT-like", "NS_dec"), labels=c("NS+", "WT-like", "NS-"))),
             size=1.8)+
  scale_color_manual("FDR=0.1", values=c("NS+"="darkblue", "WT-like"="grey", "NS-"="brown3"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line())+
  labs(x="D/E Distance from the C-terminus", y="Nucleation score")

p_remove_ed

ggsave(p_remove_ed, file="p_remove_ed.pdf", height = 3, width = 6, path=path)











