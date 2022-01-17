## FDR heatmaps for substitutions, insertions and deletions (which will include single dels, truncations and internal dels)


## required packages and functions
require(ggplot2)
require(stringr)



dir.create("_FDR heatmaps")
path="_FDR heatmaps"


##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")
load("deletions_map.RData")



#########################################################################################################################################


AB_wt<-"DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
ABseq<-unlist(strsplit(AB_wt, ""))
ABseq_pos<-paste0(ABseq, "\n", c(1:42))
ABseq_pos_v<-paste0(ABseq, c(1:42))

all_aa=c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

disease_mutations<-INDEL.df[INDEL.df$fAD=="fAD_d" & INDEL.df$dataset=="single",]$ID










####  single aa substitutions

single_subst=c()
for (i in c(1:42)){
  wt_AA<-unlist(strsplit(AB_wt, ""))[i]
  
  for (j in all_aa){
    new_mutation<-paste0(wt_AA,"-", i, "-", j)
    single_subst<-c(single_subst, new_mutation)
    
  }}


singles<-singles.df[,c("aa_seq", "ID", "nscore_c", "sigma", "p.adjust", "category_fdr")]

missing_singles<-data.frame("aa_seq"=NA,
                            "ID"=single_subst[!single_subst %in% singles$ID],
                            "nscore_c"=NA,
                            "sigma"=0,
                            "p.adjust"=NA,
                            "category_fdr"=NA)

singles<-rbind(singles, missing_singles)

singles$ID<-as.character(singles$ID)

singles$WT_AA<-""
singles$Pos<-""
singles$Mut<-""

for(i in 1:nrow(singles)){
  singles[i,]$WT_AA<-unlist(strsplit(singles[i,]$ID, "-"))[1]
  singles[i,]$Pos<-unlist(strsplit(singles[i,]$ID, "-"))[2]
  singles[i,]$Mut<-unlist(strsplit(singles[i,]$ID, "-"))[3]
}

singles$Pos<-as.numeric(singles$Pos)



#add info fAD
singles$box<-"VUS"
singles[singles$ID %in% disease_mutations,]$box<-"Dominant"


#add info syn
singles$label<-""
singles[singles$WT_AA == singles$Mut,]$label<-"*"

singles[singles$WT_AA == singles$Mut,]$nscore_c<-0



heatmap_fdr<-singles

heatmap_fdr$category<-"WT-like"
heatmap_fdr[heatmap_fdr$WT_AA == heatmap_fdr$Mut,]$category<- "syn"
heatmap_fdr[is.na(heatmap_fdr$nscore_c),]$category<-"missing"

heatmap_fdr$p.adjust<-as.numeric(heatmap_fdr$p.adjust)

heatmap_fdr[(!is.na(heatmap_fdr$aa_seq) & heatmap_fdr$p.adjust<0.25 & heatmap_fdr$nscore_c<0),]$category<- "NS- 25%"
heatmap_fdr[(!is.na(heatmap_fdr$aa_seq) & heatmap_fdr$p.adjust<0.1 & heatmap_fdr$nscore_c<0),]$category<- "NS- 10%"
heatmap_fdr[(!is.na(heatmap_fdr$aa_seq) & heatmap_fdr$p.adjust<0.05 & heatmap_fdr$nscore_c<0),]$category<- "NS- 5%"
heatmap_fdr[(!is.na(heatmap_fdr$aa_seq) & heatmap_fdr$p.adjust<0.01 & heatmap_fdr$nscore_c<0),]$category<- "NS- 1%"

heatmap_fdr[(!is.na(heatmap_fdr$aa_seq) & heatmap_fdr$p.adjust<0.25 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 25%"
heatmap_fdr[(!is.na(heatmap_fdr$aa_seq) & heatmap_fdr$p.adjust<0.1 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 10%"
heatmap_fdr[(!is.na(heatmap_fdr$aa_seq) & heatmap_fdr$p.adjust<0.05 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 5%"
heatmap_fdr[(!is.na(heatmap_fdr$aa_seq) & heatmap_fdr$p.adjust<0.01 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 1%"



ramp_colors <- c(colorRampPalette(c( "brown3", "grey95", "darkblue"))(9))
colors<-c(ramp_colors, "white", "grey30" )

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%", "syn", "missing")

p_fdr_heatmap<-ggplot(heatmap_fdr)+
  geom_tile(aes(x=Pos,y=factor(Mut, levels=rev(all_aa)),fill=factor(category, levels=levels)), size=0.1, color="white")+
  geom_tile(data=singles[singles$box=="Dominant",], aes(x=Pos,y=factor(Mut, levels=rev(all_aa))),color="black", fill=NA, size=1)+
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
  scale_fill_manual("Category (FDR)", values= colors)
p_fdr_heatmap

ggsave(p_fdr_heatmap, file="p_fdr_heatmap_subst.pdf",width=14, height=8, path=path)










#### single AA insertions

# use insertion_reps dataframe to visualise duplicated coding sequences

insertions<-insertions_reps[,c("ins_pos", "ins_aa", "nscore_c", "sigma", "p.adjust")]

my_ins<-c(paste0(insertions$ins_pos,"-", insertions$ins_aa))

#add missing insertions to heatmap
missing_ins=c()
for (i in c(1:41)){
  for (j in all_aa){
    
    this_ins<-paste0(i,"-",j)
    
    if(! this_ins %in% my_ins){  missing_ins=c(missing_ins,i,j)}
  }}


my_df<-as.data.frame(matrix(missing_ins, ncol=2, byrow = T))

missing_insertions<-data.frame("ins_pos"=my_df$V1,
                               "ins_aa"=my_df$V2,
                               "nscore_c"=NA,
                               "sigma"=NA,
                               "p.adjust"=NA)

heatmap_insertions_fdr<-rbind(insertions, missing_insertions)


heatmap_insertions_fdr$category<-"WT-like"
heatmap_insertions_fdr[is.na(heatmap_insertions_fdr$nscore_c),]$category<-"missing"

heatmap_insertions_fdr$p.adjust<-as.numeric(heatmap_insertions_fdr$p.adjust)

heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing" &  heatmap_insertions_fdr$p.adjust<0.25 & heatmap_insertions_fdr$nscore_c<0),]$category<- "NS- 25%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.1 & heatmap_insertions_fdr$nscore_c<0),]$category<- "NS- 10%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.05 & heatmap_insertions_fdr$nscore_c<0),]$category<- "NS- 5%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.01 & heatmap_insertions_fdr$nscore_c<0),]$category<- "NS- 1%"

heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.25 & heatmap_insertions_fdr$nscore_c>0),]$category<- "NS+ 25%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.1 & heatmap_insertions_fdr$nscore_c>0),]$category<- "NS+ 10%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.05 & heatmap_insertions_fdr$nscore_c>0),]$category<- "NS+ 5%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.01 & heatmap_insertions_fdr$nscore_c>0),]$category<- "NS+ 1%"




ramp_colors <- c(colorRampPalette(c( "brown3", "grey95", "darkblue"))(9))
colors<-c(ramp_colors,  "grey30" )

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%", "missing")


p_fdr_heatmap_ins<-ggplot(heatmap_insertions_fdr)+
  geom_tile(aes(x=factor(ins_pos, levels=c(1:41), labels=c(2:42)),
                y=factor(ins_aa, levels=rev(all_aa)),fill=factor(category, levels=levels)), size=0.1, color="white")+
  theme_minimal()+
  theme()+
  labs(x="Position of inserted AA", y="Inserted amino acid", fill="Category (FDR)")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x.top = element_line(),       
        plot.title =element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=15, face = "bold"),
        legend.text = element_text(size=15), 
        axis.text = element_text( size=12),
        axis.title = element_text(size = 20))+
  scale_fill_manual(values= colors)
p_fdr_heatmap_ins

ggsave(p_fdr_heatmap_ins, file="p_fdr_heatmap_insertions.pdf",width=14, height=8, path=path)









#### deletions (including truncations and single deletions)


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
missing_dels[,c("aa_seq", "nscore_c", "p.adjust")]<-NA
missing_dels$dataset<-"deletion"


deletions_fdr<-rbind(deletions_map[,colnames(missing_dels)], missing_dels)


deletions_fdr$category<-"WT-like"
deletions_fdr[is.na(deletions_fdr$aa_seq),]$category<-"missing"

deletions_fdr$p.adjust<-as.numeric(deletions_fdr$p.adjust)

deletions_fdr[(deletions_fdr$category!="missing" &  deletions_fdr$p.adjust<0.25 & deletions_fdr$nscore_c<0),]$category<- "NS- 25%"
deletions_fdr[(deletions_fdr$category!="missing"  & deletions_fdr$p.adjust<0.1 & deletions_fdr$nscore_c<0),]$category<- "NS- 10%"
deletions_fdr[(deletions_fdr$category!="missing"  & deletions_fdr$p.adjust<0.05 & deletions_fdr$nscore_c<0),]$category<- "NS- 5%"
deletions_fdr[(deletions_fdr$category!="missing"  & deletions_fdr$p.adjust<0.01 & deletions_fdr$nscore_c<0),]$category<- "NS- 1%"

deletions_fdr[(deletions_fdr$category!="missing"  & deletions_fdr$p.adjust<0.25 & deletions_fdr$nscore_c>0),]$category<- "NS+ 25%"
deletions_fdr[(deletions_fdr$category!="missing"  & deletions_fdr$p.adjust<0.1 & deletions_fdr$nscore_c>0),]$category<- "NS+ 10%"
deletions_fdr[(deletions_fdr$category!="missing"  & deletions_fdr$p.adjust<0.05 & deletions_fdr$nscore_c>0),]$category<- "NS+ 5%"
deletions_fdr[(deletions_fdr$category!="missing"  & deletions_fdr$p.adjust<0.01 & deletions_fdr$nscore_c>0),]$category<- "NS+ 1%"





ramp_colors <- c(colorRampPalette(c( "brown3", "grey95", "darkblue"))(9))
colors<-c(ramp_colors,  "grey30" )

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%", "missing")


#add info fAD
deletions_fdr$box<-"VUS"
deletions_fdr[deletions_fdr$ID %in% c("Del_k1_22-22", "Del_k6_19-24"),]$box<-"fAD"


p_deletions_fdr<-ggplot(deletions_fdr, aes(x=factor(del_start, levels=c(1:42)), y=factor(del_end, levels = c(1:42)), 
                                       fill=factor(category, levels=levels)))+
  theme_bw()+
  geom_tile(color="white")+
  geom_tile(data=deletions_fdr[deletions_fdr$box=="fAD",],color="black", fill=NA, size=1)+
  scale_fill_manual(values=colors)+
  labs(x="First deleted position", y="Last deleted position", fill="Category (FDR)")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_discrete(labels=ABseq_pos)+
  scale_y_discrete(labels=ABseq_pos_v, position = "left")

p_deletions_fdr

ggsave(p_deletions_fdr, file="p_deletions_fdr.pdf",width=9, height  = 7, path=path)











###############
# truncations #
###############

# labelled as kmers in library design
# they don't have repetitions- they are unique. those labelled also as deletions are now only truncations; just remove the deletions name


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





truncations_map$category<-"WT-like"

#no missing truncations
#truncations_map[is.na(truncations_map$aa_seq),]$category<-"missing"

truncations_map$p.adjust<-as.numeric(truncations_map$p.adjust)

truncations_map[truncations_map$p.adjust<0.25 & truncations_map$nscore_c<0,]$category<- "NS- 25%"
truncations_map[truncations_map$p.adjust<0.1 & truncations_map$nscore_c<0,]$category<- "NS- 10%"
truncations_map[truncations_map$p.adjust<0.05 & truncations_map$nscore_c<0,]$category<- "NS- 5%"
truncations_map[truncations_map$p.adjust<0.01 & truncations_map$nscore_c<0,]$category<- "NS- 1%"

truncations_map[truncations_map$p.adjust<0.25 & truncations_map$nscore_c>0,]$category<- "NS+ 25%"
truncations_map[truncations_map$p.adjust<0.1 & truncations_map$nscore_c>0,]$category<- "NS+ 10%"
truncations_map[truncations_map$p.adjust<0.05 & truncations_map$nscore_c>0,]$category<- "NS+ 5%"
truncations_map[truncations_map$p.adjust<0.01 & truncations_map$nscore_c>0,]$category<- "NS+ 1%"




ramp_colors <- c(colorRampPalette(c( "brown3", "grey95", "darkblue"))(9))
colors<-c(ramp_colors,  "grey30" )

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%", "missing")



p_truncations_fdr<-ggplot(truncations_map, 
                      aes(x=factor(k_start, levels=c(1:42)), y=factor(k_end, levels=c(1:42)), fill=factor(category, levels=levels)))+
  theme_bw()+
  geom_tile(color="black")+
  scale_fill_manual(values=colors)+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_discrete(labels=ABseq_pos, breaks=c(1:42))+
  scale_y_discrete(labels=ABseq_pos_v, breaks=c(1:42))+
  labs(x="First position of the peptide", y="Last position of the peptide", fill="Category (FDR)")
p_truncations_fdr

ggsave(p_truncations_fdr, file="p_truncations_fdr.pdf",width=8, height=6.5, path=path)






