## get chimera attribute files

## required packages
require(dplyr)


dir.create("_attribute_files")
path="_attribute_files"


##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")


#########################################################################################################################################








############ susbtitutions NS+ ############ 

singles<-singles.df

df<-as.data.frame(singles %>% group_by(Pos,category_fdr) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)))
df_inc<-df[df$category_fdr=="NS_inc",c("Pos", "freq")]


# add missing positions
seen_pos<-as.numeric(as.character(df_inc$Pos))
missing_pos<- c(1:42)[!c(1:42) %in% seen_pos]
to_add<-data.frame("Pos"=missing_pos, "freq"=0)

df_inc<-rbind(df_inc, to_add)

df_inc$row<-paste0("\t", ":", df_inc$Pos, "\t", df_inc$freq)

write_text<-c("attribute: subsinc\t", "recipient: residues\t", df_inc$row)
write.table(write_text, file=paste0(path, "/attrfile_subs_inc.txt"), quote = F, row.names = F, col.names = F)



############ susbtitutions NS- ############ 

df_dec<-df[df$category_fdr=="NS_dec",c("Pos", "freq")]

# add missing positions
seen_pos<-as.numeric(as.character(df_dec$Pos))
missing_pos<- c(1:42)[!c(1:42) %in% seen_pos]
to_add<-data.frame("Pos"=missing_pos, "freq"=0)

df_dec<-rbind(df_dec, to_add)

df_dec$row<-paste0("\t", ":", df_dec$Pos, "\t", df_dec$freq)

write_text<-c("attribute: subsdec\t", "recipient: residues\t", df_dec$row)
write.table(write_text, file=paste0(path, "/attrfile_subs_dec.txt"), quote = F, row.names = F, col.names = F)








############ insertions NS+ ############ 

insertions<-insertions_reps

df<-as.data.frame(insertions %>% group_by(ins_pos,category_fdr) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)))

df_inc<-df[df$category_fdr=="NS_inc",c("ins_pos", "freq")]

# add missing positions
seen_pos<-as.numeric(as.character(df_inc$ins_pos))
missing_pos<- c(1:41)[!c(1:41) %in% seen_pos]
to_add<-data.frame("ins_pos"=missing_pos, "freq"=0)

df_inc<-rbind(df_inc, to_add)

df_inc$row<-paste0("\t", ":", (as.numeric(df_inc$ins_pos)+1), "\t", df_inc$freq)

write_text<-c("attribute: insinc\t", "recipient: residues\t", df_inc$row)
write.table(write_text, file=paste0(path, "/attrfile_ins_inc.txt"), quote = F, row.names = F, col.names = F)




############ insertions NS- ############ 

df_dec<-df[df$category_fdr=="NS_dec",c("ins_pos", "freq")]

# add missing positions
seen_pos<-as.numeric(as.character(df_dec$ins_pos))
missing_pos<- c(1:41)[!c(1:41) %in% seen_pos]
to_add<-data.frame("ins_pos"=missing_pos, "freq"=0)

df_dec<-rbind(df_dec, to_add)

df_dec$row<-paste0("\t", ":", (as.numeric(df_dec$ins_pos)+1), "\t", df_dec$freq)

write_text<-c("attribute: insdec\t", "recipient: residues\t", df_dec$row)
write.table(write_text, file=paste0(path, "/attrfile_ins_dec.txt"), quote = F, row.names = F, col.names = F)








############ internal deletions NS+  ############ 

# exclude positions 1 and 42
#take only deletions at N-terminus

subset_df<-deletions_reps[as.numeric(deletions_reps$del_end)<29,]

my_vector<-c()
for(i in c(2:28)){
  total<-nrow(subset_df[subset_df$del_start<=i & subset_df$del_end>=i,])
  n_inc<-nrow(subset_df[subset_df$del_start<=i & subset_df$del_end>=i & subset_df$category_fdr=="NS_inc",])
  
  my_vector<-c(my_vector, i, n_inc, total, n_inc/total)
}


df<-as.data.frame(matrix(my_vector, ncol=4, byrow = T))
colnames(df)<-c("pos", "n_inc", "total", "perc")

df$row<-paste0("\t", ":", df$pos, "\t", df$perc)

write_text<-c("attribute: internaldelsinc\t", "recipient: residues\t", df$row)
write.table(write_text, file=paste0(path, "/attrfile_internaldels_inc_NT.txt"), quote = F, row.names = F, col.names = F)




############ internal deletions NS- ############ 

my_vector<-c()
for(i in c(2:27)){
  total<-nrow(subset_df[subset_df$del_start<=i & subset_df$del_end>=i,])
  n_inc<-nrow(subset_df[subset_df$del_start<=i & subset_df$del_end>=i & subset_df$category_fdr=="NS_dec",])
  
  my_vector<-c(my_vector, i, n_inc, total, n_inc/total)
}

df<-as.data.frame(matrix(my_vector, ncol=4, byrow = T))
colnames(df)<-c("pos", "n_inc", "total", "perc")

df$row<-paste0("\t", ":", df$pos, "\t", df$perc)

write_text<-c("attribute: internaldelsdec\t", "recipient: residues\t", df$row)
write.table(write_text, file=paste0(path, "/attrfile_internaldels_dec_NT.txt"), quote = F, row.names = F, col.names = F)









############ single deletions ############ 

df<-single_deletions_reps[,c("nscore_c", "del_pos")]

df$row<-paste0("\t", ":", df$del_pos, "\t", df$nscore_c)

write_text<-c("attribute: singledelsNS\t", "recipient: residues\t", df$row)
write.table(write_text, file=paste0(path, "/attrfile_single_dels_NS.txt"), quote = F, row.names = F, col.names = F)








############ N-terminal truncations NS+ ############ 

df<-truncations.df[truncations.df$side=="Nterminal",c("nscore_c", "k_start")]
df$pos<-(df$k_start)-1

df$row<-paste0("\t", ":", df$pos, "\t", df$nscore_c)

write_text<-c("attribute: truncNS\t", "recipient: residues\t", df$row)
write.table(write_text, file=paste0(path, "/attrfile_trunc_NS.txt"), quote = F, row.names = F, col.names = F)












####### CHIMERA command line ####### 

# % NS+ and % NS- (substitutions, insertions and N-terminal multi AA deletions)


##  rangecolor attrfile_inc 0 #dddddd 0.2 #C6C6D3 0.4 #9A9ABF 0.6 #424299 1 #00007C novalue #f2f2f2 
##  rangecolor attrfile_dec 0 #dddddd 0.2 #D5C9C9 0.4 #C6A2A2 0.6 #A75454 1 #911A1A novalue #f2f2f2

# example color side chains
## rangecolor attrfile_inc 0 #dddddd 0.2 #C6C6D3 0.4 #9A9ABF 0.6 #424299 1 #00007C novalue #f2f2f2; color #0c5bb0,a :lys,arg; color #ee0011,a :asp,glu; color #a9a9a9,a :ala,leu,ile,met,val; color #9a703e,a :phe,tyr; color #15983d,a :his,ser,asn,gln


# plot NS in single AA deletions and N-terminal truncations 
#min<-abs(min(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))
#max<-abs(max(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))

##  rangecolor atrrfile_NS min #CD3333 0 #dddddd max #00008B novalue #f2f2f2



