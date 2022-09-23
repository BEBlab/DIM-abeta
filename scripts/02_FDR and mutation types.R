## FDR and mutation type categories


## required packages
require(readxl)
require(dplyr)
require(stringr)


##  import required data 
load("required data/fAD.df.RData")
load("required data/indels_library.RData")
INDEL.df<-as.data.frame(read_excel("required data/MS_BL_BB_indels_processed_data.xlsx", sheet="all_variants"))


#########################################################################################################################################




# make sure the fAD are well annotated
INDEL.df$fAD<-"non-fAD"
INDEL.df[INDEL.df$aa_seq %in% fAD.df[fAD.df$dominant_recessive=="D",]$aa_seq,]$fAD<-"fAD_d"
INDEL.df[INDEL.df$aa_seq %in% fAD.df[fAD.df$dominant_recessive=="R",]$aa_seq,]$fAD<-"fAD_r"




#keep the variants with no NS estimate separately- they have input reads but no output reads
non_nuc_df<-INDEL.df[is.na(INDEL.df$nscore_c),]

INDEL.df<-INDEL.df[!is.na(INDEL.df$nscore_c),]



for(i in c("nscore_c", "nscore1_c", "nscore2_c", "nscore3_c", "sigma")){INDEL.df[[i]]<-as.numeric(INDEL.df[[i]])}

# test variants against WT at FDR=0.1
INDEL.df$zscore<-INDEL.df$nscore_c/INDEL.df$sigma
INDEL.df$p.adjust<-p.adjust(2*pnorm(-abs(INDEL.df$zscore)), method = "BH")

INDEL.df$sig_fdr<-FALSE
INDEL.df[INDEL.df$p.adjust<0.1,]$sig_fdr<-TRUE

INDEL.df$category_fdr<-"WT-like"
INDEL.df[INDEL.df$sig_fdr==T & INDEL.df$nscore_c<0,]$category_fdr<-"NS_dec"
INDEL.df[INDEL.df$sig_fdr==T & INDEL.df$nscore_c>0,]$category_fdr<-"NS_inc"






# count number of variants in the library
print(paste0("n (designed variants)= ", nrow(indels_library)))
print(as.data.frame(indels_library %>% group_by(dataset) %>% dplyr::summarise(n=n()) ))

# count number of variants with NS
print(paste0("n (scored variants)= ", nrow(INDEL.df)))
print(as.data.frame(INDEL.df %>% group_by(dataset) %>% dplyr::summarise(n=n()) ))
print(as.data.frame(INDEL.df %>% group_by(dataset, category_fdr) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq=n / sum(n)*100) ))



#divide into mutation types: single AA substitutions, single AA insertions, single AA deletions, truncations and internal deletions
singles.df<-INDEL.df[INDEL.df$dataset=="single",]
insertions.df<-INDEL.df[INDEL.df$dataset=="insertion",]
deletions.df<-INDEL.df[INDEL.df$dataset=="deletion",]
truncations.df<-INDEL.df[INDEL.df$dataset=="truncation",]
single_deletions.df<-INDEL.df[INDEL.df$dataset=="single_deletion",]




###########
# singles #
###########

singles.df$WT_AA<-""
singles.df$Pos<-""
singles.df$Mut<-""

for(i in 1:nrow(singles.df)){
  singles.df[i,]$WT_AA<-unlist(strsplit(singles.df[i,]$ID, "-"))[1]
  singles.df[i,]$Pos<-unlist(strsplit(singles.df[i,]$ID, "-"))[2]
  singles.df[i,]$Mut<-unlist(strsplit(singles.df[i,]$ID, "-"))[3]
}

singles.df$Pos<-as.numeric(singles.df$Pos)





##############
# insertions #
##############

# one coding sequence can come from more than one insertion- duplicate rows

insertions_reps<-insertions.df

insertions_reps$ins_pos<-""
insertions_reps$ins_aa<-""

for (i in 1:nrow(insertions_reps)){
  
  if(is.na(str_extract(insertions_reps[i,]$ID, ";"))){
    
    ID<-unlist(str_split(insertions_reps[i,]$ID, "_"))
    
    insertions_reps[i,]$ins_pos<-ID[3]
    insertions_reps[i,]$ins_aa<-ID[4]
    
  }else{
    
    all_ID<-unlist(str_split(insertions_reps[i,]$ID, ";"))
    num_ID<-length(all_ID)
    
    for(j in all_ID){
      
      new_row<- insertions_reps[i,]
      new_row$ID<-j
      
      ID<-unlist(str_split(j, "_"))
      
      new_row$ins_pos<-ID[3]
      new_row$ins_aa<-ID[4]
      
      insertions_reps<-rbind(insertions_reps,new_row)
      
    }
  }
}

#remove rows with more than one ID (they are already split)
insertions_reps<-insertions_reps[is.na(str_extract(insertions_reps$ID, ";")),]






#######################
# single AA deletions #
#######################

# one coding sequence can come from more than one single deletion- duplicate rows

single_deletions_reps<-single_deletions.df

single_deletions_reps$del_pos<-""

for (i in 1:nrow(single_deletions_reps)){
  
  if(is.na(str_extract(single_deletions_reps[i,]$ID, ";"))){
    
    ID<-unlist(str_split(single_deletions_reps[i,]$ID, "_"))
    
    single_deletions_reps[i,]$del_pos<-unlist(strsplit(ID[3],"-"))[1]
    
  }else{
    
    all_ID<-unlist(str_split(single_deletions_reps[i,]$ID, ";"))
    
    for(j in all_ID){
      
      new_row<- single_deletions_reps[i,]
      new_row$ID<-j
      
      ID<-unlist(str_split(j, "_"))
      
      new_row$del_pos<-unlist(strsplit(ID[3],"-"))[1]
      
      single_deletions_reps<-rbind(single_deletions_reps,new_row)
      
    }
  }
}

#remove rows with more than one ID (they are already split)
single_deletions_reps<-single_deletions_reps[is.na(str_extract(single_deletions_reps$ID, ";")),]
#also those that were kmers (=truncations)
single_deletions_reps<-single_deletions_reps[is.na(str_extract(single_deletions_reps$ID, "kmer")),]








###############
# truncations #
###############

# labelled as kmers in library design
# they don't have repetitions- they are unique. those labelled also as deletions are now only truncations; just remove the deletions name

truncations.df$k_length<-""
truncations.df$k_start<-""
truncations.df$k_end<-""

for (i in 1:nrow(truncations.df)){
  
  if(is.na(str_extract(truncations.df[i,]$ID, ";"))){
    
    ID<-unlist(str_split(truncations.df[i,]$ID, "_"))
    
    truncations.df[i,]$k_length<-gsub("kmer", "", ID[1])
    
    start_end<-unlist(str_split(ID[2], "-"))
    truncations.df[i,]$k_start<-start_end[1]
    truncations.df[i,]$k_end<-start_end[2]
    
  }else{
    
    all_ID<-unlist(str_split(truncations.df[i,]$ID, ";"))
    
    #remove IDs that are dels
    all_ID<-all_ID[is.na(str_extract(all_ID, "Del"))]
    
    num_ID<-length(all_ID)
    
    for(j in all_ID){
      
      new_row<- truncations.df[i,]
      new_row$ID<-j
      
      ID<-unlist(str_split(j, "_"))
      
      new_row$k_length<-gsub("kmer", "", ID[1])
      
      start_end<-unlist(str_split(ID[2], "-"))
      new_row$k_start<-start_end[1]
      new_row$k_end<-start_end[2]
      
      
      truncations.df<-rbind(truncations.df,new_row)
      
    }
  }
}


truncations.df$k_start<-as.numeric(truncations.df$k_start)
truncations.df$k_end<-as.numeric(truncations.df$k_end)

#remove rows with more than one ID (they are already split)
truncations.df<-truncations.df[is.na(str_extract(truncations.df$ID, "Del")),]


#also add from which side the truncation starts
truncations.df$side<-"both"
truncations.df[truncations.df$k_end==42 , ]$side<-"Nterminal"
truncations.df[truncations.df$k_start==1 , ]$side<-"Cterminal"







######################
# internal deletions #
######################


## create dataframes with repeated rows when they have more than one ID
deletions_reps<-deletions.df

for(i in 1:nrow(deletions_reps)){
  
  ID<-deletions_reps[i,]$ID
  
  if(!is.na(str_extract(ID, ";"))){
    
    all_ID<-unlist(str_split(ID, ";"))
    
    #remove IDs that are kmers
    all_ID<-all_ID[is.na(str_extract(all_ID, "kmer"))]
    
    for(j in all_ID){
      
      new_row<- deletions_reps[i,]
      new_row$ID<-j
      
      deletions_reps<-rbind(deletions_reps,new_row)
      
    }}}

#remove rows with more than one ID (they are already split)
deletions_reps<-deletions_reps[is.na(str_extract(deletions_reps$ID, ";")),]


# add info on start, end and length of deletion
deletions_reps$del_length<-0
deletions_reps$del_start<-0
deletions_reps$del_end<-0

for(i in 1:nrow(deletions_reps)){
  
  ID<-deletions_reps[i,]$ID
  
  positions<-unlist(strsplit(ID, "_"))[3]
  
  start<-as.numeric(unlist(strsplit(positions, "-"))[1])
  end<-as.numeric(unlist(strsplit(positions, "-"))[2])
  length<-end-start+1
  
  deletions_reps[i,]$del_start<-start
  deletions_reps[i,]$del_end<-end
  deletions_reps[i,]$del_length<-length
  
}







save(INDEL.df, file = "INDEL.df.RData")
save(singles.df, insertions.df, insertions_reps, single_deletions.df,  single_deletions_reps, 
     truncations.df,deletions.df,deletions_reps,  file="INDEL_datasets.RData")
save(non_nuc_df, file="required data/non_nucleating_variants.RData")


