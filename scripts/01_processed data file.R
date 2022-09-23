## to obtain MS_BL_BB_indels_processed_data file from DiMSum output 


## required packages
require(dplyr)



## working directory contains 'required data' folder

##  import required data 
load("required data/indels_library.RData") 
load("required data/fAD.df.RData")


# load DiMSum output file
all_variants_df<-all_variants
variant_data_merge_df<-variant_data_merge


#########################################################################################################################################






#ID
all_variants_df<-left_join(all_variants_df, indels_library, by=c("aa_seq"="AA_sequence"))

#centering to wt
wt_nscore<-all_variants_df[which(all_variants_df$WT==T),]$fitness

all_variants_df$nscore_c<-as.numeric(paste(as.numeric(all_variants_df$fitness)+(-wt_nscore)))
all_variants_df$nscore1_c<-as.numeric(paste(as.numeric(all_variants_df$fitness1_uncorr)+(-wt_nscore)))
all_variants_df$nscore2_c<-as.numeric(paste(as.numeric(all_variants_df$fitness2_uncorr)+(-wt_nscore)))
all_variants_df$nscore3_c<-as.numeric(paste(as.numeric(all_variants_df$fitness3_uncorr)+(-wt_nscore)))

# add info on fAD
all_variants_df$fAD<-"non-fAD"
all_variants_df[all_variants_df$aa_seq %in% fAD.df[fAD.df$dominant_recessive=="D",]$aa_seq,]$fAD<-"fAD_d"
all_variants_df[all_variants_df$aa_seq %in% fAD.df[fAD.df$dominant_recessive=="R",]$aa_seq,]$fAD<-"fAD_r"

all_variants_df<-all_variants_df[,c("aa_seq","ID","dataset", "mean_count","nscore_c", "nscore1_c", "nscore2_c","nscore3_c", "sigma", "fAD")]








#find non-nucleating variants (they have input reads but 0 output reads- NS not calculated)
#each AA variant is only resulting from one nt sequence by design (non nuc nt seq results in non nuc aa seq)
non_nuc_df<-variant_data_merge_df[variant_data_merge_df$output1_e1_s1_b1_count==0  &
                                    variant_data_merge_df$output2_e2_s1_b1_count==0  &
                                    variant_data_merge_df$output3_e3_s1_b1_count==0 ,]

non_nuc_df<-as.data.frame(non_nuc_df[,c("aa_seq")])
colnames(non_nuc_df)<-"aa_seq"
non_nuc_df<-left_join(non_nuc_df, indels_library, by=c("aa_seq"="AA_sequence"))
non_nuc_df[,c("mean_count", "nscore_c", "nscore1_c", "nscore2_c","nscore3_c", "sigma", "fAD")]<-NA


all_variants_df<-rbind(all_variants_df, non_nuc_df)





write.table(all_variants_df, file="required data/MS_BL_BB_indels_processed_data.tsv", sep="\t", quote = F, row.names = F)








