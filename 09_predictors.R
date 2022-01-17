## Correlation predictors and AA properties with NS


## required packages
require(ggplot2)
require(dplyr)
require(reshape2)
require(ggrepel)
require(pROC)
require(weights)


dir.create("_Predictors")
path="_Predictors"


##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")

load("required data/loadings_predictors.RData")
load("required data/loadings_properties.RData")


#########################################################################################################################################





colors=c("Tango"="#E41A1C", 
         "Camsol"="#377EB8",
         "Zyggregator"="#4DAF4A",
         "Waltz"="#984EA3",
         "CADD"="#F781BF",
         "Polyphen"="#FFFF33",
         "Solubility"="#A65628",
         "PC1"="#FF7F00",
         "hydrophobicity"="#999999" )





# calculate PC1 and hydrophobicity scores for each single AA substitution
df<-singles.df[,c("aa_seq","Pos", "WT_AA", "Mut", "nscore_c", "sigma", "category_fdr", "ID", "dataset")]

df$hydrophobicity<-""
df$PC1<-""

for(i in 1:nrow(df)){
  WT_score<-loadings_properties[loadings_properties$amino_acid==df[i,]$WT_AA,]$`Hydrophobicity (Kyte-Doolittle)`  
  Mut_score<-loadings_properties[loadings_properties$amino_acid==df[i,]$Mut,]$`Hydrophobicity (Kyte-Doolittle)`
  df[i,]$hydrophobicity<-paste(Mut_score-WT_score)
  
  WT_score<-loadings_properties[loadings_properties$amino_acid==df[i,]$WT_AA,]$PC1  
  Mut_score<-loadings_properties[loadings_properties$amino_acid==df[i,]$Mut,]$PC1
  df[i,]$PC1<-paste(Mut_score-WT_score)
  
}

df$hydrophobicity<-as.numeric(df$hydrophobicity)
df$PC1<-as.numeric(df$PC1)

# swap sign of PC1 so it correlates with hydrophobicity
df$PC1<-(-df$PC1)




# take predictors table, merge NS, PC1 and hydrophobicity
singlemut_predictors_df<-left_join(loadings_predictors, INDEL.df[,c("aa_seq", "ID", "nscore_c", "sigma", "category_fdr", "dataset")], by="aa_seq")
singlemut_predictors_df<-left_join(singlemut_predictors_df, df[,c("aa_seq", "PC1", "hydrophobicity")], by="aa_seq")







### for each predictor and property, correlate data with NS and build ROC curve for NS+ variants classification

# info on mutations location
subst_nt<-singles.df[singles.df$Pos %in% c(1:28),]$aa_seq 
ins_nt<-insertions_reps[insertions_reps$ins_pos %in% c(1:28),]$aa_seq
singledels_nt<-single_deletions_reps[single_deletions_reps$del_pos %in% c(1:28),]$aa_seq



#regions=c("NT", "CT", "all")
regions=c("NT", "CT")
datasets=c("single","insertion", "single_deletion")
predictors_subst=c( "Tango", "Camsol", "Zyggregator", "Solubility", "Waltz","CADD", "Polyphen", "PC1", "hydrophobicity")
predictors_indels=c("Tango", "Camsol", "Zyggregator")


for(data in datasets){
  for(reg in regions){
    
    
      if(reg=="NT" & data=="single"){ df<-singlemut_predictors_df[singlemut_predictors_df$dataset == data & singlemut_predictors_df$aa_seq %in% subst_nt,]}
      if(reg=="CT" & data=="single"){ df<-singlemut_predictors_df[singlemut_predictors_df$dataset == data & !singlemut_predictors_df$aa_seq %in% subst_nt,]}
      if(reg=="all" & data=="single"){ df<-singlemut_predictors_df[singlemut_predictors_df$dataset == data ,]}
    
      if(reg=="NT" & data=="insertion"){ df<-singlemut_predictors_df[singlemut_predictors_df$dataset == data & singlemut_predictors_df$aa_seq %in% ins_nt,]}
      if(reg=="CT" & data=="insertion"){ df<-singlemut_predictors_df[singlemut_predictors_df$dataset == data & !singlemut_predictors_df$aa_seq %in% ins_nt,]}
      if(reg=="all" & data=="insertion"){ df<-singlemut_predictors_df[singlemut_predictors_df$dataset == data ,]}
    
      if(reg=="NT" & data=="single_deletion"){ df<-singlemut_predictors_df[singlemut_predictors_df$dataset == data & singlemut_predictors_df$aa_seq %in% singledels_nt,]}
      if(reg=="CT" & data=="single_deletion"){ df<-singlemut_predictors_df[singlemut_predictors_df$dataset == data & !singlemut_predictors_df$aa_seq %in% singledels_nt,]}
      if(reg=="all" & data=="single_deletion"){ df<-singlemut_predictors_df[singlemut_predictors_df$dataset == data ,]}
        
      
  corr_vector=c()
  pval_vector=c()
  
  if(data=="single"){predictors<-predictors_subst}
  if(data!="single"){predictors<-predictors_indels}
  
  for(i in predictors){
    
    subset<-subset.data.frame(df[,c("nscore_c","sigma",i )])
    subset<-na.omit(subset)
    subset<- subset[!is.infinite(rowSums(subset)),]
    
    
    ## weighted correlation
    
    pearson<-wtd.cor(subset$nscore_c, subset[,i], weight=subset$sigma^-2)
    corr<-pearson[1,1]
    p.value<-pearson[1,4]
    
    corr_vector=c(corr_vector,corr)
    pval_vector=c(pval_vector, as.numeric(p.value))
    
    
  }


  corr_text <- data.frame(
    label = corr_vector,
    variable=predictors)
  
  pval_text <- data.frame(
    label = as.numeric(pval_vector),
    variable=predictors)


melt_predictors<-melt(df, id=c("aa_seq", "nscore_c", "sigma", "category_fdr", "ID", "dataset"))

df_plot<-melt_predictors[ melt_predictors$variable %in% predictors,]

df_plot<-na.omit(df_plot)
df_plot<- df_plot[!is.infinite(df_plot$value),]

if (data=="single_deletion"){
  p<-ggplot(df_plot,aes(x=value, y=nscore_c))+
    geom_hline(yintercept = 0, size=0.5, linetype="dashed")+
    geom_point(color="grey50", size=3)+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size=0.1))+
    facet_wrap(~variable, scales = "free_x")+
    geom_text(data=corr_text, aes(label=paste0("R=",round(label, 2)), x=-Inf, y=Inf),hjust=-0.15,vjust=1.5, size=4, colour="black")+
    geom_text(data=pval_text, aes(label=paste0("p=",format(label, digits = 2, scientific = T)), x=-Inf, y=Inf),hjust=-0.1, vjust=3, size=4, colour="black")+
    labs(x="predictor score", y="Nucleation score")
  p
  
  
}else{
  
  
  
p<-ggplot(df_plot,aes(x=value, y=nscore_c))+
    geom_hline(yintercept = 0, size=0.5, linetype="dashed")+
    stat_binhex()+
  scale_fill_viridis_c()+
  #geom_point(color="grey60")+
  theme_bw()+
    theme(panel.grid.minor  = element_blank(),
          strip.background = element_blank(),
            panel.border = element_blank(),
          axis.line = element_line(size=0.1))+
    facet_wrap(~variable, scales = "free_x")+
    geom_text(data=corr_text, aes(label=paste0("R=",round(label, 2)), x=-Inf, y=Inf),hjust=-0.15,vjust=1.5, size=4, colour="black")+
    geom_text(data=pval_text, aes(label=paste0("p=",format(label, digits = 2, scientific = T)), x=-Inf, y=Inf),hjust=-0.1, vjust=3, size=4, colour="black")+
    labs(x="predictor score", y="Nucleation score", title=paste0(data, "_", reg))
  p
  
  
} 

if(data=="single"){ggsave(p, file=paste0("p_predictors_subst_", reg, ".pdf"), height = 8, width = 8, path=path)}
if(data!="single"){ggsave(p, file=paste0("p_predictors_", data,"_", reg , ".pdf"), height = 3, width = 8, path=path)}







## and ROC to see how do they perform in classifying NS+ variants

# can't perform ROC in CT single deletions because none of them is NS+
if(data!= "single_deletion" | (data=="single_deletion" & reg %in% c("NT", "all") )){

roc_df<-df


roc_df$agg_variants<-0
roc_df[roc_df$category_fdr =="NS_inc",]$agg_variants<-1


#build ROC with pROC package

mylist = list()
my_vector_auc = c()


for(i in predictors){
  subset<-roc_df[,c("agg_variants",i)]
  subset<-na.omit(subset)
  
  glm.fit=glm(subset$agg_variants~subset[[i]], family = binomial)
  roc<-roc(subset$agg_variants, glm.fit$fitted.values)
  mylist[[i]] <- roc
  
  auc<-roc$auc
  my_vector_auc=c(my_vector_auc, auc, i)
  
  
 
  
  
}



#extract AUC
auc<-as.data.frame(matrix(my_vector_auc, ncol=2, byrow=T))
colnames(auc)<-c("AUC", "Predictor")
auc$AUC<-round(as.numeric(as.character(auc$AUC)), 2)
auc$label<-paste0(auc$Predictor," (AUC=", auc$AUC, ")" )

#reorder by decreasing AUC
auc<-auc[ order(auc$AUC, decreasing = T ),]

#also keep the order for the plot list
this_order<-as.vector(auc$Predictor)





p_agg<-ggroc(mylist[this_order],legacy.axes = T, size=0.8)+
  theme_bw()+
  geom_abline(linetype="dashed")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black'),
        axis.title = element_text(size=20),
        axis.text = element_text(size=15))+
  labs(x="False Positive Rate", y="True Positive Rate", title=paste0(data, "_", reg))+
  scale_color_manual(values=colors, labels=auc$label)
p_agg

ggsave(p_agg, path=path, file=paste0("p_ROC_", data,"_", reg , ".pdf"), width = 6.8, height = 4)


}


  }}



