## ROC curves and list of fAD candidate variants


## required packages
require(dplyr)
require(pROC)
require(ggplot2)
require(RColorBrewer)



dir.create("_ROC")
path="_ROC"



##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")
load("required data/loadings_predictors.RData")



#########################################################################################################################################










roc_df<-left_join(singles.df[,c("aa_seq", "nscore_c")],loadings_predictors, by="aa_seq")


disease_mutations<-INDEL.df[INDEL.df$fAD=="fAD_d" & INDEL.df$dataset=="single",]$aa_seq
roc_df$fAD<-0
roc_df[roc_df$aa_seq %in% disease_mutations,]$fAD<-1


print(paste0("n(fAD)=", nrow(roc_df[roc_df$fAD==1,])))
print(paste0("n(non-fAD)=",nrow(roc_df[roc_df$fAD==0,])))




predictors<-c("Camsol","Zyggregator","Tango","Solubility","Waltz","Polyphen","CADD","nscore_c")

#build ROC with pROC package

mylist_fAD12 = list()
my_vector_auc = c()

for(i in predictors){
  subset<-roc_df[,c("fAD",i)]
  subset<-na.omit(subset)
  
  glm.fit=glm(subset$fAD~subset[[i]], family = binomial)
  roc<-roc(subset$fAD, glm.fit$fitted.values)
  mylist_fAD12[[i]] <- roc
  
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

g_fAD12<-ggroc(mylist_fAD12[this_order],legacy.axes = T, size=1)+
            theme_bw()+
            theme(legend.title = element_blank(),
                    legend.text = element_text(size=12),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(), 
                    axis.line = element_line(color='black'),
                    axis.title = element_text(size=20),
                    axis.text = element_text(size=15))+
            labs(x="False Positive Rate", y="True Positive Rate")+
            geom_abline(linetype="dashed")+
            scale_color_brewer(labels=auc$label, palette = "Set1")
g_fAD12

ggsave(g_fAD12, path=path, file="p_ROC.pdf", width = 7, height = 4)






# candidate fAD variants, exlude those that are already described as fAD

exclude<-c(INDEL.df[INDEL.df$fAD!="non-fAD",]$ID)

candidate_fad<-INDEL.df[INDEL.df$category_fdr=="NS_inc" & !INDEL.df$ID %in% exclude ,c("aa_seq", "ID", "dataset", "nscore_c", "category_fdr")]
candidate_fad<-candidate_fad[order(candidate_fad$nscore_c, decreasing = T),]
print(paste0("n candidate fAD=", nrow(candidate_fad)))

write.table(candidate_fad, file=paste0(path, "/candidate_fad.tsv"), sep="\t", quote = F, row.names = F)
