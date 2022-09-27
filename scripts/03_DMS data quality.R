## Reproducibility and validation of the DMS assay


## required packages
require(ggplot2)
require(stringr)
require(ggpubr)
require(dplyr)
require(reshape2)
require(weights)
require(DescTools)
require(ggrepel)


dir.create("_DMS data quality")
path="_DMS data quality"


##  import required data 
load("INDEL.df.RData")
load("INDEL_datasets.RData")

load("required data/error_prone_pcr_singles.RData")
load("required data/small_scale_validation.RData")
load("required data/kinetics.RData")
load("required data/toxicity_validation.RData")


#########################################################################################################################################










### replicate correlations

combinations<-c("1_2", "1_3", "2_3")

my_list<-list()

for(comb in combinations){
  
  first_rep<-unlist(strsplit(comb, "_"))[1]
  second_rep<-unlist(strsplit(comb, "_"))[2]
  
  subset<-INDEL.df[!is.na(INDEL.df[[paste0("nscore", first_rep, "_c")]]),]
  subset<-subset[!is.na(subset[[paste0("nscore",second_rep, "_c")]]),]
  
  print(paste0("n(",comb,")= ",length(subset$ID)))
  
  corr<-cor.test(INDEL.df[[paste0("nscore", first_rep, "_c")]], 
                 INDEL.df[[paste0("nscore", second_rep, "_c")]], use="complete.obs")
  R<-corr$estimate
  p<-corr$p.value
  
  p_corr<-ggplot(INDEL.df, aes(x=.data[[paste0("nscore", first_rep, "_c")]], 
                               y=.data[[paste0("nscore", second_rep, "_c")]] ))+
    stat_binhex()+
    theme_bw()+
    labs(x=paste0("Replicate ", first_rep), y=paste0("Replicate ", second_rep))+
    theme(  panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color='black'),
            axis.title  = element_text(size = 25),
            axis.text = element_text(size=18))+
    annotate("text", x = -5, y = 3, label = paste0("R=", round(R, 2)), size=8)+
    annotate("text", x = -5, y = 2, label = paste0("p=",format(p, digits = 2, scientific = T)), size=8)+
    scale_fill_gradient(high="grey30", low="grey90")
  p_corr
  
  my_list[[comb]]<-p_corr
  
}

p_all<-ggarrange(my_list$`1_2`, my_list$`1_3`, my_list$`2_3`, ncol=3, common.legend = T, legend = "right")
p_all

ggsave(p_all, file="p_corr_reps.pdf", width = 12, height = 3.5, path=path)







### correlation synthetic library and error-prone PCR library

synth_singles<-singles.df[,c("ID", "nscore_c", "sigma")]

singles_all<-full_join(synth_singles, ep_singles, by="ID")
colnames(singles_all)<-c("ID", "nscore_c_synth", "sigma_synth","nscore_c_ep", "sigma_ep" )


n=length(singles_all[!is.na(singles_all$nscore_c_ep) & !is.na(singles_all$nscore_c_synth),]$ID)
print(paste0("n= ", n))

corr<-cor.test(singles_all$nscore_c_synth,singles_all$nscore_c_ep)$estimate
p_value<-cor.test(singles_all$nscore_c_synth,singles_all$nscore_c_ep)$p.value


p_corr<-ggplot(singles_all, aes(x=nscore_c_synth, y=nscore_c_ep))+
  theme_bw()+
  stat_binhex()+
  labs(x="Synthetic library", y="epPCR library")+
  theme(  plot.title =element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'))+
  annotate("text", label=paste0("R=", round(corr, 2)), x=-4, y=2, size=6)+
  annotate("text", label=paste0("p=", format(p_value, digits = 2, scientific = T)), x=-4, y=1.2, size=6)+
  scale_x_continuous(limits=c(-6,2.5))+
  scale_y_continuous(limits=c(-6,2.5))+
  scale_fill_gradient(high="grey30", low="grey90")
p_corr

ggsave(p_corr, file="p_corr_libraries.pdf", width=3.8, height=3, path=path)








### correlation large and small scale

ind_var$mean_ind_score<-apply(ind_var[,c("rep1", "rep2", "rep3")], 1, mean)
ind_var$sd_ind_score<-apply(ind_var[,c("rep1", "rep2", "rep3")], 1, sd)

ind_var<-left_join(ind_var, INDEL.df[,c("ID","aa_seq", "nscore_c", "sigma")], by="aa_seq")
ind_var[ind_var$ID_small_scale=="AB42",c("nscore_c", "sigma")]<-2.2e-16


corr<-cor.test(ind_var$nscore_c, ind_var$mean_ind_score)$estimate
p_value<-cor.test(ind_var$nscore_c, ind_var$mean_ind_score)$p.value


p_small_large_scale<-ggplot(ind_var, aes(x=mean_ind_score, y=nscore_c))+
  geom_smooth(method = "lm", se=F, color="grey", linetype="dashed")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=nscore_c-1.96*sigma, ymax=nscore_c+1.96*sigma))+
  geom_errorbar(aes(xmin=mean_ind_score-1.96*sd_ind_score, xmax=mean_ind_score+1.96*sd_ind_score))+
  theme_bw()+
  labs(x="Small scale" ,y="Large scale")+
  theme(axis.title = element_text(size=20),
        axis.text= element_text(size=16),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color='black'))+
  annotate("text", x = 35, y = (-2), label=paste0("R=",round(corr, 2)), size=8)+
  annotate("text", x = 35, y = -2.8, label =paste0("p=",format(p_value, digits = 2, scientific = T)), size=8)
p_small_large_scale


ggsave(p_small_large_scale, file="p_small_large_scale.pdf", width = 4.5, height = 4, path=path)








### in vitro measurements of A-beta aggregation


kinetics<-left_join(kinetics, INDEL.df[,c("ID", "nscore_c", "sigma")], by="ID")
kinetics[kinetics$ID=="WT", c("nscore_c", "sigma")]<-2.2e-16
kinetics[kinetics$ID=="I-41-*",]$nscore_c<-INDEL.df[INDEL.df$aa_seq=="DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV",]$nscore_c
kinetics[kinetics$ID=="I-41-*",]$sigma<-INDEL.df[INDEL.df$aa_seq=="DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV",]$sigma


corr_vector=c()
pval_vector=c()
measurements<-c("k+kn_tht","k+k2_tht", "lambda_HDXMS", "kappa_HDXMS")

for (i in measurements){
  
  #weighted pearson correlation
  pearson<-wtd.cor(kinetics$nscore_c, log10(as.numeric(kinetics[[i]])), weight = kinetics$sigma^-2)
  
  corr<-pearson[1,1]
  corr_vector=c(corr_vector,corr)
  
  p.value<-pearson[1,4]
  pval_vector=c(pval_vector, as.numeric(p.value))
  
}

corr_text<-data.frame(label = corr_vector,variable=measurements)
pval_text<-data.frame(label = as.numeric(pval_vector),variable=measurements)



melt_kinetics<-melt(kinetics, id=c("nscore_c", "sigma", "ID"))

tht_measurements<-c("k+kn_tht", "k+k2_tht")

p_kinetics<-ggplot(melt_kinetics[melt_kinetics$variable %in% tht_measurements,], aes(x=as.numeric(value), y=as.numeric(nscore_c)))+
  facet_wrap(vars(factor(variable, levels=c( "k+kn_tht", "k+k2_tht"),
                         labels=c("primary (k+kn)","secondary (k+k2)"))),
             scales = "free", ncol=5)+
             #ncol=5)+
             
  geom_smooth(method='lm',linetype = 2, size=1,  color = "black", fill="grey75")+
  scale_x_continuous(trans="log10")+
  geom_point(size=3, color="black")+
  labs(y="Nucleation score")+
  theme_bw()+
  geom_text(data=corr_text[corr_text$variable %in% tht_measurements,], aes(label=paste0("R=",round(label, 2)), x=Inf, y=Inf),hjust=3,vjust=4, size=4, colour="black")+
  geom_text(data=pval_text[pval_text$variable %in% tht_measurements,], aes(label=paste0("p=",format(label, digits = 2, scientific = T)), x=Inf, y=Inf),hjust=3, vjust=6, size=4, colour="black")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.line = element_line(color="black", size=0.15),
        strip.background=element_rect( color=NA, fill=NA))
p_kinetics  

ggsave(p_kinetics, file="p_kinetics_NS.pdf", path=path, width = 6, height = 3)





lambda_kappa_measurements<-c("lambda_HDXMS", "kappa_HDXMS")

p_kinetics<-ggplot(melt_kinetics[melt_kinetics$variable %in% lambda_kappa_measurements,], aes(x=as.numeric(value), y=as.numeric(nscore_c)))+
 
   facet_wrap(vars(factor(variable, levels=c( "lambda_HDXMS", "kappa_HDXMS"),
                         labels=c("lambda (primary)","kappa (secondary)"))),
             scales = "free", ncol=5)+

  geom_smooth(method='lm',linetype = 2, size=1,  color = "black", fill="grey75")+
  scale_x_continuous(trans="log10")+
  scale_y_continuous(limits = c(-2, 4))+
  geom_point(size=3, color="black")+
  labs(y="Nucleation score")+
  theme_bw()+
  geom_text(data=corr_text[corr_text$variable %in% lambda_kappa_measurements,], aes(label=paste0("R=",round(label, 2)), x=Inf, y=Inf),hjust=3,vjust=4, size=4, colour="black")+
  geom_text(data=pval_text[pval_text$variable %in% lambda_kappa_measurements,], aes(label=paste0("p=",format(label, digits = 2, scientific = T)), x=Inf, y=Inf),hjust=3, vjust=6, size=4, colour="black")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.line = element_line(color="black", size=0.15),
        strip.background=element_rect( color=NA, fill=NA))
p_kinetics  

ggsave(p_kinetics, file="p_kinetics_NS_lambda_kappa.pdf", path=path, width = 6, height = 3)



p<-ggplot(kinetics, aes(x=log10(lambda_HDXMS),
                        y=log10(kappa_HDXMS)))+
  geom_point()+
  geom_text_repel(aes(label=ID))+
  theme_bw()+
  geom_abline()+
  scale_x_continuous(limits = c(-9,-1))+
  scale_y_continuous(limits = c(-9,-1))+
  labs(x="lambda (primary)_HDXMS (log10)", y="kappa (secondary)_HDXMS (log10)")+
  coord_fixed(ratio = 1)
p

ggsave(p, file="p_kappa_lambda.pdf", path=path,width = 3.5, height = 3.5)







### toxicity validation

toxicity_df$ID_condition_bio_rep<-paste(toxicity_df$ID, toxicity_df$condition, toxicity_df$bio_rep, sep=".")

my_df<-data.frame("ID_condition_bio_rep"=as.character(unique(toxicity_df$ID_condition_bio_rep)),
                  "mean_growth_rate"= unlist(lapply(unique(toxicity_df$ID_condition_bio_rep), function(x){
  mean(toxicity_df[toxicity_df$ID_condition_bio_rep==x,]$growth_rate, na.rm=T)
}))
)


my_df<-left_join(my_df, toxicity_df[,c("ID", "condition", "bio_rep", "ID_condition_bio_rep")], by="ID_condition_bio_rep")
my_df<-my_df[!duplicated(my_df),]

my_df$ID_condition<-paste(my_df$ID, my_df$condition, sep=".")



levels_test<-c(unique(my_df$ID_condition))
levels_test<-c("supN.cu", levels_test[levels_test!="supN.cu"])

duntest<-DunnettTest(my_df$mean_growth_rate, factor(my_df$ID_condition, levels=levels_test), na.rm=T)

duntest<-as.data.frame(duntest$supN.cu)
duntest$ID<-str_remove(row.names(duntest), "-supN.cu")

#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
duntest$sig<-""
tryCatch({ duntest[duntest$pval<0.05,]$sig<-"*" }, error=function(e){})
tryCatch({ duntest[duntest$pval<0.01,]$sig<-"**" }, error=function(e){})
tryCatch({ duntest[duntest$pval<0.001,]$sig<-"***" }, error=function(e){})


my_df<-left_join(my_df, duntest[,c("ID", "sig")], by=c("ID_condition"="ID"))



levels_id<-c("supN","AB42","a2v","k28q","g37l","g38v","AB11_42","e22d","e22g","l34i")
labels_id=c("SupN", "AB42","A2V","K28Q",  "G37L", "G38V","AB11-42", "E22d", "E22G", "L34I")


p_tox<-ggplot(my_df, aes(x=factor(ID, levels = levels_id, labels = labels_id), y=mean_growth_rate, 
                          color=factor(condition, levels=c("no_ind", "cu"), labels=c("no Cu2+", "Cu2+"))))+
  #geom_boxplot(size=0.2)+
  geom_point(aes(group=factor(condition, levels=c( "no_ind", "cu"), labels=c("no Cu2+", "Cu2+"))), 
             position=position_dodge(width = 0.5) , size=2)+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.line.y = element_line(),
        #axis.line.x = element_line(),
        axis.text.x=element_text(angle=45, vjust = 1, hjust=1))+
  scale_y_continuous(limits = c(0.1, 0.4), expand = c(0,0))+

  labs(x="", y="Growth rate", color="", group="")+
  scale_color_manual(values=c( "no Cu2+"="grey50", "Cu2+"="darkblue"))
p_tox


ggsave(p_tox, file="p_toxicity_validation.pdf",width = 6, height = 3, path=path)





