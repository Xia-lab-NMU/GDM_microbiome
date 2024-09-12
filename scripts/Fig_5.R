source_dir_root1<-"D:/Project/pes_GM_GDM/Data_review"
setwd(source_dir_root1)


# Figure 7. prediction performance----------------------------------------------------------------

library(openxlsx)
library(randomForest)
library(ggplot2)
library(ggsci)
library(openxlsx)
library(ggplot2)
library(DALEX)
library(iBreakDown)
library(nnet)
library(questionr)
library(dplyr)
library(caret)
library(ROCR)
library(pROC)
library(car)

#setwd("D:/project/pes_GM_GDM")
dir.create("Result/Figure7")

# load data ---------------------------------------------------------------


load("Data/microbialFeatures_standardized.RData")

# read meta data 
load("Data/metadata.RData")
load("Data/data_ers.RData")
rownames(metadata.meta)<-metadata.meta$SampleID
meta<-cbind(metadata.meta,metagTaxa.features.st)

library(dplyr)
dat_ers<-data_ers[,c("number","ers")]

meta<-left_join(meta,dat_ers,by="number")

meta<-cbind(meta,metagMod.features.st)

Species_name_clr<-paste0(Species_name,"_clr")

# 一、genus_level_selected1 ------------------------------------------------------------

dat_linear<-read.csv("Result/GM_GDM/metagenomic/Reg/linear_group.csv",header = T,sep = ",")

dat_sig<-subset(dat_linear,pvalue2<0.05)
dat_sig<-subset(dat_sig,genus%in%Species_name_clr)

taxonomy_selected<-as.character(dat_sig$genus)

#pesticide

dat_linear_pes<-read.table("Result/pes_GDM/linear/result_bino2.txt",header = T,sep = "\t")

dat_sig_pes<-subset(dat_linear_pes,p_3i<0.05)
dat_sig_pes<-subset(dat_sig_pes,i%in%pes_log_new)
dat_sig_pes<-subset(dat_sig_pes,j=="GDM")

pes_selected<-as.character(dat_sig_pes$i)


#pathway
dat_linear_path<-read.csv("Result/GM_GDM/metagenomic/Func/reg/linear_group.csv",header = T,sep = ",")
library(dplyr)

###read pathway names
path_func<-"Data/Sequencing/metagenomic/Data"

dat_metacyc<-read.csv(paste0(path_func,"/MetaCyc_pathway_filtered.csv"),header = T,sep = ",",row.names = 1)

dat_metacyc<-as.data.frame(t(dat_metacyc))

dat_metacyc<-dat_metacyc[,-c(1:2)]
metacyc_name<-colnames(dat_metacyc)
metacyc_name_new<-paste0("MetaCyc_",seq(1,length(metacyc_name)))

metacyc_n<-data.frame(metacyc_name_new,metacyc_name)

metacyc_n<- metacyc_n %>% separate(metacyc_name, c("pathway","description"), ":")


names(dat_metacyc)<-metacyc_name_new

dat_KO<-read.csv(paste0(path_func,"/KO_pathway_filtered.csv"),header = T,sep = ",",row.names = 1)


dat_KO<-as.data.frame(t(dat_KO))
dat_KO<-dat_KO[,-c(1:2)]
KO_name<-colnames(dat_KO)

dat_KEGG<-read.csv(paste0(path_func,"/KEGG_pathway_filtered.csv"),header = T,sep = ",",row.names = 1)

dat_KEGG<-as.data.frame(t(dat_KEGG))
#dat_KEGG<-dat_KEGG[,-c(1:2)]
#KEGG_name<-colnames(dat_KEGG)
KEGG_name<-colnames(dat_KEGG)
KEGG_name_new<-paste0("KEGG_",seq(1,length(KEGG_name)))

KEGG_n<-data.frame(KEGG_name_new,KEGG_name)
names(dat_KEGG)<-KEGG_name_new


all_name_func<-c(metacyc_name_new,KO_name,KEGG_name_new)


##filter pathway

dat_sig_path<-subset(dat_linear_path,pvalue2<0.05&outcome=="GDM")

KEGG_name_new_clr<-paste0(KEGG_name_new,"_clr")
dat_sig_path<-subset(dat_sig_path,genus%in%KEGG_name_new_clr)

pathway_selected<-as.character(dat_sig_path$genus)

cova<-c("pre_BMI","Age","Educational_level","weightgain","Parity","passsmokHis_1y")

dat_sig<-read.csv("Result/mediation/Result_combine2.csv",header = T,sep = ",")

taxonomy_selected<-subset(dat_sig,phenotype%in%c("GDM")&Group=="Bacteria metagenome")
taxonomy_selected<-unique(taxonomy_selected$taxon)

KEGG_selected<-subset(dat_sig,phenotype%in%c("GDM")&Group=="KEGG Module")
KEGG_selected<-unique(KEGG_selected$taxon)

MetaCyc_selected<-subset(dat_sig,phenotype%in%c("GDM")&Group=="MetaCyc Pathway")
MetaCyc_selected<-unique(MetaCyc_selected$taxon)

pes_selected<-c("ers","Dimethenamid_log","Dimethoate_log")

cova<-covariate


path_func<-"Data/Sequencing/metagenomic/Data"
dat_metacyc_name<-read.csv(paste0(path_func,"/metacyc_n.csv"),header = T,sep = ",",row.names = 1)
names(dat_metacyc_name)<-c("ID","Pathway","Description")
dat_KEGG_name<-read.csv(paste0(path_func,"/KEGG_n.csv"),header = T,sep = ",",row.names = 1)
dat_KEGG_name<-dat_KEGG_name[,c(1,3,2)]
names(dat_KEGG_name)<-c("ID","Pathway","Description")
dat_annotation<-rbind(dat_metacyc_name,dat_KEGG_name)

#1.taxonomy_selected -------------------------------------------------------


# 
# #1. genus alone ------------------------------------------------------------
dir_root<-"D:/project/pes_GM_GDM/Data_review/Result/Figure7/ML2"
dir.create(dir_root)
setwd(dir_root)
meta$GDM<-as.factor(meta$GDM)

meta$IGT<-as.factor(meta$IGT)

rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
  dir.create(i)
  dir.create(paste0(i,"/taxonomy_selected"))
  varia<-c(i,taxonomy_selected)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  ratio=j
  for (b in 1:m){
    
    Train <- createDataPartition(meta[,i], p=ratio, list=FALSE)
    training <- meta[ Train, varia]
    testing <- meta[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=stats::predict(rf,testing)
    df_b = table(observed = testing[,i], predict = pred_b)
    
    df_b<-as.data.frame(df_b)
    sensi_b<-df_b[4,3]/(df_b[2,3]+df_b[4,3]) 
    speci_b<-df_b[1,3]/(df_b[1,3]+df_b[3,3])
    
    pred_rf_b<-stats::predict(rf,testing,type = "prob")
    pred_2b_b<-prediction(pred_rf_b[,2],testing[,i])
    perf_b<-performance(pred_2b_b,"tpr","fpr")
    auc_b<-performance(pred_2b_b,"auc")@y.values
    auc_b<-as.numeric(auc_b)
    
    predictions2 <- stats::predict(rf, testing, type = "class")
    result_matrix<-confusionMatrix(predictions2, testing[,i])
    accura_b<-result_matrix$overall
    
    accura_b<-c(accura_b[1],accura_b[3],accura_b[4],auc_b,b)
    names(accura_b)<-c("Accuracy","Lower","Upper","AUC","Times")
    accura<-rbind(accura,accura_b)
    
    sensi<-result_matrix$byClass
    sensi<-as.data.frame(sensi)
    sensi<-as.data.frame(t(sensi))
    times<-rep(b,times=dim(sensi)[1])
    result_b<-data.frame(sensi,times)
    result_m<-rbind(result_m,result_b)
    
    pred_2b<-data.frame(pred_b,pred_rf_b)
    times<-rep(b,times=dim(pred_2b)[1])
    
    groundtruth<-testing[,i]
    pred_2b<-data.frame(pred_2b,groundtruth,times)
    pred<-rbind(pred,pred_2b)
    # rownames(pred)<-rownames(testing)
    #ROC曲线
    err_b<-rf$err.rate[500,1]
    result_b<-c(b,sensi_b,speci_b,auc_b,err_b)
    result_b<-as.data.frame(result_b) 
    result_b<-as.data.frame(t(result_b)) 
    result<-rbind(result,result_b) 
    
    sensi_value_b<-perf_b@y.values[[1]]
    speci_value_b<-perf_b@x.values[[1]]
    plot_data_b<-data.frame(sensi_value_b,speci_value_b)
    group<-rep(b,times=dim(plot_data_b)[1])
    plot_data_b<-data.frame(plot_data_b,group)
    plot_data<-rbind(plot_data,plot_data_b)
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("times","sensitivity","specifity","auc","error")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/taxonomy_selected/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/taxonomy_selected/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/taxonomy_selected/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/taxonomy_selected/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/taxonomy_selected/pred_",ratio,".csv"),row.names = T)
}



for (i in c("GDM")) {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}


# 2.genus+pesticide -------------------------------------------------------------

rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
  dir.create(i)
  dir.create(paste0(i,"/taxonomy_selected_pes"))
  varia<-c(i,taxonomy_selected,pes_selected)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  ratio=j
  for (b in 1:m){
    
    Train <- createDataPartition(meta[,i], p=ratio, list=FALSE)
    training <- meta[ Train, varia]
    testing <- meta[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=stats::predict(rf,testing)
    df_b = table(observed = testing[,i], predict = pred_b)
    
    df_b<-as.data.frame(df_b)
    sensi_b<-df_b[4,3]/(df_b[2,3]+df_b[4,3]) 
    speci_b<-df_b[1,3]/(df_b[1,3]+df_b[3,3])
    
    pred_rf_b<-stats::predict(rf,testing,type = "prob")
    pred_2b_b<-prediction(pred_rf_b[,2],testing[,i])
    perf_b<-performance(pred_2b_b,"tpr","fpr")
    auc_b<-performance(pred_2b_b,"auc")@y.values
    auc_b<-as.numeric(auc_b)
    
    predictions2 <- stats::predict(rf, testing, type = "class")
    result_matrix<-confusionMatrix(predictions2, testing[,i])
    accura_b<-result_matrix$overall
    
    accura_b<-c(accura_b[1],accura_b[3],accura_b[4],auc_b,b)
    names(accura_b)<-c("Accuracy","Lower","Upper","AUC","Times")
    accura<-rbind(accura,accura_b)
    
    sensi<-result_matrix$byClass
    sensi<-as.data.frame(sensi)
    sensi<-as.data.frame(t(sensi))
    times<-rep(b,times=dim(sensi)[1])
    result_b<-data.frame(sensi,times)
    result_m<-rbind(result_m,result_b)
    
    pred_2b<-data.frame(pred_b,pred_rf_b)
    times<-rep(b,times=dim(pred_2b)[1])
    
    groundtruth<-testing[,i]
    pred_2b<-data.frame(pred_2b,groundtruth,times)
    pred<-rbind(pred,pred_2b)
    # rownames(pred)<-rownames(testing)
    #ROC曲线
    err_b<-rf$err.rate[500,1]
    result_b<-c(b,sensi_b,speci_b,auc_b,err_b)
    result_b<-as.data.frame(result_b) 
    result_b<-as.data.frame(t(result_b)) 
    result<-rbind(result,result_b) 
    
    sensi_value_b<-perf_b@y.values[[1]]
    speci_value_b<-perf_b@x.values[[1]]
    plot_data_b<-data.frame(sensi_value_b,speci_value_b)
    group<-rep(b,times=dim(plot_data_b)[1])
    plot_data_b<-data.frame(plot_data_b,group)
    plot_data<-rbind(plot_data,plot_data_b)
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("times","sensitivity","specifity","auc","error")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/taxonomy_selected_pes/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/taxonomy_selected_pes/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/taxonomy_selected_pes/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/taxonomy_selected_pes/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/taxonomy_selected_pes/pred_",ratio,".csv"),row.names = T)
}



for (i in "GDM") {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}


# 3.genus+pesticide+pathway -------------------------------------------------------------

#meta<-merge(meta[,c("X",all_name_meta_clr)],data_total_func,by.x = "X",by.y = "X")
#meta$GDM<-as.factor(meta$GDM)
pathway_selected<-c(KEGG_selected,MetaCyc_selected)
rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
  dir.create(i)
  dir.create(paste0(i,"/taxonomy_selected_pes_path"))
  varia<-c(i,taxonomy_selected,pes_selected,pathway_selected)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  ratio=j
  for (b in 1:m){
    
    Train <- createDataPartition(meta[,i], p=ratio, list=FALSE)
    training <- meta[ Train, varia]
    testing <- meta[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=stats::predict(rf,testing)
    df_b = table(observed = testing[,i], predict = pred_b)
    
    df_b<-as.data.frame(df_b)
    sensi_b<-df_b[4,3]/(df_b[2,3]+df_b[4,3]) 
    speci_b<-df_b[1,3]/(df_b[1,3]+df_b[3,3])
    
    pred_rf_b<-stats::predict(rf,testing,type = "prob")
    pred_2b_b<-prediction(pred_rf_b[,2],testing[,i])
    perf_b<-performance(pred_2b_b,"tpr","fpr")
    auc_b<-performance(pred_2b_b,"auc")@y.values
    auc_b<-as.numeric(auc_b)
    
    predictions2 <- stats::predict(rf, testing, type = "class")
    result_matrix<-confusionMatrix(predictions2, testing[,i])
    accura_b<-result_matrix$overall
    
    accura_b<-c(accura_b[1],accura_b[3],accura_b[4],auc_b,b)
    names(accura_b)<-c("Accuracy","Lower","Upper","AUC","Times")
    accura<-rbind(accura,accura_b)
    
    sensi<-result_matrix$byClass
    sensi<-as.data.frame(sensi)
    sensi<-as.data.frame(t(sensi))
    times<-rep(b,times=dim(sensi)[1])
    result_b<-data.frame(sensi,times)
    result_m<-rbind(result_m,result_b)
    
    pred_2b<-data.frame(pred_b,pred_rf_b)
    times<-rep(b,times=dim(pred_2b)[1])
    
    groundtruth<-testing[,i]
    pred_2b<-data.frame(pred_2b,groundtruth,times)
    pred<-rbind(pred,pred_2b)
    # rownames(pred)<-rownames(testing)
    #ROC曲线
    err_b<-rf$err.rate[500,1]
    result_b<-c(b,sensi_b,speci_b,auc_b,err_b)
    result_b<-as.data.frame(result_b) 
    result_b<-as.data.frame(t(result_b)) 
    result<-rbind(result,result_b) 
    
    sensi_value_b<-perf_b@y.values[[1]]
    speci_value_b<-perf_b@x.values[[1]]
    plot_data_b<-data.frame(sensi_value_b,speci_value_b)
    group<-rep(b,times=dim(plot_data_b)[1])
    plot_data_b<-data.frame(plot_data_b,group)
    plot_data<-rbind(plot_data,plot_data_b)
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("times","sensitivity","specifity","auc","error")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/taxonomy_selected_pes_path/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/taxonomy_selected_pes_path/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/taxonomy_selected_pes_path/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/taxonomy_selected_pes_path/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/taxonomy_selected_pes_path/pred_",ratio,".csv"),row.names = T)
}



for (i in "GDM") {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}



# 4.genus+pesticide+pathway+cova -------------------------------------------------------------

rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
  dir.create(i)
  dir.create(paste0(i,"/taxonomy_selected_pes_path_cov"))
  varia<-c(i,taxonomy_selected,pes_selected,pathway_selected,cova)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  ratio=j
  for (b in 1:m){
    
    Train <- createDataPartition(meta[,i], p=ratio, list=FALSE)
    training <- meta[ Train, varia]
    testing <- meta[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=stats::predict(rf,testing)
    df_b = table(observed = testing[,i], predict = pred_b)
    
    df_b<-as.data.frame(df_b)
    sensi_b<-df_b[4,3]/(df_b[2,3]+df_b[4,3]) 
    speci_b<-df_b[1,3]/(df_b[1,3]+df_b[3,3])
    
    pred_rf_b<-stats::predict(rf,testing,type = "prob")
    pred_2b_b<-prediction(pred_rf_b[,2],testing[,i])
    perf_b<-performance(pred_2b_b,"tpr","fpr")
    auc_b<-performance(pred_2b_b,"auc")@y.values
    auc_b<-as.numeric(auc_b)
    
    predictions2 <- stats::predict(rf, testing, type = "class")
    result_matrix<-confusionMatrix(predictions2, testing[,i])
    accura_b<-result_matrix$overall
    
    accura_b<-c(accura_b[1],accura_b[3],accura_b[4],auc_b,b)
    names(accura_b)<-c("Accuracy","Lower","Upper","AUC","Times")
    accura<-rbind(accura,accura_b)
    
    sensi<-result_matrix$byClass
    sensi<-as.data.frame(sensi)
    sensi<-as.data.frame(t(sensi))
    times<-rep(b,times=dim(sensi)[1])
    result_b<-data.frame(sensi,times)
    result_m<-rbind(result_m,result_b)
    
    pred_2b<-data.frame(pred_b,pred_rf_b)
    times<-rep(b,times=dim(pred_2b)[1])
    
    groundtruth<-testing[,i]
    pred_2b<-data.frame(pred_2b,groundtruth,times)
    pred<-rbind(pred,pred_2b)
    # rownames(pred)<-rownames(testing)
    #ROC曲线
    err_b<-rf$err.rate[500,1]
    result_b<-c(b,sensi_b,speci_b,auc_b,err_b)
    result_b<-as.data.frame(result_b) 
    result_b<-as.data.frame(t(result_b)) 
    result<-rbind(result,result_b) 
    
    sensi_value_b<-perf_b@y.values[[1]]
    speci_value_b<-perf_b@x.values[[1]]
    plot_data_b<-data.frame(sensi_value_b,speci_value_b)
    group<-rep(b,times=dim(plot_data_b)[1])
    plot_data_b<-data.frame(plot_data_b,group)
    plot_data<-rbind(plot_data,plot_data_b)
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("times","sensitivity","specifity","auc","error")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/taxonomy_selected_pes_path_cov/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/taxonomy_selected_pes_path_cov/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/taxonomy_selected_pes_path_cov/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/taxonomy_selected_pes_path_cov/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/taxonomy_selected_pes_path_cov/pred_",ratio,".csv"),row.names = T)
}



for (i in "GDM") {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}

#五、 ROC_95%CI ---------------------------------------------------------------

dir_root<-"D:/project/pes_GM_GDM/Data_review/Result/Figure7/ML2"
setwd(dir_root)


# Figure7c ---------------------------------------------------------------



i<-"GDM"
j<-"taxonomy_selected_pes_path_cov"
k<-"0.8"

dir1<-paste0(dir_root,"/",i,"/",j)

result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")

result1<-result[order(result$auc,decreasing = T),]
imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")

library(pROC)
time<-result1[1,2]
set.seed(315)
pred1<-pred[pred$times==time,]
  # proc_obj_signature<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
  # sensi.ci<-ci.se(proc_obj_signature,specificities = seq(0, 1, 0.01))
  # plot(proc_obj_signature, type="shape", col="lightblue") # plot as a blue shape

ci.sp.obj <- ci.sp(rocobj, sensitivities=seq(0, 1, .01), boot.n=1000)
rocobj <- plot.roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
plot(ci.sp.obj, type="shape", col="lightblue")

ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_top.pdf"),width=5,height=5)
  #write.csv(data_ci,file = paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_top.csv"))

#   ##importance
# ggplot(imp,aes(x=variable,y=MeanDecreaseGini))+
#     geom_boxplot()+
#     theme_bw()+
#     # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
#     coord_flip()
# ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_varimp.pdf"),width=10,height=10)


#2. varImportance -----------------------------------------------------------


# plot_imp<-function(i,j,k){
#   
#   dir1<-paste0(dir_root,"/",i,"/",j)
#   
#   #result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
#   #plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
#   #pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
#   
#   imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")
#   
#   imp<-left_join(imp,dat_annotation,by=c("variable"="ID"))
#   
#   imp$variable_new<-ifelse(imp$variable%in%dat_annotation$ID,imp$Description,imp$variable)
#   
#   
#   imp$variable_new<-as.character(imp$variable_new)
#   t<-tapply(imp$MeanDecreaseGini,imp$variable_new,mean)##get mean value of imp
#   t<-t[order(t,decreasing = T)]
#   var_name<-names(t)
#   
#   # mean_values<-aggregate(MeanDecreaseGini~variable_new,data = imp,FUN=mean)
#   # mean_values$Above_5<-mean_values$MeanDecreaseGini>5
#   ggplot(imp,aes(x=variable_new,y=MeanDecreaseGini))+
#     geom_boxplot(fill="#95C6C6")+
#     theme_bw()+
#     scale_x_discrete(limits=rev(var_name))+
#     # scale_fill_manual(values = ifelse(mean_values$Above_5, "#E98DAF", "#95C6C6")) +  # 根据均值是否大于5填充颜色
#     # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
#     guides(fill="none")+
#     coord_flip()
#   
#   
#   ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_varimp.pdf"),width=10,height=10)
#   
# }
# 
# outome<-"GDM"
# var2<-c("taxonomy_selected","taxonomy_selected_pes","taxonomy_selected_pes_path","taxonomy_selected_pes_path_cov")
# ratio<-c(0.8)
# for (i in outome){
#   for (j in var2) {
#     for (k in ratio) {
#       plot_imp(i,j,k)
#     }
#   }
# }


# Figure7d ----------------------------------------------------------------
setwd(source_dir_root1)
i<-"GDM"
j<-"taxonomy_selected"
k<-"0.8"

#plot_imp<-function(i,j,k){
  
  dir1<-paste0(dir_root,"/",i,"/",j)
  
  #result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
  #plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
  #pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
  
  imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")
  
  imp<-left_join(imp,dat_annotation,by=c("variable"="ID"))
  
  imp$variable_new<-ifelse(imp$variable%in%dat_annotation$ID,imp$Description,imp$variable)
  
  
  imp$variable_new<-as.character(imp$variable_new)
  t<-tapply(imp$MeanDecreaseGini,imp$variable_new,mean)##get mean value of imp
  t<-t[order(t,decreasing = T)]
  var_name<-names(t)
  
  # mean_values<-aggregate(MeanDecreaseGini~variable_new,data = imp,FUN=mean)
  # mean_values$Above_5<-mean_values$MeanDecreaseGini>5
  ggplot(imp,aes(x=variable_new,y=MeanDecreaseGini))+
    geom_boxplot(fill="#95C6C6")+
    theme_bw()+
    scale_x_discrete(limits=rev(var_name))+
    # scale_fill_manual(values = ifelse(mean_values$Above_5, "#E98DAF", "#95C6C6")) +  # 根据均值是否大于5填充颜色
    # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    guides(fill="none")+
    coord_flip()
  
  ggsave("Result/Figure7/Figure7d.pdf",width = 10,height = 10)
  #ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_varimp.pdf"),width=10,height=10)
  
#}



# FIgure7a-----------------------------------------------------

# i<-"Hypercholesterolemia_third"
# j<-"Family_test"

outome<-"GDM"
var2<-c("taxonomy_selected","taxonomy_selected_pes","taxonomy_selected_pes_path","taxonomy_selected_pes_path_cov")

imp_total<-data.frame()
plot_signature_top<-data.frame()
plot_signature_median<-data.frame()
auc_median<-data.frame()
auc_top<-data.frame()
data_ci_median<-data.frame()
data_ci_top<-data.frame()
for (i in outome) {
  for (j in var2) {
    for (k in c(0.8)) {
      dir1<-paste0(dir_root,"/",i,"/",j)
      
      result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
      plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
      pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
      #data_ci<-read.csv(paste0(dir1,"/",i,"_",j,"_",k,"_roc_ci_top.csv"),header = T,sep = ",")
      ####imp
      imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")
      Group<-rep(i,times=nrow(imp))
      variab<-rep(j,times=nrow(imp))
      rati<-rep(k,times=nrow(imp))
      imp_i<-data.frame(imp,Group,variab,rati)
      imp_total<-rbind(imp_total,imp_i)
      
      result1<-result[order(result$auc,decreasing = T),]
      
      library(pROC)
      time<-result1[50,2]
      set.seed(315)
      pred1<-pred[pred$times==time,]
      proc_obj_signature<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
      sensi.ci<-ci.se(proc_obj_signature,specificities = seq(0, 1, .01))
      
      auc<-auc(proc_obj_signature)[1]
      auc_low<-ci(proc_obj_signature,of="auc")[1]
      auc_high<-ci(proc_obj_signature,of="auc")[3]
      auc_signature<-paste(j," AUC:",round(auc,digits = 3),"(",
                           round(auc_low,digits = 3),",",round(auc_high,digits = 3),")",sep = "")
      
      auc_median1<-c(auc,auc_low,auc_high,i,j,k)
      auc_median1<-as.data.frame(auc_median1)
      auc_median1<-as.data.frame(t(auc_median1))
      auc_median<-rbind(auc_median,auc_median1)
      
      data_ci<-sensi.ci[1:101,1:3]
      data_ci<-as.data.frame(data_ci)
      y1<-data_ci[,1]
      y2<-data_ci[,3]
      x<-seq(0,1,0.01)
      data_ci_signature<-data.frame(y1,y2,x)
      Group<-rep(i,times=nrow(data_ci))
      variab<-rep(j,times=nrow(data_ci))
      rati<-rep(k,times=nrow(data_ci))
      data_ci<-data.frame(data_ci,Group,variab,rati)
      data_ci_median<-rbind(data_ci_median,data_ci)
      
      
      sensitivity<-proc_obj_signature$sensitivities
      specificity<-proc_obj_signature$specificities
      Group<-rep(i,times=length(sensitivity))
      variab<-rep(j,times=length(sensitivity))
      rati<-rep(k,times=length(sensitivity))
      plot_signature_i<-data.frame(sensitivity,specificity,Group,variab,rati)
      plot_signature_median<-rbind(plot_signature_median,plot_signature_i)
      #dev.off()
      ###top value
      library(pROC)
      time<-result1[1,2]
      set.seed(315)
      pred1<-pred[pred$times==time,]
      proc_obj_signature<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
      sensi.ci<-ci.se(proc_obj_signature,specificities = seq(0, 1, .01))
      
      auc<-auc(proc_obj_signature)[1]
      auc_low<-ci(proc_obj_signature,of="auc")[1]
      auc_high<-ci(proc_obj_signature,of="auc")[3]
      auc_signature<-paste(j," AUC:",round(auc,digits = 3),"(",
                           round(auc_low,digits = 3),",",round(auc_high,digits = 3),")",sep = "")
      
      auc_high1<-c(auc,auc_low,auc_high,i,j,k)
      auc_high1<-as.data.frame(auc_high1)
      auc_high1<-as.data.frame(t(auc_high1))
      auc_top<-rbind(auc_top,auc_high1)
      
      data_ci<-sensi.ci[1:101,1:3]
      data_ci<-as.data.frame(data_ci)
      y1<-data_ci[,1]
      y2<-data_ci[,3]
      x<-seq(0,1,0.01)
      data_ci_signature<-data.frame(y1,y2,x)
      Group<-rep(i,times=nrow(data_ci))
      variab<-rep(j,times=nrow(data_ci))
      rati<-rep(k,times=nrow(data_ci))
      data_ci<-data.frame(data_ci,Group,variab,rati)
      data_ci_top<-rbind(data_ci_top,data_ci)
      
      
      sensitivity<-proc_obj_signature$sensitivities
      specificity<-proc_obj_signature$specificities
      Group<-rep(i,times=length(sensitivity))
      variab<-rep(j,times=length(sensitivity))
      rati<-rep(k,times=length(sensitivity))
      plot_signature_1i<-data.frame(sensitivity,specificity,Group,variab,rati)
      plot_signature_top<-rbind(plot_signature_top,plot_signature_1i)
      
      
    }
  }
}

setwd(paste0(source_dir_root1,"/Result/Figure7"))
write.csv(plot_signature_top,file = "plot_signature_top.csv")
write.csv(plot_signature_median,file = "plot_signature_median.csv")
write.csv(imp_total,file = "imp_total.csv")
write.csv(data_ci_median,file = "data_ci_median.csv")
write.csv(data_ci_top,file = "data_ci_top.csv")
write.csv(auc_top,file = "auc_top.csv")
write.csv(auc_median,file = "auc_median.csv")

imp_total<-read.csv("imp_total.csv",header = T,sep = ",")
plot_signature_median<-read.csv("plot_signature_median.csv",header = T,sep = ",")
plot_signature_top<-read.csv("plot_signature_top.csv",header = T,sep = ",")

# data_ci_median<-read.csv("data_ci_median.csv",header = T,sep = ",")

#pdf(file="Hypercholesterolemia/family_test/roc_ci_2.pdf",width = 5,height = 5)

require(ggplot2)
library(ggsci)
for (i in outome) {
  for (j in c(0.8)) {
    plot_data<-subset(plot_signature_top,Group==i&rati==j)
    
    plot_data<-subset(plot_signature_top,Group==i&rati==j)
    p<-ggplot(plot_data, aes(x = 1-specificity, y=sensitivity)) +
      geom_path(aes(color = variab), size=1) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                   colour='grey', linetype = 'dotdash') +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.justification=c(1, 0), legend.position=c(.95, .05),
            legend.title=element_blank(), 
            legend.background = element_rect(fill=NULL, size=0.5, 
                                             linetype="solid", colour ="black"))+
      labs(title = i)+
      scale_color_npg()
   # ggsave(p,filename = paste0(i,"_",j,".pdf"),width = 8,height = 8)
    ggsave(p,filename = "Figure7a.pdf",width = 8,height = 8)
    
  }
  
}


# 
# ggplot(imp,aes(x=variable,y=MeanDecreaseGini))+
#   geom_boxplot()+
#   theme_bw()+
#   scale_x_discrete(limits=rev(var_name))+
#   # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
#   coord_flip()
# ggsave(paste0(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_varimp.pdf")),width=10,height=10)
# 





#Figure7b --------------------------------------------------------------------

dir_root<-"D:/project/pes_GM_GDM/Data_review/Result/Figure7/ML2"
setwd(dir_root)
outome<-"GDM"
var2<-c("taxonomy_selected","taxonomy_selected_pes","taxonomy_selected_pes_path","taxonomy_selected_pes_path_cov")
ratio<-c(0.8)

i<-outome[1]
j<-var2[1]


data_full<-data.frame()
data_full_total<-data.frame()
for (i in outome) {
  for (j in var2) {
    dir1<-paste0(dir_root,"/",i,"/",j)
    files2<-list.files(dir1,pattern = "accura*.")
    setwd(dir1)
    # data_full<-data.frame()
    for (k in files2){
      dat<-read.csv(k,header = T,sep = ",")
      group<-rep(j,times=dim(dat)[1])
      group2<-rep(i,times=dim(dat)[1])
      group3<-rep(k,times=dim(dat)[1])
      dat<-data.frame(dat,group,group2,group3)
      data_full<-rbind(data_full,dat)
    }
    # category<-rep(i,times=nrow(data_full))
    # data_full2<-data.frame(data_full,category)
    data_full_total<-rbind(data_full_total,data_full)
    #write.csv(data_full,file = "data_full.csv",row.names = F)
  }
  
}

data_full_2<-subset(data_full_total,group3=="accura_0.8.csv")
library(ggpubr)
setwd("D:/project/pes_GM_GDM/Data_review/Result/Figure7")
#dir.create("combine")
#setwd("D:/project/pes_GM_GDM/Data_review/Result/ML2/combine")
# for (i in var2) {
#   plot_dat<-subset(data_full_2,group==i)
#   p1<-ggboxplot(plot_dat,x="group2",y="AUC",#add="iitter",
#                 color = "group2")+
#     # stat_compare_means(method = "kruskal.test")+
#     theme_classic()+
#     theme(legend.position = "none")+
#     theme(panel.border = element_rect(color = "black",fill = NA,size = 1))+
#     
#     scale_color_npg()+
#     labs(x="",y="AUC")
#   #  facet_wrap(~group_2,nrow = 2)
#   p2<-ggboxplot(plot_dat,x="group2",y="Accuracy",#add="iitter",
#                 color = "group2")+
#     # stat_compare_means(method = "kruskal.test")+
#     theme_classic()+
#     theme(legend.position = "none")+
#     theme(panel.border = element_rect(color = "black",fill = NA,size = 1))+
#     
#     scale_color_npg()+
#     # scale_fill_manual(values = cbPalette)+
#     
#     labs(x="",y="Accuracy")
#   ggsave(p1,filename =paste(i,"_AUC.pdf"),width =6,height = 4 )
#   ggsave(p2,filename =paste(i,"_Accuracy.pdf"),width = 6,height = 4 )
#   
# }




p3<-ggboxplot(data_full_2,x="group",y="AUC",#add="iitter",
              color = "group")+
  # stat_compare_means(method = "kruskal.test")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black",fill = NA,size = 1))+
  
  scale_color_npg()+
  labs(x="",y="AUC")
ggsave(p3,filename ="Figure7b1.pdf",width = 8,height = 4 )


p4<-ggboxplot(data_full_2,x="group",y="AUC",#add="iitter",
              color = "group")+
  # stat_compare_means(method = "kruskal.test")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black",fill = NA,size = 1))+
  
  scale_color_npg()+
  labs(x="",y="Accuracy")
ggsave(p4,filename ="Figure7b2.pdf",width = 8,height = 4 )

# 
# ggsave(p3,filename ="AUC_combine.png",width = 8,height = 4 )
# 
# ggsave(p4,filename ="Accuracy_combine.png",width = 8,height = 4 )

# Figure S9. Performance of random forest classifier predicting glucose --------


#setwd("D:/BaiduNetdiskDownload/amplicon/33MachineLearning/RF_classification")
library(openxlsx)
library(ggplot2)

setwd("D:/project/pes_GM_GDM/Data_review")
dir.create("Result/FigureS9")

dir.create("Result/FigureS9/ML2")

# load data ---------------------------------------------------------------


load("Data/microbialFeatures_standardized.RData")

# read meta data 
load("Data/metadata.RData")
load("Data/data_ers.RData")
rownames(metadata.meta)<-metadata.meta$SampleID
meta<-cbind(metadata.meta,metagTaxa.features.st)

library(dplyr)
dat_ers<-data_ers[,c("number","ers")]

meta<-left_join(meta,dat_ers,by="number")

meta<-cbind(meta,metagMod.features.st)


#二、 taxonomy selected -------------------------------------------------------------

#随机森林分类
library(randomForest)
library(caret)
library(ggplot2)
library(bootstrap)
library(ROCR)
library(pROC)
set.seed(315)

dir_root<-"D:/project/pes_GM_GDM/Data_review/Result/FigureS9/ML2"
setwd(dir_root)

# 
# #1. genus alone ------------------------------------------------------------

rf_classify<-function(i){
  # i<-"CHOL_third"
  dir.create(i)
  dir.create(paste0(i,"/taxonomy_selected"))
  varia<-c(i,taxonomy_selected)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  result<-data.frame()
  b=1
  df<-data.frame()
  for (b in 1:m){
    data_meta<-na.omit(meta[,varia])
    Train <- createDataPartition(data_meta[,i], p=0.8, list=FALSE)
    training <- na.omit(data_meta[ Train, varia])
    testing <- na.omit(data_meta[ -Train, varia])
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=stats::predict(rf,testing,type="response")
    df_b = data.frame(observed = testing[,i], predict = pred_b)
    times2<-rep(b,times=nrow(df_b))
    
    
    cor = cor.test(df_b[,1], df_b[,2], method = "spearman")
    df2 = df_b
    colnames(df2) = c("x", "y")
    m = lm(y ~ x, df2)
    rho<-round(cor$estimate, digits = 3)
    pvalue<-round(cor$p.value, digits = 3)
    R2<-round(summary(m)$r.squared, digits = 3)
    result_i<-c(rho,pvalue,R2)
    result<-rbind(result,result_i)
    
    df_b<-data.frame(df_b,times2)
    df<-rbind(df,df_b)
    # p = ggplot(df, aes(predict, observed)) +
    #   geom_point() +
    #   geom_smooth(method = "lm") +
    #   labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    #   theme_bw()
    # p
    # ggsave(paste("test_predict.pdf", sep=""), p, width = 4, height = 2.5)
    # 
    # imp= as.data.frame(rf$importance)
    # imp = imp[order(imp[,1],decreasing = T),]
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("rho","pvalue","R2")
  # names(result)<-c("times","sensitivity","specifity","auc")
  # names(plot_data)<-c("sensitivity","specificity","Group")
  # names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  # 
  # write.csv(accura,file=paste0(i,"/genus_selected2/accura.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/taxonomy_selected/imp.csv"),row.names = T)
  write.csv(df,file=paste0(i,"/taxonomy_selected/df.csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/taxonomy_selected/result.csv"),row.names = T)
  # write.csv(plot_data,file=paste0(i,"/genus_selected2/plot_data.csv"),row.names = T)
  # write.csv(pred,file=paste0(i,"/genus_selected2/pred.csv"),row.names = T)
}


for (i in c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")) {
  rf_classify(i)
}

# #2.  genus+pesticide ------------------------------------------------------------

rf_classify<-function(i){
  # i<-"CHOL_third"
  dir.create(i)
  dir.create(paste0(i,"/taxonomy_selected_pes"))
  varia<-c(i,taxonomy_selected,pes_selected)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  result<-data.frame()
  b=1
  df<-data.frame()
  for (b in 1:m){
    data_meta<-na.omit(meta[,varia])
    Train <- createDataPartition(data_meta[,i], p=0.8, list=FALSE)
    training <- na.omit(data_meta[ Train, varia])
    testing <- na.omit(data_meta[ -Train, varia])
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=stats::predict(rf,testing,type="response")
    df_b = data.frame(observed = testing[,i], predict = pred_b)
    times2<-rep(b,times=nrow(df_b))
    
    
    cor = cor.test(df_b[,1], df_b[,2], method = "spearman")
    df2 = df_b
    colnames(df2) = c("x", "y")
    m = lm(y ~ x, df2)
    rho<-round(cor$estimate, digits = 3)
    pvalue<-round(cor$p.value, digits = 3)
    R2<-round(summary(m)$r.squared, digits = 3)
    result_i<-c(rho,pvalue,R2)
    result<-rbind(result,result_i)
    
    df_b<-data.frame(df_b,times2)
    df<-rbind(df,df_b)
    # p = ggplot(df, aes(predict, observed)) +
    #   geom_point() +
    #   geom_smooth(method = "lm") +
    #   labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    #   theme_bw()
    # p
    # ggsave(paste("test_predict.pdf", sep=""), p, width = 4, height = 2.5)
    # 
    # imp= as.data.frame(rf$importance)
    # imp = imp[order(imp[,1],decreasing = T),]
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("rho","pvalue","R2")
  # names(result)<-c("times","sensitivity","specifity","auc")
  # names(plot_data)<-c("sensitivity","specificity","Group")
  # names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  # 
  # write.csv(accura,file=paste0(i,"/genus_selected2/accura.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/taxonomy_selected_pes/imp.csv"),row.names = T)
  write.csv(df,file=paste0(i,"/taxonomy_selected_pes/df.csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/taxonomy_selected_pes/result.csv"),row.names = T)
  # write.csv(plot_data,file=paste0(i,"/genus_selected2/plot_data.csv"),row.names = T)
  # write.csv(pred,file=paste0(i,"/genus_selected2/pred.csv"),row.names = T)
}


for (i in c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")) {
  rf_classify(i)
}

# #3. genus+pesticide+path ------------------------------------------------------------

pathway_selected<-c(KEGG_selected,MetaCyc_selected)
rf_classify<-function(i){
  # i<-"CHOL_third"
  dir.create(i)
  dir.create(paste0(i,"/taxonomy_selected_pes_path"))
  varia<-c(i,taxonomy_selected,pes_selected,pathway_selected)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  result<-data.frame()
  b=1
  df<-data.frame()
  for (b in 1:m){
    data_meta<-na.omit(meta[,varia])
    Train <- createDataPartition(data_meta[,i], p=0.8, list=FALSE)
    training <- na.omit(data_meta[ Train, varia])
    testing <- na.omit(data_meta[ -Train, varia])
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=stats::predict(rf,testing,type="response")
    df_b = data.frame(observed = testing[,i], predict = pred_b)
    times2<-rep(b,times=nrow(df_b))
    
    
    cor = cor.test(df_b[,1], df_b[,2], method = "spearman")
    df2 = df_b
    colnames(df2) = c("x", "y")
    m = lm(y ~ x, df2)
    rho<-round(cor$estimate, digits = 3)
    pvalue<-round(cor$p.value, digits = 3)
    R2<-round(summary(m)$r.squared, digits = 3)
    result_i<-c(rho,pvalue,R2)
    result<-rbind(result,result_i)
    
    df_b<-data.frame(df_b,times2)
    df<-rbind(df,df_b)
    # p = ggplot(df, aes(predict, observed)) +
    #   geom_point() +
    #   geom_smooth(method = "lm") +
    #   labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    #   theme_bw()
    # p
    # ggsave(paste("test_predict.pdf", sep=""), p, width = 4, height = 2.5)
    # 
    # imp= as.data.frame(rf$importance)
    # imp = imp[order(imp[,1],decreasing = T),]
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("rho","pvalue","R2")
  # names(result)<-c("times","sensitivity","specifity","auc")
  # names(plot_data)<-c("sensitivity","specificity","Group")
  # names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  # 
  # write.csv(accura,file=paste0(i,"/genus_selected2/accura.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/taxonomy_selected_pes_path/imp.csv"),row.names = T)
  write.csv(df,file=paste0(i,"/taxonomy_selected_pes_path/df.csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/taxonomy_selected_pes_path/result.csv"),row.names = T)
  # write.csv(plot_data,file=paste0(i,"/genus_selected2/plot_data.csv"),row.names = T)
  # write.csv(pred,file=paste0(i,"/genus_selected2/pred.csv"),row.names = T)
}


for (i in c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")) {
  rf_classify(i)
}

# #4. genus+pesticide+pathway+cova ------------------------------------------------------------

rf_classify<-function(i){
  # i<-"CHOL_third"
  dir.create(i)
  dir.create(paste0(i,"/taxonomy_selected_pes_path_cov"))
  varia<-c(i,taxonomy_selected,pes_selected,pathway_selected,cova)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  result<-data.frame()
  b=1
  df<-data.frame()
  for (b in 1:m){
    data_meta<-na.omit(meta[,varia])
    Train <- createDataPartition(data_meta[,i], p=0.8, list=FALSE)
    training <- na.omit(data_meta[ Train, varia])
    testing <- na.omit(data_meta[ -Train, varia])
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=stats::predict(rf,testing,type="response")
    df_b = data.frame(observed = testing[,i], predict = pred_b)
    times2<-rep(b,times=nrow(df_b))
    
    
    cor = cor.test(df_b[,1], df_b[,2], method = "spearman")
    df2 = df_b
    colnames(df2) = c("x", "y")
    m = lm(y ~ x, df2)
    rho<-round(cor$estimate, digits = 3)
    pvalue<-round(cor$p.value, digits = 3)
    R2<-round(summary(m)$r.squared, digits = 3)
    result_i<-c(rho,pvalue,R2)
    result<-rbind(result,result_i)
    
    df_b<-data.frame(df_b,times2)
    df<-rbind(df,df_b)
    # p = ggplot(df, aes(predict, observed)) +
    #   geom_point() +
    #   geom_smooth(method = "lm") +
    #   labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    #   theme_bw()
    # p
    # ggsave(paste("test_predict.pdf", sep=""), p, width = 4, height = 2.5)
    # 
    # imp= as.data.frame(rf$importance)
    # imp = imp[order(imp[,1],decreasing = T),]
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("rho","pvalue","R2")
  # names(result)<-c("times","sensitivity","specifity","auc")
  # names(plot_data)<-c("sensitivity","specificity","Group")
  # names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  # 
  # write.csv(accura,file=paste0(i,"/genus_selected2/accura.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/taxonomy_selected_pes_path_cov/imp.csv"),row.names = T)
  write.csv(df,file=paste0(i,"/taxonomy_selected_pes_path_cov/df.csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/taxonomy_selected_pes_path_cov/result.csv"),row.names = T)
  # write.csv(plot_data,file=paste0(i,"/genus_selected2/plot_data.csv"),row.names = T)
  # write.csv(pred,file=paste0(i,"/genus_selected2/pred.csv"),row.names = T)
}


for (i in c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")) {
  rf_classify(i)
}




#1. Figure S9----------------------------------------------------

library(pROC)
dir_root<-"D:/project/pes_GM_GDM/Data_review/Result/FigureS9/ML2"
setwd(dir_root)
dir.create("combine")

setwd("D:/project/pes_GM_GDM/Data_review/Result/FigureS9")
plot_result<-function(i,j){
  dir1<-paste0(dir_root,"/",i,"/",j)
  
  result<-read.csv(paste0(dir1,"/result.csv"),header = T,sep = ",")
  # plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
  pred<-read.csv(paste0(dir1,"/df.csv"),header = T,sep = ",")
  
  result1<-result[order(result$rho,decreasing = T),]
  imp<-read.csv(paste0(dir1,"/imp.csv"),header = T,sep = ",")
  
   time<-result1[1,1]
  set.seed(315)
  pred1<-pred[pred$times2==time,]
  cor = cor.test(pred1$predict, pred1$observed, method = "spearman")
  df2 = pred1[,c("predict","observed")]
  colnames(df2) = c("x", "y")
  m = lm(y ~ x, df2)
  p = ggplot(pred1, aes(predict, observed)) +
    geom_point() +
    geom_smooth(method = "lm",color="darkblue") +
    labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    theme_bw()
  p
  #ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_reg_top.pdf"),width=5,height=5)  
  dir.create(i)
  ggsave(paste0(i,"/",j,"_reg.pdf"),width=5,height=5)  
  ##importance

  
  
}

#i<-"TG_third"
#j<-"Family_alone"

outome<-c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")
var2<-c("taxonomy_selected","taxonomy_selected_pes","taxonomy_selected_pes_path","taxonomy_selected_pes_path_cov")

for (i in outome) {
  for (j in var2) {
    plot_result(i,j)
  }
}

# # 2.combine reg linear plot & varimportance plot-----------------------------------------------------
# 
# #i<-"TG_third"
# #j<-"Family_alone"
# 
# outome<-c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")
# var2<-c("taxonomy_selected","taxonomy_selected_pes","taxonomy_selected_pes_path","taxonomy_selected_pes_path_cov")
# 
# imp_total<-data.frame()
# pred_top<-data.frame()
# pred_median<-data.frame()
# 
# for (i in outome) {
#   for (j in var2) {
#     dir1<-paste0(dir_root,"/",i,"/",j)
#     
#     result<-read.csv(paste0(dir1,"/result.csv"),header = T,sep = ",")
#     # plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
#     pred<-read.csv(paste0(dir1,"/df.csv"),header = T,sep = ",")
#     
#     result1<-result[order(result$rho,decreasing = T),]
#     ####imp
#     imp<-read.csv(paste0(dir1,"/imp.csv"),header = T,sep = ",")
#     Group<-rep(i,times=nrow(imp))
#     variab<-rep(j,times=nrow(imp))
#     imp_i<-data.frame(imp,Group,variab)
#     imp_total<-rbind(imp_total,imp_i)
#     time<-result1[50,1]
#     set.seed(315)
#     pred1<-pred[pred$times2==time,]
#     
#     Group<-rep(i,times=nrow(pred1))
#     variab<-rep(j,times=nrow(pred1))
#     imp_i<-data.frame(pred1,Group,variab)
#     
#     pred_median<-rbind(pred_median,pred1)
#     
#     time<-result1[1,1]
#     set.seed(315)
#     pred1<-pred[pred$times2==time,]
#     Group<-rep(i,times=nrow(pred1))
#     variab<-rep(j,times=nrow(pred1))
#     imp_i<-data.frame(pred1,Group,variab)
#     pred_top<-rbind(pred_top,pred1)
#     
#     
#     
#   }
# }
# 
# 
# write.csv(pred_top,file = "combine//reg_pred_top.csv")
# write.csv(pred_median,file = "combine//reg_pred_median.csv")
# write.csv(imp_total,file = "combine//reg_imp_total.csv")
# 
# library(ggsci)
# 
# for (j in var2) {
#   imp_plot<-subset(imp_total,variab==j)
#   ggplot(imp_plot,aes(x=variable,y=X.IncMSE,fill=Group))+
#     geom_boxplot()+
#     theme_bw()+
#     scale_fill_nejm()+
#     # facet_wrap(~variab,ncol = 4)+
#     #  scale_x_discrete(limits=rev(var_name))+
#     # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
#     coord_flip()
#   ggsave(paste0("combine/varimp_reg_",j,".pdf"),width=10,height=12)
#   
# }




