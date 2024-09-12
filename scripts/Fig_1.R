
source_dir_root1<-"D:/Project/pes_GM_GDM/Data_review"
setwd(source_dir_root1)
package_list <- c("readxl","pastecs","carData","car",
                  "compareGroups","ggsci",
                  "tidyr","reshape2",
                  "dplyr","ggplot2","bkmr","gWQS",#"DSA",
                  "caret","glmnet","parallel","doParallel",
                  "rexposome","qgcomp","openxlsx","knitr",
                  "ggpubr","psych","missForest")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    # install.packages(p, repos=site)
    install.packages(p,dependencies = TRUE)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


path_metadata<-"Data/"
path_16s<-"Data/Sequencing/16S/Data/"
path_meta<-"Data/Sequencing/metagenomic/Data/"


data_combine<-read.xlsx("Data/data_full2.xlsx",sheet = 1)

data_combine$bloodsample[data_combine$blood.sample==""]<-0
data_combine$bloodsample[data_combine$blood.sample!=0]<-1

data_combine$PTB[data_combine$Delivery_week<37]<-1
data_combine$PTB[data_combine$Delivery_week>=37]<-0

data_combine$LBW[data_combine$Birthweight<2500]<-1
data_combine$LBW[data_combine$Birthweight>=2500]<-0

data_combine$Birthweight<-as.numeric(data_combine$Birthweight)
data_combine$OGTT0_24<-as.numeric(data_combine$OGTT0_24)

data_combine$OGTT1_24<-as.numeric(data_combine$OGTT1_24)
data_combine$OGTT2_24<-as.numeric(data_combine$OGTT2_24)
data_combine$HbAlc_24<-as.numeric(data_combine$HbAlc_24)

# data_combine$GDM[data_combine$OGTT0_24>=5.1]<-1
# data_combine$GDM[data_combine$OGTT1_24>=10]<-1
# data_combine$GDM[data_combine$OGTT2_24>=8.5]<-1
# 
# data_combine$GDM[data_combine$OGTT0_24<5.1&data_combine$OGTT1_24<10&data_combine$OGTT2_24<8.5]<-0
# #data_combine$GDM[data_combine$FBG_32>=5.1]<-1
data_combine$IGT<-ifelse(data_combine$OGTT0_24<7 & data_combine$OGTT2_24>=7.8 & data_combine$OGTT2_24<=11.1,1,0)

# 一、Inclusion criterion ---------------------------------------------------


# Fig 1b. ------------------------------------------------------------------

#having pesticides detected at 24week
pes_name<-colnames(data_combine[,175:203])
pes<-pes_name

#data2<-subset(data_combine,FBG_24<5.1)#1491
data2<-subset(data_combine,Diabetes_family_his==0&pre_DM==0)#1466
data2<-subset(data2,Tumor==0)#1446
data2<-subset(data2,antibioticUse_3month!=1)#1428
#ART
data2<-subset(data2,ART==0)#1398
#twins
data2<-subset(data2,Twins==0)#1378
#age
#data2<-subset(data2,Age<=35)
#data2<-subset(data2,Tumor==0&Eclampsia==0&HBV==0&ICP==0)#358
# 
# data2<-subset(data2,OGTT1_24!="NA")
# data2<-subset(data2,OGTT0_24!="NA")
# data2<-subset(data2,OGTT2_24!="NA")#1355

data2<-data2[complete.cases(data2[,pes_name]),]#852
data2<-subset(data2,`24wFaeces`==1) #852

#2. imputation(missForest) --------------------------------------------------------------
root_dir<-"D:/Project/pes_GM_GDM/Data_review"
#创建当前工作日期目录
#filename<-paste("result.",Sys.Date(),sep = "")
# dir.create("Result")
filename<-"Result"
dir.create(filename)
setwd(paste(root_dir,filename,sep = "/"))

dir2<-"imputation"
dir.create(dir2)
setwd(paste(source_dir_root1,filename,dir2,sep = "/"))
cov_factor<-c("pre_BMI_group","Income_group","preDrinking","SmokingHis","passsmokHis_1y",
              "Educational_level","Parity")
cov_linear<-c("Age","pre_BMI","Delivery_week","weightgain","Birthweight")

for (i in cov_factor){
  data2[,i]<-as.factor(data2[,i])
}


library(mice)
data_metadata<-data2[,c("number",cov_linear,cov_factor,"GDM")]
pdf(file = "missing_plot.pdf",height = 10,width = 8)
md.pattern(data_metadata,rotate.names = T)
dev.off()
#
library(missForest)
set.seed(123)
data_impute<-missForest(data_metadata)
data_impute$OOBerror
#NRMSE is normalized mean squared error. It is used to represent error derived from imputing continuous values. PFC (proportion of falsely classified) is used to represent error derived from imputing categorical values.
data_impute$ximp
data_imputed<-data_impute$ximp
dat<-data2[,-which(names(data2)%in%c(cov_linear,cov_factor,"GDM"))]
data3<-cbind(data_imputed,dat)
dim(data3)

#dir.create("Data")
write.csv(data3,file = "data_imputed.csv",row.names = F)

# 二、General characteristic ------------------------------------------------


# Table S1- pesticide exposure level----------------------------------------------------------------


setwd(paste(root_dir,filename,sep = "/"))
dir.create("linear")
dir.create("TableS1")
library(tableone)
library(EnvStats)
library(ggplot2)
#library( survJamda)
#compute detection rate
detect<-function(x){
  length(which(x > min(x)))*100/length(x)
}
quant<-function(x){
  quantile(x,probs = c(0,0.05,0.50,0.75,0.90,1),na.rm = T)
}

m<-function(x){
  mean(x,na.rm = T)
}

s<-function(x){
  sd(x,na.rm = T)
}

geo_m<-function(x){
  geoMean(x,na.rm = T)
}

geo_s<-function(x){
  geoSD(x,na.rm = T)
}

iqr_s<-function(x){
  IQR(x,na.rm = T)
}

dr<-apply(data2[pes_name], 2, detect)
pes_name2<-data2[,pes_name]
m1<-apply(data2[pes_name],2,m)
sd1<-apply(data2[pes_name],2,s)
iqr1<-apply(data2[pes_name], 2, iqr_s)
library(EnvStats)
gm1<-apply(pes_name2,2,geo_m)
gmd1<-apply(pes_name2,2,geo_s)
qua<-apply(pes_name2,2,quant)
qua<-t(qua)
result<-cbind(dr,m1,sd1,gm1,gmd1,qua,iqr1)
write.csv(result,file="linear/TableS1.csv")

# Table 1-general characteristic -----------------------------------------------------------

dir.create("Table1")
bas.t <- descrTable(GDM~Age+pre_BMI+weightgain+pre_BMI_group+Parity+
                      Educational_level+Income_group+#income+
                      preDrinking+SmokingHis+passsmokHis_1y+
                      OGTT1_24+OGTT1_24+OGTT2_24+
                      SEX+Delivery_week+weightgain+PTB+Birthweight+LBW,data=data2,show.p.overall = F,show.p.trend = TRUE)
export2csv(bas.t,file = 'Table1/table1_1.csv')
tab1.1 <- read.csv('Table1/table1_1.csv')



bas.t2 <- descrTable(~Age+pre_BMI+weightgain+pre_BMI_group+Parity+
                       Educational_level+Income_group+#income+
                       preDrinking+SmokingHis+passsmokHis_1y+
                       OGTT1_24+OGTT1_24+OGTT2_24+
                       SEX+Delivery_week+weightgain+PTB+Birthweight+LBW,data=data2,show.p.overall = F,show.p.trend = TRUE)
export2csv(bas.t2,file = 'Table1/table1_2.csv')
tab1.2 <- read.csv('Table1/table1_2.csv')


tab1 <- cbind(tab1.2[,1:2],tab1.1[,2:4])
colnames(tab1) <- c('Characteristic','Overall','GDM-','GDM+','P-value')
tab1
write.csv(tab1,file = 'Table1/Table1.csv',row.names = F)


# Figure 3a -----------------------------------------------------------------
setwd(source_dir_root1)
dir.create("Result/Figure3")
pes_new<-c("Mirex", "BHC_beta", "Chlorpyrifos", "Dimethoate", 
           "Parathion_methyl", "Dimethylvinphos", "Dicrotophos", 
           "Methidathion", "Monocrotophos", "Phosphamidon", 
           "Atrazine", "Metribuzin", "Myclobutanil", "Paclobutrazole", 
           "Dimethipin", "Quinoxyfen", "Diphenamid", "Vinclozolin", 
           "Monolinuron", "Dimethenamid", "Clomazone", "Propyzamide", 
           "Dicloran", "Pyrimethanil", "Furalaxyl", "Nuarimol", 
           "Isocarbamide", "Alachlor")
pes_log_new<-c("Mirex_log", "BHC_beta_log", "Chlorpyrifos_log", "Dimethoate_log", 
               "Parathion_methyl_log", "Dimethylvinphos_log", "Dicrotophos_log", 
               "Methidathion_log", "Monocrotophos_log", "Phosphamidon_log", 
               "Atrazine_log", "Metribuzin_log", "Myclobutanil_log", "Paclobutrazole_log", 
               "Dimethipin_log", "Quinoxyfen_log", "Diphenamid_log", "Vinclozolin_log", 
               "Monolinuron_log", "Dimethenamid_log", "Clomazone_log", "Propyzamide_log", 
               "Dicloran_log", "Pyrimethanil_log", "Furalaxyl_log", "Nuarimol_log", 
               "Isocarbamide_log", "Alachlor_log")


data3<-read.csv("Data/data_imputed.csv",header = T,sep = ",")
library(psych)
library(dplyr)
library(reshape2)
library(ggplot2)
library(reshape2)
library(openxlsx)
data4<-data3


description<-read.xlsx("Data/description.xlsx",sheet=1)

setwd(paste0(source_dir_root1,"/Result"))
description<-subset(description,DR>20)
data4[pes_log_new]<-log(data3[pes_new])

pes_name<-as.character(description$Pesticides)

numb<-seq(1,852,by=1)
plot_dat<-data4[,pes_log_new]
plot_dat<-data.frame(numb,plot_dat)
plot_dat<-melt(plot_dat,id.vars = "numb")

library(stringr)

t<-str_split_fixed(plot_dat$variable, "_log", 2)
t<-as.data.frame(t)

plot_dat<-data.frame(plot_dat,t)

plot_dat$V1<-factor(plot_dat$V1,levels = pes_name)
#plot_dat$variable<-factor(plot_dat$variable)
library(circlize)

circos.clear()
circos.par("track.height" = 0.25,cell.padding = c(0.02, 0, 0.02, 0))
#circos.par("track.height" = 0.25)

circos.initialize(plot_dat$V1, x = plot_dat$value)
#plot_dat$value_scale<-scale(plot_dat$value)
##add point
circos.track(plot_dat$V1, y = plot_dat$numb,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$cell.ylim[2] + mm_y(5),
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
library(xtermStyle)

col=c(rep("#96c9a5",8),rep("#dcbf2a",10),rep("#da6c7d",10))
# 
library(RColorBrewer)
circos.trackPoints(plot_dat$V1, plot_dat$value,plot_dat$numb,  col = col,
                   pch = 16, cex = 0.5)
circos.trackHist(plot_dat$V1, plot_dat$value, #bin.size = 0.2,
                 bg.col = "white", col =col)

circos.clear()



#ggsave(filename = "Figure3/Figure3a.pdf",width = 20,height = 20)


# Figure 3b ---------------------------------------------------------------


# 1.calculate pesticide risk score ----------------------------------------

###########################################################################
### Using adaptive elastic net to create environmental risk score (ERS) ###
###########################################################################

covariate<-c("pre_BMI","Age","Educational_level","weightgain","Parity","passsmokHis_1y")

d1<-data3[,c(pes_new,covariate,"GDM")]
##data proprocessing
d1$Educational_level<-as.factor(d1$Educational_level)
d1$Parity<-as.factor(d1$Parity)
d1$passsmokHis_1y<-as.factor(d1$passsmokHis_1y)

dummyVar<-model.matrix(~pre_BMI+Age+Educational_level+weightgain+Parity+passsmokHis_1y,d1)
dummyVar<-dummyVar[,-1]
ind.data<-data.frame(d1[,pes_new],dummyVar)

depen.data<-d1$GDM
preds<-28
covars<-7
a<-29:35

library(glmnet)
library(gcdnet)

aENET_fun = function(ind.data, depen.data, preds, covars, a){
  ind.data<-data.matrix(ind.data)
  depen.data<-data.matrix(depen.data)
  #set a seed of random numbers
  set.seed(111) 
  
  #five fold cross-validation list for estimating lambda 2 tuning parameter
  num.obs=nrow(ind.data)
  foldid=sample(rep(seq(5),length=nrow(ind.data))) 
  lambda2.list=seq(0,1,by=0.01) #candidate lambda2 values
  num.lambda2=length(lambda2.list)
  min.cv=numeric(num.lambda2)
  
  #creates two vectors of penalty factors that allows shrinkage of exposures, but not covariates
  pf_en<-c(rep(1,preds),rep(0,covars)) 
  pf_en_2<-c(rep(1,preds),rep(0,covars))
  
  # loop that runs each lambda2 value in the cv.gcdnet function to test for the optimal (minimum CV error) lambda2 value  
  for(i in 1:num.lambda2){
    # cv.en=cv.gcdnet(ind.data, depen.data, foldid=foldid, lambda2=lambda2.list[i], 
    #                 method="hhsvm", pred.loss="misclass", pf=pf_en, pf2=pf_en_2, standardize=T)  
    cv.en=cv.gcdnet(ind.data, depen.data, foldid=foldid, lambda2=lambda2.list[i], 
                    method="logit", pred.loss="loss",pf=pf_en, pf2=pf_en_2,standardize=T)  
    # lambda (for L1 penalty) will be automatically generated by default setting
    min.cv[i]=min(cv.en$cvm)  # collects minimum cv error for each lambda2 candidate
  }
  
  lambda2.opt=lambda2.list[which.min(min.cv)] #define optimal lambda2 that has the minimum CV error
  print(lambda2.opt)
  # now creating potential lambda1 tuning parameter 
  lambda.list=cv.en$lambda #candidate lambda (for L1 penalty)
  
  #cross validation for lambda given the optimal lambda2
  par(mar=c(5,5,4,2))
  cv.en=cv.gcdnet(ind.data, depen.data, foldid=foldid, lambda2=lambda2.opt, 
                  method="logit", pred.loss="loss", pf=pf_en, pf2=pf_en_2, standardize=T) 
  plot(cv.en)
  minval=cv.en$lambda.min
  fit.en = gcdnet(ind.data, depen.data, lambda=minval, standardize=T, 
                  lambda2=lambda2.opt, method="logit",pf=pf_en, pf2=pf_en_2)
  beta.en = coef(fit.en)
  
  ### Adaptive ENET
  v <- log(ncol(ind.data))/log(num.obs)
  gamma <- ceiling(2*v/(1-v))+1
  x.sd <- colSds(as.matrix(ind.data))*((num.obs-1)/num.obs)^0.5 # sample sd for each predictor 
  beta.enet.star <- beta.en[-1,]*x.sd
  
  ada.wts <- (abs(beta.enet.star)+1/num.obs)^(-gamma)
  ada.wts[a] <- 0
  min.cv <- numeric(num.lambda2)
  # CV for each candidate lambda2
  for(i in 1:num.lambda2){
    cv.adenet = cv.gcdnet(ind.data, depen.data, foldid=foldid, lambda2=lambda2.list[i],
                          method="logit", pf=ada.wts, pred.loss="loss", standardize=T)
    # lambda (for L1 penalty) will be automatically generated by default setting
    min.cv[i] = min(cv.adenet$cvm)  # collect minimum cv error for each lambda2 candidate
  }
  
  lambda2.adenet.opt=lambda2.list[which.min(min.cv)] 
  
  #cross validation for lambda given the optimal lambda2
  cv.adenet <- cv.gcdnet(ind.data, depen.data, foldid=foldid, lambda2=lambda2.adenet.opt,
                         method="logit", pf=ada.wts, pred.loss="loss", standardize=T) 
  plot(cv.adenet)
  minval = cv.adenet$lambda.min
  fit.all <- gcdnet(ind.data, depen.data, lambda=minval, 
                    lambda2=lambda2.adenet.opt, pf=ada.wts, method="logit", standardize=T ) 
  beta.all <- as.matrix(coef(fit.all))# get the coefficients
  beta.all
  
  ##standard error and p-value for each non-zero estimated coefficient
  index.a=which(beta.all[-1]!=0)#the coordinate # of the non-zero beta.adenet except the intercept
  p.ne0=length(index.a)#the number of non-zero estimated-effect predicators
  x.wi=as.matrix(cbind(1,ind.data)) #the design matrix with first column being 1
  sigma.2=t(depen.data-x.wi%*%beta.all)%*%(depen.data-x.wi%*%beta.all)/(num.obs-p.ne0-1)#estimated sigma squared
  
  x.mean=colMeans(ind.data)
  x.stdi=matrix(0,num.obs,ncol(ind.data))#standardized x
  colnames(x.stdi)=colnames(ind.data)
  for(i in 1:ncol(ind.data)) {
    x.stdi[,i]=(ind.data[,i]-x.mean[i])/x.sd[i]  # standardized x
  }
  x.stdi.a=x.stdi[,index.a]  #standardized x matrix for non-zero estimated-effect predicators
  Sig.a.stdi=t(x.stdi.a)%*%x.stdi.a  #Sigma_A in the Theorem 3.3 of Zou and Zhang (2009)      
  
  ##variance matrix for the non-zero estimated coefficents for standardized predicators
  var.a.stdi=as.numeric(sigma.2)*(1+lambda2.adenet.opt/2)^2*solve(Sig.a.stdi+diag(rep(num.obs*lambda2.adenet.opt,p.ne0))+ (num.obs*lambda2.adenet.opt/2)^2*solve(Sig.a.stdi))
  
  x.sd.a.inv=diag(1/x.sd[index.a])
  var.a.ori=x.sd.a.inv%*%var.a.stdi%*%t(x.sd.a.inv)   #variance matrix for the non-zero estimated coefficents for original predicators
  var.int=t(as.matrix(x.mean[index.a]))%*%var.a.ori%*%as.matrix(x.mean[index.a]) # variance for the intercept
  
  name.x.a=c("(Intercept)",colnames(ind.data)[index.a])# the names of variabls have non-zero estimated coefficient including intercept
  se.a=data.frame(se=sqrt(c(var.int,diag(var.a.ori))),row.names=name.x.a) #Standard error of non-zero beta.adenet including the intercept
  
  pvalue.a=data.frame(pvalue=2-2*apply(abs(beta.all[c(1,index.a+1)])/se.a,1,pnorm),row.names=name.x.a)
  
  result.adenet=cbind(numeric(ncol(ind.data)+1),matrix(NA,ncol(ind.data)+1,2))
  colnames(result.adenet)=c("beta","SE","p-value")
  rownames(result.adenet)=c("(Intercept)",colnames(ind.data))
  result.adenet[name.x.a,1]=beta.all[name.x.a,]
  result.adenet[name.x.a,2]=se.a[name.x.a,]
  result.adenet[name.x.a,3]=pvalue.a[name.x.a,]
  result.adenet = as.data.frame(result.adenet)
  result.adenet$loci = result.adenet$beta - 1.96*result.adenet$SE
  result.adenet$hici = result.adenet$beta + 1.96*result.adenet$SE
  return(result.adenet)
}

results<-aENET_fun(ind.data, depen.data, preds, covars, a)

#Extract enet beta coefficients for weights
nexp<-28 #define the number of exposures
weights <- (results[,1])[2:(nexp+1)]

#isolate your log-transformed or standardized exposure variables
exp_mat <- as.matrix(d1[,colnames(d1) %in% rownames(results)[2:(nexp+1)]])

#Construct environmental risk score by applying linear combination of weights and exposure matrix
ers <- exp_mat%*%as.matrix(weights, ncol = 1)

#Combine the ERS to your dataset
d1 <- cbind(d1, ers)

summary(glm(GDM~ers,family = "binomial",data=d1))
#dir.create(paste0(paste0(path_med_16S,"/ERS")))
#dir.create("linear/ERS")
write.csv(d1,file = "Figure3/data_ers.csv",row.names = F)
write.csv(result.adenet,file = "Figure3/result.adenet.csv")

##variable selection plot
dat_plot<-read.csv(paste0("Figure3/result.adenet.csv"),header = T,sep = ",")
dat_plot<-dat_plot[-1,]
dat_plot$Group<-rep("Coefficient",times=nrow(dat_plot))
dat_plot$beta<-ifelse(dat_plot$beta==0,NA,dat_plot$beta)
order<-as.character(dat_plot$X)
#str(sig.mat)
ggplot(dat_plot, aes(x = X, y = Group, fill = beta)) +
  # geom_tile(color = "black") +
  geom_tile() +
  # geom_text(
  #   aes(label = sig.mat),
  #   color = "white", size = 6) +
  scale_fill_gradient(low="white",high = "#43997A",na.value = 'white')+
  
  theme_bw()+
  scale_x_discrete(limits=order)+
  theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
  labs(x="",y="")
# coord_fixed() 
#dir.create("plot/vari_select")
ggsave(paste0("Figure3/Figure3b.pdf"),height = 4,width = 10)  


###errorbar

data_ers<-read.csv(paste0("Figure3/data_ers.csv"))

# metals_log_new<-c("Fe_log", "Cu_log", "Zn_log", "Mo_log", "Co_log", "V_log", 
#                   "As_log", "Cd_log", "Pb_log", "Hg_log", "Sr_log", "Rb_log"
# )
#dat_plot<-subset(dat_plot,i %in% pes_log_new)



data_ers$GDM<-as.factor(data_ers$GDM)
p<-ggviolin(data = data_ers,x = "GDM",y = "ers",fill="GDM") +
  #stat_compare_means(method = "t.test")+
  scale_fill_manual(values=c("#d95f02", "#1b9e77"))+
  #geom_jitter(size=0.3) +
  labs(x="",y="Environmental Risk Score")+
  scale_x_discrete(limits=c("0","1"),labels=c("Control","GDM"))+
  theme_bw()+
  theme(legend.position = "none")
ggsave(paste0("Figure3/boxplot.pdf"),height = 4,width = 4)  



# Figure 3c ---------------------------------------------------------------
## 1.1 logistic regression -----------------------------------------------------
#all_variable<-c(pes_new,pes_log_new)
#pes
#data3$GDM<-as.factor(data3$GDM)

d2<-data.frame(d1,data3[,c("IGT","OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")])

library(fdrtool)
result<-data.frame()
result2<-data.frame()
for (j in c("GDM")){
  for(i in "ers"){
    #unadjusted
    fit_i<-glm(as.formula(paste(j,"~",i)),family = "binomial",data=d2)
    p_i<-sprintf("%.3f",tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"}))
    b_i<-sprintf("%.3f",tryCatch({exp(stats::coef(fit_i)[2])},error=function(e){"NA"}))
    blow_i<-sprintf("%.3f",tryCatch({exp(confint(fit_i)[2,1])},error=function(e){"NA"}))
    bhigh_i<-sprintf("%.3f",tryCatch({exp(confint(fit_i)[2,2])},error=function(e){"NA"}))
    b1<-paste(b_i," (",blow_i,",",bhigh_i,")",sep = "")
    
    # xnam_2i<-paste(i,"pre_BMI","Age","Educational_level","Parity",sep = "+")
    # fit_2i<-glm(as.formula(paste("GDM~",xnam_2i)),family = "binomial",data=data3)
    # p_2i<-sprintf("%.3f",tryCatch({summary(fit_2i)$coefficients[2,4]},error=function(e){"NA"}))
    # b_2i<-sprintf("%.3f",tryCatch({exp(coef(fit_2i)[2])},error=function(e){"NA"}))
    # blow_2i<-sprintf("%.3f",tryCatch({exp(confint(fit_2i)[2,1])},error=function(e){"NA"}))
    # bhigh_2i<-sprintf("%.3f",tryCatch({exp(confint(fit_2i)[2,2])},error=function(e){"NA"}))
    # b2<-paste(b_2i," (",blow_2i,",",bhigh_2i,")",sep = "")
    
    xnam_3i<-paste(j,i,"pre_BMI","Age","Educational_level","weightgain","Parity","passsmokHis_1y",
                   sep = "+")
    fit_3i<-glm(as.formula(paste(j,"~",xnam_3i)),family = "binomial",data=d2)
    p_3i<-sprintf("%.3f",tryCatch({summary(fit_3i)$coefficients[2,4]},error=function(e){"NA"}))
    b_3i<-sprintf("%.3f",tryCatch({exp(stats::coef(fit_3i)[2])},error=function(e){"NA"}))
    blow_3i<-sprintf("%.3f",tryCatch({exp(confint(fit_3i)[2,1])},error=function(e){"NA"}))
    bhigh_3i<-sprintf("%.3f",tryCatch({exp(confint(fit_3i)[2,2])},error=function(e){"NA"}))
    b3<-paste(b_3i," (",blow_3i,",",bhigh_3i,")",sep = "")
    
    result_i<-data.frame(j,i,b1,p_i,#b2,p_2i,
                         b3,p_3i)  
    result<-rbind(result,result_i)
    
    result_2i<-data.frame(j,i,b_i,blow_i,bhigh_i,p_i,
                          #b_2i,blow_2i,bhigh_2i,p_2i,
                          b_3i,blow_3i,bhigh_3i,p_3i)
    result2<-rbind(result2,result_2i)
  }
}


names(result)<-c('outcomes','Pesticides',"OR1","p1","OR3","p3")
# write.table(result,file="linear/result_bino.txt",sep="\t",row.names = F)
# write.table(result2,file="linear/result_bino2.txt",sep="\t",row.names = F)

write.csv(result,file=paste0(path_med_16S,"/ERS/result_ers.csv"),row.names = F)
write.csv(result2,file=paste0(path_med_16S,"/ERS/result_ers2.csv"),row.names = F)


# Figure S4c-Quantile g-computation -----------------------------------------------------------
setwd(paste(root_dir,filename,sep = "/"))

dir.create("FigureS4")
mod2<-c("pre_BMI","Age","Educational_level","Parity","weightgain","passsmokHis_1y")
#outcome<-c("HbAlc_24","OGTT0_24","OGTT1_24","OGTT2_24","GDM","IGT")

pes_new<-description$Pesticides[description$DR>20]
pes_upper<-description$Pesticides[description$DR>80]
pes_log_new<-paste(pes_new,"_log",sep = "")
data3[pes_log_new]<-log(data3[pes_new])

#dir.create("G_computation")
package_list <- c("qgcomp","openxlsx","ggplot2","knitr")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    # install.packages(p, repos=site)
    install.packages(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# we save the names of the mixture variables in the variable "Xnm"
Xnm <- pes_new
Xnm_log<-paste(Xnm,"_log",sep = "")
covariate<-c("pre_BMI","Age","Educational_level","weightgain","Parity","passsmokHis_1y")

gcom<-function(i){
  qc.fit <- qgcomp.noboot(as.formula(paste(i,"~.",sep = "")),dat=data3[,c(Xnm, i)],
                          family=binomial(), q=4)
  
  b<-sprintf("%.3f",exp(qc.fit$psi))
  #b_ci<-qc.fit$ci.coef[2,]
  b_ci<-sprintf("%.3f",exp(qc.fit$ci))
  p_value<-qc.fit$pval[2]
  
  result<-c(b,b_ci,p_value)
  result<-as.data.frame(result)
  result<-as.data.frame(t(result))
  names(result)<-c("beta","LCI","UCI","pvalue")
  
  weight_p<-qc.fit$pos.weights
  pes_pos<-names(weight_p)
  Group<-rep("Positive",times=length(pes_pos))
  weight_pos<-data.frame(pes_pos,weight_p,Group)
  names(weight_pos)<-c("pes","weight","Group")
  
  weight_n<-qc.fit$neg.weights
  pes_neg<-names(weight_n)
  Group<-rep("Negative",times=length(pes_neg))
  weight_neg<-data.frame(pes_neg,weight_n,Group)
  
  names(weight_neg)<-c("pes","weight","Group")
  
  weight_combine<-rbind(weight_pos,weight_neg)
  
  write.csv(weight_combine,file = paste0("FigureS4/",i,"_weight.csv"),row.names = F)
  write.csv(result,file = paste0("FigureS4/",i,"_result.csv"),row.names = F)
  pdf(file =paste0("FigureS4/",i,"_weight.pdf"),width = 6,height = 6 )
  p = plot(qc.fit)
  dev.off()
  #ggsave(p,filename = "G_computation/weight_plot.pdf",width = 6,height = 8)
  
  qcboot.fit5 <- qgcomp(as.formula(paste(i,"~.+.^2")),
                        expnms=Xnm,
                        na.omit(data3[,c(Xnm, i)]), family=binomial(), q=4, degree=2,
                        B=10, rr=FALSE, seed=125)
  # qgcomp::pointwisebound.boot(qcboot.fit5)
  pdf(file =paste0("FigureS4/",i,"_joint_CI.pdf"),width = 6,height = 6 )
  plot(qcboot.fit5)
  dev.off()
}



for (i in c("GDM")) {
  # gcom_log(i)
  gcom(i)
}



# Figure S2a --------------------------------------------------------------
setwd(paste(root_dir,filename,sep = "/"))

dir.create("FigureS2")

#library(ALL)
# data(ALL)
# d=exprs(ALL)
data <- data3[pes_log_new]
data <- na.omit(data)
#data <- scale(data)
d<-t(data)
#筛选前5000标准差的基因
#mads=apply(d,1,mad)
#d=d[rev(order(mads))[1:5000],]

#sweep函数减去中位数进行标准化
d = sweep(d,1, apply(d,1,median,na.rm=T))
d<-as.matrix(d)
#一步完成聚类
library(ConsensusClusterPlus)
dir.create("Concensus_hc")
title="Concensus_hc"
results = ConsensusClusterPlus(d,maxK=10,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="spearman",seed=1262118388.71279,plot="pdf")

lable<-as.data.frame(results[[2]][["consensusClass"]])


library(ggplot2)
library(reshape2)
#data<-read.table("clipboard",header=T)
d2<-data
names(lable)<-"Clust"
d2<-cbind(d2,lable)
data_new<-melt(d2,id="Clust")
write.csv(d2,file = "FigureS2/pes_sample_cluster.csv",row.names = T)


###group_compare

library(tableone)
library(compareGroups)
data_new<-d2
vars2<-paste(pes_log_new,collapse = "+")
bas.t <- descrTable(as.formula(paste("Clust~",vars2,sep = "")),data=data_new,show.p.overall = F,show.p.trend = TRUE)
export2csv(bas.t,file = 'FigureS2/exposure_k2.csv')


library(reshape2)

data_total<-cbind(data3,d2)

summary(glm(GDM~Clust,data=data_total))
summary(lm(OGTT0_24~Clust,data=data_total))
summary(lm(OGTT1_24~Clust,data=data_total))
summary(lm(OGTT2_24~Clust,data=data_total))
summary(lm(HbAlc_24~Clust,data=data_total))

plot_data=melt(data_total[,c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24","Clust")],id.vars = "Clust")
#Group<-paste("Group",plot_data$N_cluster,sep = "")
plot_data$Clust<-as.factor(plot_data$Clust)

library(ggplot2)
library(ggpubr)

# Figure S3 ---------------------------------------------------------------
setwd(paste(root_dir,filename,sep = "/"))

dir.create("FigureS3")

#boxplot

plot_data=melt(data_total[,c(pes_log_new,"Clust")],id.vars = "Clust")
#Group<-paste("Group",plot_data$N_cluster,sep = "")
plot_data$Clust<-as.factor(plot_data$Clust)
# one box per variety


library(ggplot2)
library(ggpubr)

pdf(file = "FigureS3/FigureS3.pdf",height = 20,width = 12)
ggviolin(plot_data, x="Clust", y="value",color = "Clust",add="boxplot",add.params = list(fill = "white"))+
  stat_compare_means(aes(group=Clust),label = "p.signif")+
  labs(x="",y="Level (mmol/L)")+
  #scale_color_discrete(labels=c("Low", "Medium", "High"))+
  #  scale_fill_npg()+
  scale_color_brewer(palette='Dark2')+
  #  scale_x_discrete(limits=id_name,label=id_name)+
  #theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
  theme(panel.background = element_rect(color = "black"),
        legend.title=element_blank())+
  facet_wrap(~variable,scales = "free_y",ncol = 4)
# scale_color_lancet()
dev.off()


# ##cluster and GDM


# Figure S4b --------------------------------------------------------------
setwd(paste(root_dir,filename,sep = "/"))

dir.create("FigureS4")


## LASSO regression ------------------------------------------------------
method="binomial"

#dir.create("linear/LASSO")
#log
Lasso_selection_log<-function(outcome){
  data4<-data3[complete.cases(data3[,c(outcome,pes_log_new)]),]
  set.seed(2021)
  training.samples <- data4[,outcome] %>%
    createDataPartition(p = 0.8, list = FALSE)
  train.data  <- data4[training.samples, ]
  test.data <- data4[-training.samples, ]
  
  # Predictor variables
  x <- as.matrix(train.data[,pes_log_new])
  # Outcome variable
  y <- train.data[,outcome]
  # Find the best lambda using cross-validation
  set.seed(2021) 
  cv <- cv.glmnet(x, y, alpha = 1,family=method)
  # Fit the final model on the training data
  model1 <- glmnet(x, y, alpha = 1, lambda = cv$lambda.min)
  model2 <- glmnet(x, y, alpha = 1, lambda = cv$lambda.1se)
  # Dsiplay regression coefficients
  result_coef<-stats::coef(model2)
  # re<-data.frame(result_coef@Dimnames[[1]])
  # value<-result_coef@x
  name<-result_coef@Dimnames[[1]]
  result_summary<-summary(result_coef)#delete constant variable
  remain_pes<-name[result_summary$i]
  remain_pes<-as.data.frame(remain_pes)
  group<-rep(outcome,times=dim(remain_pes)[1])
  result<-data.frame(remain_pes,group)
  write.table(result,file = paste("FigureS4/",outcome,"_",method,"_log.txt"))
  
}

for (i in "GDM"){
  Lasso_selection_log(i)
  #wqs_function(i)
}


## combine result ------------------------------------------------------

method<-"FigureS4"
library(dplyr)
library(tidyr)
#  file_function<-function(time){
source_dir<-paste(source_dir_root1,filename,method,sep = "/")
files=list.files(path = source_dir, pattern = "*_log.txt", all.files = FALSE,
                 full.names = FALSE, recursive = FALSE,
                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
setwd(source_dir)
data<-data.frame()
for (i in files){
  item_i<-read.table(i,sep="\t",header = F)
  item_2i <- item_i %>% separate(V1, c("index","pes","outcome"), "[ ]")
  item_i<-item_2i[-c(1:2),-1]
  #result_i<-as.data.frame(item_i)
  #result_i<-as.data.frame(t(result_i))
  data<-rbind(data,item_i)
}
data$method<-rep(method,times=dim(data)[1])
#  data$Time<-rep(time,times=dim(data)[1])
names(data)<-c("pes","outcome","method")

write.csv(data,file = "result_LASSO.csv")




data_lasso<-read.csv("result_LASSO.csv",header = T,sep = "," )
# data3<-read.csv("general/data3.csv",header = T,sep = "," )
# dir.create("result")
# dir.create("result/lasso")
data_lasso<-data_lasso[,-1]
data_lasso<-data_lasso[!duplicated(data_lasso),]
result1<-data.frame()



i<-"GDM"
result1<-data.frame()
for (i in c("GDM")){
  new_dat<-subset(data_lasso,outcome==i)
  cova<-paste(new_dat$pes,collapse = "+")
  cova2<-c(as.character(new_dat$pes),"pre_BMI","Age","Educational_level","weightgain","Parity","passsmokHis_1y")
  cov2<-paste(cova2,collapse = "+")
  # data3[,i]<-as.numeric(as.character(data3[,i]))
  fit_i<-glm(as.formula(paste(i,"~",cova)),family = binomial(),data=data3)
  
  p_i<-tryCatch({summary(fit_i)$coefficients[2:(dim(new_dat)[1]+1),4]},error=function(e){"NA"})
  OR_i<-tryCatch({exp(stats::coef(fit_i)[2:(dim(new_dat)[1]+1)])},error=function(e){"NA"})
  ORlow_i<-tryCatch({exp(confint(fit_i)[2:(dim(new_dat)[1]+1),1])},error=function(e){"NA"})
  ORhigh_i<-tryCatch({exp(confint(fit_i)[2:(dim(new_dat)[1]+1),2])},error=function(e){"NA"})
  
  fit_2i<-glm(as.formula(paste(i,"~",cov2)),family = binomial(),data=data3)
  p_2i<-tryCatch({summary(fit_2i)$coefficients[2:(dim(new_dat)[1]+1),4]},error=function(e){"NA"})
  OR_2i<-tryCatch({exp(stats::coef(fit_2i)[2:(dim(new_dat)[1]+1)])},error=function(e){"NA"})
  ORlow_2i<-tryCatch({exp(confint(fit_2i)[2:(dim(new_dat)[1]+1),1])},error=function(e){"NA"})
  ORhigh_2i<-tryCatch({exp(confint(fit_2i)[2:(dim(new_dat)[1]+1),2])},error=function(e){"NA"})
  
  name<-rep(i,times=dim(new_dat)[1])
  result_i<-data.frame(name,OR_i,ORlow_i,ORhigh_i,p_i,
                       OR_2i,ORlow_2i,ORhigh_2i,p_2i)  
  result1<-rbind(result1,result_i) 
  write.csv(result1,file = paste("lasso_bivar.csv"),row.names = T)
  
}


## DSA -------------------------------------------------------------------
setwd(paste0(source_dir_root1,"/",filename))
library(DSA)
# devtools::install_github("romainkp/modelUtils")
# devtools::install_github("romainkp/DSA")
dir.create("FigureS4")
dat2<-data3
exposure<-dat2[pes_log_new]
covariate<-c('pre_BMI',"Age","weightgain","Educational_level","Parity","passsmokHis_1y")

#setwd('C:/Users/18896/Desktop/neuro/DSA')
#data3<-dat[,c(cov,pes,outcome_variable)]
setwd(paste(getwd(),"FigureS4",sep = "/"))
library('DSA')
library('parallel')
library('doParallel')
write(paste0("prepare output -DSA") , stderr())
# preparing output
# execute 100 DSA selection among all imputed dataset
t1_50_pre<-Sys.time() 
write("bucle", stderr())
## Calculate the number of cores
##set the number of iterations and number of cores
num_iter <- 50 + 1 
no_cores <- 6


data4<-dat2[,c("GDM","Age","pre_BMI","Educational_level","weightgain","Parity",
               "passsmokHis_1y",pes_log_new)]
library('DSA')
library('parallel')
library('doParallel')
# res <- DSA(GDM ~ M.age+BMI, data = data4, maxsize = 10, maxsumofpow = 1, maxorderint = 1)
# residuals(res)
# coefficients(res)
# summary(res)

write(paste0("prepare output -DSA") , stderr())
# preparing output

# execute 100 DSA selection among all imputed dataset
t1_50_pre<-Sys.time() 

write("bucle", stderr())

## Calculate the number of cores

##set the number of iterations and number of cores
num_iter <- 50 + 1 
no_cores <- 6


# Initiate cluster
cl <- makeCluster(no_cores)

stacked_d<-c(covariate,"GDM",pes_log_new)
clusterExport(cl=cl, varlist=c("stacked_d"))
clusterEvalQ(cl, library('DSA'))


# tests.dsa<-parLapply(cl, 1:num_iter, function(iter_pre, stacked_d){
#   DSA(GDM~M.age+BMI,  family = binomial,data=stacked_d, maxsize=30, maxorderint=1, maxsumofpow=1, id=stacked_d[,"newid"])
# }, stacked_d)
data4$Educational_level<-as.numeric(data4$Educational_level)
data4$Parity<-as.numeric(data4$Parity)
data4$passsmokHis_1y<-as.numeric(data4$passsmokHis_1y)

tests.dsa<-parLapply(cl, 1:num_iter, function(iter_pre, data4){
  DSA(GDM~Age+pre_BMI+weightgain+Educational_level+Parity+passsmokHis_1y, family = binomial(),data=data4, maxsize=30, maxorderint=1, maxsumofpow=1)
}, data4)

stopCluster(cl)

#save.image(paste0('models_DSA_', format(Sys.Date(), "%d%b%Y"),'.Rdata'))


##we take the first test as reference
mod_GDM1_pre<- tests.dsa[[1]]

tot_var_GDM_pre<-mod_GDM1_pre$candidate.vars
tot_var_GDM2_pre<-as.data.frame(tot_var_GDM_pre)

# set nb to 0 for all variables
tot_var_GDM2_pre$nb_all<-0


##conversión of all the models to text to get the frequency of the variables, except the first test.
test.models<-lapply(tests.dsa[-1], function(x){ as.character(x$model.selected)[3]})
## frequency count
tot_var_GDM2_pre$nb_all<-sapply(tot_var_GDM2_pre$tot_var_GDM_pre, function(var, test.models){ length(grep(var,test.models))}, test.models)

t2_50_pre<-Sys.time() 
t2_50_pre-t1_50_pre  
#look at the results (ie, variable selected at least in one out of 50 DSA)

# sink(paste0('result_pal_hs_mvpa_GDM', format(Sys.Date(), "%d%b%Y"), '.txt'))
# tot_var_GDM2_pre[tot_var_GDM2_pre$nb_all>0,]
# sink()

varia<-rep("GDM",times=nrow(tot_var_GDM2_pre))
tot_var_GDM2_pre<-data.frame(tot_var_GDM2_pre,varia)
write.table(tot_var_GDM2_pre,file = "GDM.txt",sep = "\t")
# 


#  combine result 
#source_dir_root1<-"E:/study2/microbiome/GDM2/pes_GM_GDM/result.2021-09-21"
method<-"DSA"
library(dplyr)
library(tidyr)
#file_function<-function(time){
source_dir<-paste(source_dir_root1,filename,method,sep = "/")
files=list.files(path = source_dir, pattern = "*.txt", all.files = FALSE,
                 full.names = FALSE, recursive = FALSE,
                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
setwd(source_dir)
data<-data.frame()
for (i in files){
  item_i<-read.table(i,sep="\t",header = T)
  names(item_i)<-c("Var","n","Outcome")
  # item_i$tot_var_bmi_pre.nb_all<-as.character(item_i$tot_var_bmi_pre.nb_all)
  # item_2i <- item_i %>% separate(tot_var_bmi_pre.nb_all, c("number","pes","count"), "[ ]")
  # item_i<-item_2i[-c(1:2),-1]
  #result_i<-as.data.frame(item_i)
  #result_i<-as.data.frame(t(result_i))
  data<-rbind(data,item_i)
}
data$method<-rep(method,times=dim(data)[1])
#data$Time<-rep(time,times=dim(data)[1])
names(data)<-c("variable","num","outcome","method")

#setwd(paste(source_dir_root1,filename,sep = "/"))
write.csv(data,file = "result_DSA.csv")

#  setwd("E:/neuro/result.2021-08-16")
library(openxlsx)
#data<-read.table("clipboard",header = T,sep = "\t")
data_DSA<-read.csv("result_DSA.csv",header = T,sep = "," )
# data3<-read.csv("general/data3.csv",header = T,sep = "," )
# dir.create("result")
# dir.create("result/enet")
result1<-data.frame()

#outcomes_bivar<-c("GDM","IGT")
# t<-paste0(as.character(unique(data_DSA$outcome)),"_24")
# outcome_variable<-intersect(outcome_linear,t)

result1<-data.frame()
outcome_bivar2<-c("GDM")

for (i in outcome_bivar2[1]){
  new_dat<-subset(data_DSA,outcome==i)
  new_dat<-subset(new_dat,num>=5)
  met<-intersect(new_dat$variable,pes_log_new)
  cova<-paste(met,collapse = "+")
  cova2<-c(met,"pre_BMI","Age","Educational_level","weightgain","Parity","passsmokHis_1y")
  cov2<-paste(cova2,collapse = "+")
  fit_i<-glm(as.formula(paste(i,"~",cova)),family = binomial(),data=data3)
  
  p_i<-tryCatch({summary(fit_i)$coefficients[2:(length(met)+1),4]},error=function(e){"NA"})
  b_i<-tryCatch({exp(stats::coef(fit_i)[2:(length(met)+1)])},error=function(e){"NA"})
  blow_i<-tryCatch({exp(confint(fit_i)[2:(length(met)+1),1])},error=function(e){"NA"})
  bhigh_i<-tryCatch({exp(confint(fit_i)[2:(length(met)+1),2])},error=function(e){"NA"})
  
  fit_2i<-glm(as.formula(paste(i,"~",cov2)),family = binomial(),data=data3)
  p_2i<-tryCatch({summary(fit_2i)$coefficients[2:(length(met)+1),4]},error=function(e){"NA"})
  b_2i<-tryCatch({exp(stats::coef(fit_2i)[2:(length(met)+1)])},error=function(e){"NA"})
  blow_2i<-tryCatch({exp(confint(fit_2i)[2:(length(met)+1),1])},error=function(e){"NA"})
  bhigh_2i<-tryCatch({exp(confint(fit_2i)[2:(length(met)+1),2])},error=function(e){"NA"})
  
  name<-rep(i,times=length(b_i))
  result_i<-data.frame(name,b_i,blow_i,bhigh_i,p_i,
                       b_2i,blow_2i,bhigh_2i,p_2i)  
  result1<-rbind(result1,result_i) 
  write.csv(result1,file = paste("DSA_bino.csv"),row.names = T)
  
}

setwd(paste0(source_dir_root1,"/",filename))
##LASSO
LASSO_linear<-read.csv("FigureS4/lasso_bivar.csv",header = T,sep = ",")
LASSO_linear<-data.frame(LASSO_linear,str_split_fixed(LASSO_linear$X, "_", 2))

LASSO_linear<-LASSO_linear[,c(2,11,7:10)]
names(LASSO_linear)[1:2]<-c("outcome","pesticides")


dat_m<-data.frame(pes_new,rep(0,times=length(pes_new)))

names(dat_m)<-c("pesticides","value")
dat_full<-data.frame()
for (i in "GDM"){
  dat_plot<-subset(LASSO_linear,outcome==i)
  
  dat_plot_full<-left_join(dat_m,dat_plot,suffix=c("pesticides","pesticides"))
  outcomes<-rep(i,times=nrow(dat_plot_full))
  dat_plot_full<-data.frame(outcomes,dat_plot_full)
  dat_full<-rbind(dat_full,dat_plot_full)
  
}



##combine result
dsa_linear<-read.csv("FigureS4/DSA_bino.csv",header = T,sep = ",")
dsa_linear<-data.frame(dsa_linear,str_split_fixed(dsa_linear$X, "_", 2))

dsa_linear<-dsa_linear[,c(2,11,7:10)]
names(dsa_linear)[1:2]<-c("outcome","pesticides")


#dat_m<-data.frame(metals,rep(0,times=length(metals)))

dat_full_dsa<-data.frame()
for (i in "GDM"){
  dat_plot<-subset(dsa_linear,outcome==i)
  
  dat_plot_full<-left_join(dat_m,dat_plot,suffix=c("pesticides","pesticides"))
  outcomes<-rep(i,times=nrow(dat_plot_full))
  dat_plot_full<-data.frame(outcomes,dat_plot_full)
  dat_full_dsa<-rbind(dat_full_dsa,dat_plot_full)
  
}
names(dat_full_dsa)<-colnames(dat_full)
dat_comb<-rbind(dat_full,dat_full_dsa)
Group<-c(rep("LASSO",times=nrow(dat_full)),rep("DSA",times=nrow(dat_full_dsa)))

dat_comb<-data.frame(dat_comb,Group)


dat_comb$logp<--log(dat_comb$p_2i)

# install.packages("ggplot2")
library(ggplot2)
#for (i in c("BW","gesweek")) {
#dat_plot<-subset(dat_comb,outcomes=)
dat_comb$pesticides<-as.character(dat_comb$pesticides)
dat_comb$Group<-as.character(dat_comb$Group)
dat_comb$outcomes<-as.character(dat_comb$outcomes)

getSig <- function(dc) {
  sc <-''
  if (is.na(dc) ) sc<-''
  else if (dc < 0.001) sc <- '***'
  else if (dc < 0.01) sc <- '**'
  else if (dc < 0.05) sc <- '*'
  sc
}

sig.mat<-purrr::map_chr(dat_comb$p_2i, getSig)
#str(sig.mat)
ggplot(dat_comb, aes(x = pesticides, y = Group, fill = logp)) +
  # geom_tile(color = "black") +
  geom_tile() +
  geom_text(
    aes(label = sig.mat),
    color = "white", size = 6) +
  scale_fill_gradient(low="white",high = "red",na.value = 'white')+
  facet_wrap(~outcomes,nrow = 1)+
  
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))
# coord_fixed() 
#dir.create("plot/vari_select")
ggsave( "FigureS4/FigureS4b.pdf",height = 4,width = 10)  



# Figure S4a --------------------------------------------------------------

#setwd("E:/study2/microbiome/GDM/pes_GM_GDM")
setwd(paste(source_dir_root1,filename,sep = "/"))
#data3<-read.csv("imputation/data_imputed.csv",header = T,sep = ",")

#data3[pes_log_new]<-log(data3[pes_new])

library(rexposome)
library(ggplot2)
outcomes<-c("HbAlc_24","OGTT0_24","OGTT1_24","OGTT2_24","GDM","IGT")
covariate<-c("pre_BMI","Age","weightgain","Educational_level","Parity","passsmokHis_1y")

dat2<-data3[,c("number",covariate,outcomes,pes_log_new)]
colnames(dat2)<-c("number",covariate,outcomes,pes_new)

dat2$GDM<-as.numeric(as.character(dat2$GDM))

loadExposome_plain <- function(data, data_id, sep = ",", pheno_cols, 
                               na.strings = c("NA", "-", "?", " ", ""),
                               families = NULL, exposures.asFactor = 5, warnings = TRUE){
  
  if(class(data) == "data.frame"){is_path <- FALSE}
  else if(class(data) == "character"){is_path <- TRUE}
  
  if(is_path){
    data <- utils::read.table(data, header = TRUE,
                              row.names = data_id, sep = sep, na.strings = na.strings)
  }
  
  exposures <- data[, !(colnames(data) %in% pheno_cols)]
  phenotype <- data[, pheno_cols]
  
  if(is.null(families)){
    description <- data.frame(Family = colnames(exposures), Exposure = colnames(exposures), Name = NA)
  }
  else{
    # Check that all exposures are present on the supplied families list,
    # if they are not, the exposomeSet creation will fail
    fams <- unlist(families)
    if(!all(fams %in% colnames(exposures))){
      stop("Missing exposures on the families list")
    }
    
    description <- reshape2::melt(do.call(cbind, families))[,2:3]
    colnames(description) <- c("Family", "Exposure")
    description <- cbind(description, Name = NA)
  }
  rownames(description) <- description$Exposure
  
  exposome <- loadExposome( exposures = exposures, description = description,
                            phenotype = phenotype, description.famCol = "Family", exposures.asFactor,
                            warnings )
  return(exposome)
} 

exp_com <- loadExposome_plain(dat2, data_id = "number",
                              pheno_cols = c("number",covariate,outcomes))




def_ewas_binomial<-function(outcome){
  fl_ew <- exwas(exp_com, formula = as.formula(paste(outcome,"~ Age + pre_BMI + Educational_level+
                                                     Parity+passsmokHis_1y")), family = "binomial")
  #head(extract(fl_ew))
  clr <- rainbow(length(familyNames(exp_com)))
  names(clr) <- familyNames(exp_com)
  
  p_effect<-plotExwas(fl_ew,color = clr) + 
    ggtitle(paste("Univariate Approach-",outcome,sep = ""))+
    geom_vline(xintercept = 1.3,color="red",linetype=2)
  ggsave(p_effect,filename = paste("FigureS4/",outcome,"_P_value.pdf",sep = ""),height = 8,width=12)
  
  p<-plotEffect(fl_ew)+geom_vline(xintercept = 0,color="red",linetype=2)+
    ggtitle(paste("Univariate Approach-",outcome,sep = ""))
  ggsave(p,filename = paste("FigureS4/",outcome,"_P_effect_FigureS4a.pdf",sep = ""),height = 8,width=8)
  
  exp_std <- standardize(exp_com, method = "normal")
  bl_mew <- mexwas(exp_std, phenotype = outcome, family = "binomial")
  
  p<-plotExwas(bl_mew) + ylab("") + ggtitle(paste("Multivariate Approach-",outcome,sep = ""))
  ggsave(p,filename = paste("FigureS4/",outcome,"_ENET.pdf",sep = ""),height = 8,width=8)
}

def_ewas_binomial('GDM')


