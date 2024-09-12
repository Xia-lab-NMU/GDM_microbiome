
source_dir<-"D:/Project/pes_GM_GDM/Data_review"
setwd(source_dir)
#load("Result/mediation/Med.RData")
#setwd(source_dir_root1)

#data3<-read.csv("imputation/data_imputed.csv")

pes<-c("Carbophenothion_sulfone", "Chlorpyrifos",  "Dimethoate", "Parathion_methyl", 
       "Dicrotophos", "Dimethylvinphos", "Methidathion", 
       "Monocrotophos",  "Phosphamidon",  
       "Mirex", "BHC_beta", "Myclobutanil", "Paclobutrazole", "Atrazine", 
       "Metribuzin", "Dimethipin", "Diphenamid", "Quinoxyfen", 
       "Isocarbamide", 
       "Vinclozolin", "Dicloran", "Clomazone", "Dimethenamid", "Propyzamide", 
       "Pyrimethanil", "Furalaxyl", "Nuarimol", "Alachlor","Monolinuron")
pes_new<-pes[-1]
pes_log_new<-paste0(pes_new,"_log")

library(gcdnet)
library(GMCM)
library(matrixStats)
library(stats)  
library(epicalc)
library(Hmisc)
library(Epi)
library(mgcv)
library(openxlsx)
path_16s<-"Data/Sequencing/16S/Data/"
path_meta<-"Data/Sequencing/metagenomic/Data/"

dir.create("Result/TableS5")
# path_med_16S<-"Result/mediation/16S"
# path_med_metagen<-"Result/mediation/metagenomic"
# path_med_func<-"Result/mediation/func"

# dir.create(path_med_16S)
# dir.create(path_med_metagen)
# dir.create(path_med_func)
data_total_16S<-read.csv(paste0(path_16s,"data_taxonomy_transformed.csv"),header = T,sep = ",")

data_total_meta<-read.csv(paste0(path_meta,"data_taxonomy_transformed.csv"),header = T,sep = ",")

data_total_func<-read.csv(paste0(path_meta,"data_func_transformed.csv"),header = T,sep = ",")

data_total_16S[pes_log_new]<-log(data_total_16S[pes_new])
data_total_meta[pes_log_new]<-log(data_total_meta[pes_new])
data_total_func[pes_log_new]<-log(data_total_func[pes_new])



# Table S4 ----------------------------------------------------------------

load("Data/microbialFeatures_standardized.RData")

# read meta data 
load("Data/metadata.RData")

load("Data/data_ers.RData")

dir.create("TableS4")

#1. amplicon 
rownames(amplicon.metadata)<-amplicon.metadata$SampleID
rownames(data_ers)<-data_ers$SampleID
meta <- cbind(amplicon.metadata,bact.features.st) 
meta<-cbind(meta,data_ers[,c("SampleID","ers")])

# Y: Glucose levels
Ys1 <- c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")
Ys2<-c("GDM","IGT")

#Y_df <- meta %>% select(SampleID, all_of(c(Ys1,Ys2)))

Y_df<-meta[,c("SampleID",Ys1,Ys2)]
# Exposure: 
Expos <- c("ers",pes_log_new)
#Expos_df <- meta %>% select(SampleID, all_of(Expos))
Expos_df<-meta[,c("SampleID",Expos)]
Expos_df <- cbind.data.frame("SampleID" = Expos_df$SampleID,
                             sapply(Expos_df[,-1], function(x) as.numeric(sub( "Y", 1, sub("N",0,x))) ),
                             stringsAsFactors=F) 


# df of covariables, rowname=sample, colnames=variable 
covariables <- covariate
rownames(meta)<-NULL
#covar <- meta %>% select(SampleID, all_of(covariables)) %>% tibble::column_to_rownames("SampleID") 
covar<-meta[,c("SampleID",covariables)]%>% tibble::column_to_rownames("SampleID") 


##GDM

Pvalue_res <- NULL

##linear regression
for(y in Ys2){
  # y=Ys[1]
  writeLines(paste0("glucose level:" , y, " , ", which(Ys2 == y), " out of ", length(Ys2)))
  for(ep in Expos){
    # ep = Expos[1]  
    writeLines(paste0("exposure: " , ep, " , ", which(Expos == ep), " out of ", length(Expos)))
    
    # meta.tmp <- merge(meta %>% select(SampleID, all_of(y)),
    #                   Expos_df %>% select(SampleID, all_of(ep)),
    #                   by="SampleID")
    
    meta.tmp <- merge(select(meta,all_of(c("SampleID",y))),
                      select(Expos_df,all_of(c("SampleID",ep))),
                      by="SampleID")
    
    for(df.name in c("bact.features.st")){
      # df.name = "bact.features.st"
      
      microbFtr_df <- eval(parse(text = df.name)) %>% as.data.frame()
      colnames(microbFtr_df) <- gsub("\\W", "_", colnames(microbFtr_df))
      writeLines(paste0("Microbial feature type: " , df.name))
      
      for(i in 1:ncol(microbFtr_df)){
        # i=1
        # writeLines(paste0(i, " out of ", ncol(microbFtr_df), " microbial features."))
        
        mf = colnames(microbFtr_df)[i]
        mf_df.tmp <- microbFtr_df %>% select(mf)
        
        dat <- merge(merge(meta.tmp, mf_df.tmp, by.x = "SampleID", by.y = 0),
                     covar,
                     by.x="SampleID",by.y=0)
        
        dat <- dat[complete.cases(dat),]
        
        fml <- paste0(y, " ~ ", paste(colnames(covar), collapse = " + "), " + ", ep, " + ", mf, " + ",
                      ep, " * ", mf )
        m <- glm(as.formula(fml), family = "binomial",data = dat)
        
        #m.res <- summary(m)$coefficients %>% as.data.frame()
        #m.res$`Pr(>|t|)`[i.pval]
        
        # an1 = anova(m)
        an1 = summary(m)
        # F.stat = an1[10,4]
        # 
        # Pvalue = an1[10,5]  
        Pvalue<-tryCatch({summary(m)$coefficients[10,4]},error=function(e){"NA"})
        OR<-tryCatch({exp(coef(m)[10])},error=function(e){"NA"})
        OR_l<-tryCatch({exp(confint(m)[10,1])},error=function(e){"NA"})
        OR_H<-tryCatch({exp(confint(m)[10,2])},error=function(e){"NA"})
        
        res_vec <- c("Glucose" = y, "Exposure" = ep,"MicrobType" = df.name, "MicrobFeature" = mf,
                     "OR.value" = OR,"OR_lower"=OR_l,"OR_higher"=OR_H,"P.interaction" = Pvalue)
        Pvalue_res <- bind_rows(Pvalue_res, res_vec)
        
      }
      
    } # loop through microbe features
  }# loop through exposures
}# loop through lung functions



dir.create("Result/TableS4")

write.csv(Pvalue_res, file = "Result/TableS4/Interaction_exposure.microbial_GDM.csv", quote = F, row.names = F)



#2. metagen 

rownames(metadata.meta)<-metadata.meta$SampleID
meta <- cbind(metadata.meta,metagTaxa.features.st) 


meta<-left_join(meta,data_ers[,c("number","ers")])

# Y: Glucose levels
Ys1 <- c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")
Ys2<-c("GDM","IGT")

Y_df <- meta %>% select(SampleID, all_of(c(Ys1,Ys2)))


# Exposure: 
Expos <- c("ers",pes_log_new)
Expos_df <- meta %>% select(SampleID, all_of(Expos))
Expos_df <- cbind.data.frame("SampleID" = Expos_df$SampleID,
                             sapply(Expos_df[,-1], function(x) as.numeric(sub( "Y", 1, sub("N",0,x))) ),
                             stringsAsFactors=F) 


# df of covariables, rowname=sample, colnames=variable 
covariables <- covariate
rownames(meta)<-NULL
covar <- meta %>% select(SampleID, all_of(covariables)) %>% tibble::column_to_rownames("SampleID") 



##GDM

Pvalue_res <- NULL

##linear regression
for(y in Ys2){
  # y=Ys[1]
  writeLines(paste0("glucose level:" , y, " , ", which(Ys2 == y), " out of ", length(Ys2)))
  for(ep in Expos){
    # ep = Expos[1]  
    writeLines(paste0("exposure: " , ep, " , ", which(Expos == ep), " out of ", length(Expos)))
    
    # meta.tmp <- merge(meta %>% select(SampleID, all_of(y)),
    #                   Expos_df %>% select(SampleID, all_of(ep)),
    #                   by="SampleID")
    
    meta.tmp <- merge(select(meta,all_of(c("SampleID",y))),
                      select(Expos_df,all_of(c("SampleID",ep))),
                      by="SampleID")
    
    for(df.name in c("metagTaxa.features.st","metagMod.features.st")){
      # df.name = "bact.features.st"
      
      microbFtr_df <- eval(parse(text = df.name)) %>% as.data.frame()
      colnames(microbFtr_df) <- gsub("\\W", "_", colnames(microbFtr_df))
      writeLines(paste0("Microbial feature type: " , df.name))
      
      for(i in 1:ncol(microbFtr_df)){
        # i=1
        # writeLines(paste0(i, " out of ", ncol(microbFtr_df), " microbial features."))
        
        mf = colnames(microbFtr_df)[i]
        mf_df.tmp <- microbFtr_df %>% select(mf)
        
        dat <- merge(merge(meta.tmp, mf_df.tmp, by.x = "SampleID", by.y = 0),
                     covar,
                     by.x="SampleID",by.y=0)
        
        dat <- dat[complete.cases(dat),]
        
        fml <- paste0(y, " ~ ", paste(colnames(covar), collapse = " + "), " + ", ep, " + ", mf, " + ",
                      ep, " * ", mf )
        m <- glm(as.formula(fml), family = "binomial",data = dat)
        
        #m.res <- summary(m)$coefficients %>% as.data.frame()
        #m.res$`Pr(>|t|)`[i.pval]
        
        # an1 = anova(m)
        an1 = summary(m)
        # F.stat = an1[10,4]
        # 
        # Pvalue = an1[10,5]  
        Pvalue<-tryCatch({summary(m)$coefficients[10,4]},error=function(e){"NA"})
        OR<-tryCatch({exp(coef(m)[10])},error=function(e){"NA"})
        OR_l<-tryCatch({exp(confint(m)[10,1])},error=function(e){"NA"})
        OR_H<-tryCatch({exp(confint(m)[10,2])},error=function(e){"NA"})
        
        res_vec <- c("Glucose" = y, "Exposure" = ep,"MicrobType" = df.name, "MicrobFeature" = mf,
                     "OR.value" = OR,"OR_lower"=OR_l,"OR_higher"=OR_H,"P.interaction" = Pvalue)
        Pvalue_res <- bind_rows(Pvalue_res, res_vec)
        
      }
      
    } # loop through microbe features
  }# loop through exposures
}# loop through lung functions

#write.csv(Pvalue_res, file = "Result/TableS4/Interaction_dimethoate.microbial_pathway_GDM_new.csv", quote = F, row.names = F)


#
#dir.create("Result/mediation/metagenomic/interaction")

write.csv(Pvalue_res, file = "Result/TableS4/Interaction_exposure.microbial_GDM_metagen.csv", quote = F, row.names = F)




# Table S5 ----------------------------------------------------------------
dir.create("Result/TableS5")

load("Data/microbialFeatures_standardized.RData")

# read meta data 
load("Data/metadata.RData")
load("Data/data_ers.RData")
####16s

library(mediation)
library(dplyr)
Mediation.forward <- NULL
Mediation.reverse <- NULL
Mediation.reverse2<- NULL
meta<-data_total_16S
for(exp in Exposures){
  #exp=Exposures[1]
  meta.expo <- meta %>% dplyr::select(X, all_of(exp))
  
  for(lf in c("GDM")){
    #lf = LungFuns[1]
    meta.lf = meta %>% dplyr::select(X, all_of(lf))
    
    for(mf.dn in c("bact.features.st")){
      #mf.dn = "bact.features.st"
      
      mf.data <- eval(parse(text = mf.dn))
      
      for(i in 1:ncol(mf.data)){
        mf.dat.tmp <- mf.data[,i, drop=F]
        mf = colnames(mf.data)[i]
        
        dat <- merge(merge(merge(meta.expo, meta.lf, by = "X"),
                           covar_df, by="X"),
                     mf.dat.tmp, by.x="X", by.y = 0)
        
        dat <- dat[complete.cases(dat),]
        
        # forward mediation --------------------
        id = paste(exp, mf, lf, sep = "|")
        
        model.m = lm(as.formula( paste(mf, " ~ ", exp, " + ", 
                                       paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                       sep = "") ), data = dat)
        
        model.y = glm(as.formula( paste(lf, " ~ ", exp," + ", mf," + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                        sep = "") ),family = "binomial", data = dat)
        
        summary = summary(mediate(model.m,model.y,treat=exp,mediator=mf,boot=F,sims=1000))
        res <- capture.output(summary,append=FALSE)
        
        #sub( "^()\\s", "\\1", res[7])
        # tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
        # ACME.p <- tmp[length(tmp)]
        ACME.p<-summary$d1.p
        # tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T) ]
        # ADE.p <- tmp[length(tmp)]
        ADE.p<-summary$z1.p
        
        
        # tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" ]
        # i_str = which(grepl("Mediated", tmp))
        # prop.mediated <- tmp[(i_str + 1)]
        prop.mediated<-summary$n.avg
        
        forw_vec = c(id, ACME.p, ADE.p, prop.mediated)
        names(forw_vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
        
        Mediation.forward <- bind_rows(Mediation.forward, forw_vec)
        
        # reverse mediation --------------------
        id = paste(exp, lf, mf, sep = "|")
        
        model.m = glm(as.formula( paste(lf, " ~ ", exp, " + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                        sep = "") ),family="binomial", data = dat)
        
        model.y = lm(as.formula( paste(mf, " ~ ", exp," + ", lf," + ", 
                                       paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                       sep = "") ), data = dat)
        
        summary = summary(mediate(model.m,model.y,treat=exp,mediator=lf,boot=F,sims=1000))
        res <- capture.output(summary,append=FALSE)
        
        #sub( "^()\\s", "\\1", res[7])
        # tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
        # ACME.p <- tmp[length(tmp)]
        ACME.p<-summary$d1.p
        # tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T) ]
        # ADE.p <- tmp[length(tmp)]
        ADE.p<-summary$z1.p
        
        
        # tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" ]
        # i_str = which(grepl("Mediated", tmp))
        # prop.mediated <- tmp[(i_str + 1)]
        prop.mediated<-summary$n.avg
        
        revers_vec = c(id, ACME.p, ADE.p, prop.mediated)
        names(revers_vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
        
        Mediation.reverse <- bind_rows(Mediation.reverse, revers_vec)
        
        # reverse mediation 2--------------------
        id = paste(mf,exp, lf,  sep = "|")
        
        model.m = lm(as.formula( paste(exp, " ~ ", mf, " + ", 
                                       paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                       sep = "") ), data = dat)
        
        model.y = glm(as.formula( paste(lf, " ~ ", mf," + ", exp," + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                        sep = "") ), family="binomial",data = dat)
        
        summary = summary(mediate(model.m,model.y,treat=mf,mediator=exp,boot=F,sims=1000))
        res <- capture.output(summary,append=FALSE)
        
        #sub( "^()\\s", "\\1", res[7])
        #sub( "^()\\s", "\\1", res[7])
        # tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
        # ACME.p <- tmp[length(tmp)]
        ACME.p<-summary$d1.p
        # tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T) ]
        # ADE.p <- tmp[length(tmp)]
        ADE.p<-summary$z1.p
        
        
        # tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" ]
        # i_str = which(grepl("Mediated", tmp))
        # prop.mediated <- tmp[(i_str + 1)]
        prop.mediated<-summary$n.avg
        
        revers_vec2 = c(id, ACME.p, ADE.p, prop.mediated)
        names(revers_vec2) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
        
        Mediation.reverse2 <- bind_rows(Mediation.reverse2, revers_vec2)
        
      }# loop through individual microbial features
    }# loop through microbial data frames
  }# loop through lung functions
}# loop through exposures

dir.create("Result/TableS5/16S")
write.csv(Mediation.forward,file="Result/TableS5/16S/Mediation.forward_GDM_with_pvalue.csv",row.names = F)

write.csv(Mediation.reverse,file="Result/TableS5/16S/Mediation.reverse_GDM_with_pvalue.csv",row.names = F)
write.csv(Mediation.reverse2,file="Result/TableS5/16S/Mediation.reverse2_GDM_with_pvalue.csv",row.names = F)




###metagenomics
Species<-read.csv(paste0(path_meta,"Classes\\Species_filtered.csv"),header=T,sep=",")
Species_name<-Species$Species
#tree<-Species[,c(1:7)]
Species_name_log<-paste0(Species_name,"_log")
dir.create(paste0(path_med_metagen,"/bidirection"))
covariate<-c("pre_BMI","Age","Educational_level","weightgain","Parity","passsmokHis_1y")

d1<-data_total_meta[,c(pes_new,covariate,"GDM")]

vars_toAnalyze<-pes_log_new


# exposures are the treators: 
Exposures <- c("ers",pes_log_new)

# outcome linear:
glu_levels <- c("OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")
outcome_binary<-c("GDM","IGT")

# Covariates:
Covars <- c("pre_BMI","Age","Educational_level","weightgain","Parity","passsmokHis_1y")
covar_df <- data_total_meta %>% dplyr::select(X, all_of(Covars)) #X means sampleID

metagTaxa.features.st<-data_total_meta[,Species_name_log]
rownames(metagTaxa.features.st)<-data_total_meta$X
# mediation analysis

# binomial levels 
library(mediation)
library(dplyr)
Mediation.forward <- NULL
Mediation.reverse <- NULL
Mediation.reverse2<- NULL
meta<-data_total_meta
for(exp in Exposures){
  #exp=Exposures[1]
  meta.expo <- meta %>% dplyr::select(X, all_of(exp))
  
  for(lf in c("GDM")){
    #lf = LungFuns[1]
    meta.lf = meta %>% dplyr::select(X, all_of(lf))
    
    for(mf.dn in c("metagTaxa.features.st")){
      #mf.dn = "bact.features.st"
      
      mf.data <- eval(parse(text = mf.dn))
      
      for(i in 1:ncol(mf.data)){
        mf.dat.tmp <- mf.data[,i, drop=F]
        mf = colnames(mf.data)[i]
        
        dat <- merge(merge(merge(meta.expo, meta.lf, by = "X"),
                           covar_df, by="X"),
                     mf.dat.tmp, by.x="X", by.y = 0)
        
        dat <- dat[complete.cases(dat),]
        
        # forward mediation --------------------
        id = paste(exp, mf, lf, sep = "|")
        
        model.m = lm(as.formula( paste(mf, " ~ ", exp, " + ", 
                                       paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                       sep = "") ), data = dat)
        
        model.y = glm(as.formula( paste(lf, " ~ ", exp," + ", mf," + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                        sep = "") ),family = "binomial", data = dat)
        
        summary = summary(mediate(model.m,model.y,treat=exp,mediator=mf,boot=F,sims=1000))
        res <- capture.output(summary,append=FALSE)
        
        # tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
        # ACME.p <- tmp[length(tmp)]
        ACME.p<-summary$d1.p
        # tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T) ]
        # ADE.p <- tmp[length(tmp)]
        ADE.p<-summary$z1.p
        
        
        # tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" ]
        # i_str = which(grepl("Mediated", tmp))
        # prop.mediated <- tmp[(i_str + 1)]
        prop.mediated<-summary$n.avg
        forw_vec = c(id, ACME.p, ADE.p, prop.mediated)
        names(forw_vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
        
        Mediation.forward <- bind_rows(Mediation.forward, forw_vec)
        
        # reverse mediation --------------------
        id = paste(exp, lf, mf, sep = "|")
        
        model.m = glm(as.formula( paste(lf, " ~ ", exp, " + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                        sep = "") ),family="binomial", data = dat)
        
        model.y = lm(as.formula( paste(mf, " ~ ", exp," + ", lf," + ", 
                                       paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                       sep = "") ), data = dat)
        
        summary = summary(mediate(model.m,model.y,treat=exp,mediator=lf,boot=F,sims=1000))
        res <- capture.output(summary,append=FALSE)
        
        #sub( "^()\\s", "\\1", res[7])
        # tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
        # ACME.p <- tmp[length(tmp)]
        ACME.p<-summary$d1.p
        # tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T) ]
        # ADE.p <- tmp[length(tmp)]
        ADE.p<-summary$z1.p
        
        
        # tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" ]
        # i_str = which(grepl("Mediated", tmp))
        # prop.mediated <- tmp[(i_str + 1)]
        prop.mediated<-summary$n.avg
        
        revers_vec = c(id, ACME.p, ADE.p, prop.mediated)
        names(revers_vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
        
        Mediation.reverse <- bind_rows(Mediation.reverse, revers_vec)
        
        # reverse mediation 2--------------------
        id = paste(mf,exp, lf,  sep = "|")
        
        model.m = lm(as.formula( paste(exp, " ~ ", mf, " + ", 
                                       paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                       sep = "") ), data = dat)
        
        model.y = glm(as.formula( paste(lf, " ~ ", mf," + ", exp," + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("X",lf, mf, exp))], collapse = " + "),
                                        sep = "") ), family="binomial",data = dat)
        
        summary = summary(mediate(model.m,model.y,treat=mf,mediator=exp,boot=F,sims=1000))
        res <- capture.output(summary,append=FALSE)
        
        #sub( "^()\\s", "\\1", res[7])
        # tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
        # ACME.p <- tmp[length(tmp)]
        ACME.p<-summary$d1.p
        # tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" & tmp!="."]
        # tmp <- tmp[!grepl("*",tmp,fixed = T) ]
        # ADE.p <- tmp[length(tmp)]
        ADE.p<-summary$z1.p
        
        
        # tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
        # tmp <- tmp[tmp != "" ]
        # i_str = which(grepl("Mediated", tmp))
        # prop.mediated <- tmp[(i_str + 1)]
        prop.mediated<-summary$n.avg
        
        revers_vec2 = c(id, ACME.p, ADE.p, prop.mediated)
        names(revers_vec2) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
        
        Mediation.reverse2 <- bind_rows(Mediation.reverse2, revers_vec2)
        
      }# loop through individual microbial features
    }# loop through microbial data frames
  }# loop through lung functions
}# loop through exposures
dir.create("Result/TableS5/metagen")
write.csv(Mediation.forward,file="Result/TableS5/metagen/Mediation.forward_GDM_with_pvalue.csv",row.names = F)

write.csv(Mediation.reverse,file="Result/TableS5/metagen/Mediation.reverse_GDM_with_pvalue.csv",row.names = F)
write.csv(Mediation.reverse2,file="Result/TableS5/metagen/Mediation.reverse2_GDM_with_pvalue.csv",row.names = F)

