
source_dir_root1<-"D:/Project/pes_GM_GDM/Data_review"
setwd(source_dir_root1)

# Figure 2a ----------------------------------------------------------------
# 导入所需的包
library(ape)
library(ggtree)
library(tidyverse)
library(openxlsx)
library(itol.toolkit)
# 读取物种丰度表
# abundance_table <- read.xlsx("taxo_abundance.xlsx", sheet = 1)
# rownames(abundance_table)<-abundance_table$label
# # 创建物种树对象
# tree <- read.tree("phylo.tre")

otu_filter<-read.table("Data/Sequencing/16S/1OTU/decontam/otu_mat_filter.txt",header = T,sep = "\t")
tax_filter<-read.table("Data/Sequencing/16S/1OTU/decontam/tax_mat_filter.txt",header = T,sep = "\t")
###
tree=read.tree("Data/Sequencing/16S/1OTU/decontam/rooted-tree-filter.nwk") # 读取nwk文件

# 1.taxonomy filter 

# otu_filter<-abundance_table[,-c(1:8)]
# tax_filter<-abundance_table[,1:8]

otutab<-otu_filter
###filter taxonomy
abundance = 0.01

norm = as.data.frame(t(t(otutab)/colSums(otutab,na=T)*100))
# 丰度由大到小排序
idx = order(rowMeans(norm), decreasing = T)
norm = norm[idx,]
# 按丰度筛选
idx = rowMeans(norm) > abundance
filtered_otutab = norm[idx,]

detectionRate=5

detect<-function(x){
  length(which(x > 0))*100/length(x)
}

DR<-apply(filtered_otutab, 1, detect)
idx3= apply(filtered_otutab, 1, detect) >detectionRate

filtered_otutab2<-filtered_otutab[idx3,]

# 按数量筛选
#filtered_otutab = head(norm, number)
# 添加均值并保留4位小数
filtered_otutab2 = round(cbind(rowMeans(filtered_otutab2), filtered_otutab2), digits = 4)
colnames(filtered_otutab2)[1] = "Mean"
# 对应过滤物种注释
idx = rownames(filtered_otutab2) %in% rownames(tax_filter)
filtered_otutab = filtered_otutab2[idx,]
filtered_taxonomy = tax_filter[rownames(filtered_otutab),]

##prune tree
data=fortify(tree)

pruned.tree<-drop.tip(tree,tree$tip.label[-match(filtered_taxonomy$Sequencing, tree$tip.label)])

#filtered_taxonomy$Phylum_new<-ifelse(filtered_taxonomy$Phylum%in%c("p__Campilobacterota","p__Desulfobacterota","p__Deferribacterota","p__Verrucomicrobiota"),"Others",filtered_taxonomy$Phylum)



#pruned.tree$tip.label<-filtered_taxonomy[match(pruned.tree$tip.label,filtered_taxonomy$label),]$X2

dir.create("Result/Figure2")
write.tree(pruned.tree,file = "Result/Figure2/Figure2a/pruned_tree.nwk")


write.csv(filtered_taxonomy,file = "Result/Figure2/Figure2a/filtered_taxonomy_16S.csv") ##手动添加宏基因组数据对应的平均物种丰度
write.csv(filtered_otutab,file = "Result/Figure2/Figure2a/filtered_otutab_16S.csv")


##metagenomic
Genus<-read.csv("Data/Sequencing/metagenomic/Data/Classes/Genus.csv",header = T,sep = ",")
rownames(Genus)<-Genus$Genus
Genus<-Genus[,-c(1:6)]

Genus = round(cbind(rowMeans(Genus), Genus), digits = 4)
colnames(Genus)[1] = "Mean"
write.csv(Genus,file = "Result/Figure2/Figure2a/genus_metagen.csv")


# tree plot---

###
tree2<-pruned.tree

hub <- create_hub(tree = tree2)

###调整参数
groups<-filtered_taxonomy[,c("Sequencing","Phylum","Genus")]
#count<-norm2[asv_id,]

count<-cbind(groups,filtered_otutab)
rownames(count)<-as.character(count$Sequencing)
count<-count[,-c(1:3)]
rownames(groups)<-groups$Sequencing

# relabel by genus，改名字
groups$Genus<-as.character(groups$Genus)
groups$Genus[groups$Genus=="Unassigned"]<-"g__Unassigned"
#rename genus name
groups2<-data.frame(groups,str_split_fixed(groups$Genus,"g__",2))

df_rename2 <- data.frame(id = groups2$Sequencing, 
                         new_label = groups2$X2)

unit_1 <- create_unit(data = df_rename2, 
                      key = "E3_rename", 
                      type = "LABELS",
                      tree = pruned.tree)

dir.create("Result/Figure2/Figure2a")
write_unit(unit_1,paste0("Result/Figure2","/Figure2a/itol_labels_Genus.txt"))




# tree_colors range by phylum 按照门填充分组颜色
unit_2 <- create_unit(data = filtered_taxonomy %>% select(Sequencing, Phylum),
                      key = "itol_3al_2_range",
                      type = "TREE_COLORS",
                      subtype = "range",
                      tree = pruned.tree)
write_unit(unit_2,paste0("Result/Figure2","/Figure2a/itol_range_Phylum.txt"))




# color_strip by class 按class分组圆环
set.seed(123)
unit_3 <- create_unit(data = filtered_taxonomy %>% select(Sequencing, Class),
                      key = "itol_3al_3_strip",
                      type = "DATASET_COLORSTRIP",
                      color = "wesanderson",
                      tree = pruned.tree)
unit_3@common_themes$basic_theme$margin <- 50
write_unit(unit_3,paste0("Result/Figure2","/Figure2a/itol_strip_Class.txt"))


dat_bar<-read.xlsx("Result/Figure2/Figure2a/tree_abundance.xlsx",sheet=1,rowNames = T)
names(dat_bar)<-c("label","Genus","amplicon","metagen","metagen_presence")
filtered_taxonomy2<-cbind(filtered_taxonomy,dat_bar)
# simple_bar by NS 柱状图


## full 

# filtered_taxonomy2$pre<-apply(filtered_taxonomy2[,16:19],1,sum)
# filtered_taxonomy2$mid<-apply(filtered_taxonomy2[,12:15],1,sum)
# MULTIBAR_bar by OS NS，多个数据bar
unit_5 <- create_unit(data = filtered_taxonomy2 %>% select(Sequencing, amplicon,metagen),
                      key = "itol_3al_5_multibar",
                      type = "DATASET_MULTIBAR",
                      tree = pruned.tree)
unit_5@specific_themes$basic_plot$size_max <- 300
write_unit(unit_5,paste0("Result/Figure2","/Figure2a/itol_multibar.txt"))

filtered_taxonomy2$amplicon_presence<-rep(1,times=nrow(filtered_taxonomy2))
##绘制热图

unit_7 <- create_unit(data = filtered_taxonomy2 %>% select(Sequencing,amplicon_presence,metagen_presence), 
                      key = "itol_7_heatmap",
                      type = "DATASET_HEATMAP", 
                      tree = pruned.tree)
write_unit(unit_7,paste0("Result/Figure2","/Figure2a/itol_heatmap.txt"))


# color_strip by class 按class分组圆环
set.seed(123)
unit_8 <- create_unit(data = filtered_taxonomy2 %>% select(Sequencing, metagen_presence),
                      key = "itol_metagen",
                      type = "DATASET_COLORSTRIP",
                      color = "table2itol",
                      tree = pruned.tree)
write_unit(unit_8,paste0("Result/Figure2","/Figure2a/itol_metagen.txt"))

# color_strip by class 按class分组圆环
set.seed(123)
unit_9 <- create_unit(data = filtered_taxonomy2 %>% select(Sequencing, amplicon_presence),
                      key = "itol_metagen",
                      type = "DATASET_COLORSTRIP",
                      color = "wesanderson",
                      tree = pruned.tree)
write_unit(unit_9,paste0("Result/Figure2","/Figure2a/itol_amplicon.txt"))

#使用Fig2a文件夹中的数据，itol网站中复现：https://itol.embl.de/

# Figure 2b ---------------------------------------------------------------
path_meta<-"Data/Sequencing/metagenomic/Data/"
path_16s<-"Data/Sequencing/16S/Data/"
### import data 
library(dplyr)
metadata = read.csv("Result/linear/Concensus_hc/pes_sample_cluster.csv", header = T,sep = ",")
#metadata<-subset(metadata,X%in%ID_new)
#ID_new<-paste0("trim.",metadata$X,".fq")
rownames(metadata)<-metadata$X
dat_cluster<-metadata[,c("Clust")]

dat_cluster2<-metadata[,c("X","Clust")]
metadata_meta<-read.csv(paste0(path_meta,"metadata.csv"),header = T,sep = ",",row.names = 1)


data_16S<-read.csv(paste0(path_16s,"data_taxonomy_transformed.csv"),header = T,sep = ",",row.names = 1)
data_16S<-cbind(data_16S,dat_cluster2)
data_16S[pes_log_new]<-log(data_16S[pes_new])
data_meta<-read.csv(paste0(path_meta,"data_taxonomy_transformed.csv"),header = T,sep = ",",row.names = 1)

data_meta$X<-paste0("trim.",data_meta$number,".fq")
data_meta<-left_join(data_meta,dat_cluster2, by = c("X" = "X"))
data_meta[pes_log_new]<-log(data_meta[pes_new])
dat_func<-read.csv(paste0(path_meta,"data_func_transformed.csv"),header = T,sep = ",",row.names = 1)

dat_func$X<-paste0("trim.",dat_func$number,".fq")
dat_func<-left_join(dat_func,dat_cluster2, by = c("X" = "X"))
dat_func[pes_log_new]<-log(dat_func[pes_new])

#dat_func<-left_join(dat_func,dat_cluster2, by = c("number" = "X"))

library(egg)
covariate<-c("pre_BMI","Age","Clust","GDM","Educational_level","weightgain","Parity",
             "passsmokHis_1y",pes_log_new
)

confounder.anova= data_16S[covariate]
#confounder.anova$age = as.numeric(confounder.anova$age)
#confounder.anova$bmi = as.numeric(confounder.anova$bmi)

names = covariate
names


library(RColorBrewer)
pbrbl <- colorRampPalette(brewer.pal(9, "Blues"), interpolate = "spline")

sample_ID<-as.character(data_16S$X)
adonis.res = list()
method= "unifrac"
# run PERMANOVA test on weighted uniFrac distance

library(openxlsx)
source_dir<-"D:/Project/pes_GM_GDM"
setwd(source_dir)
#dir.create("Result/pes_GM")
#setwd("D:/Project/pes_GM_GDM/Sequencing/metagenomic")


def_permoanova<-function(j){
  data<-read.table(paste0(path_16s,j,"_distance_matrix.tsv"),header=TRUE,row.names=1)
  
  data<-data[sample_ID,sample_ID]
  adonis.res = list()
  for(i in 1:length(names)){
    # confounder2 = confounder.anova[confounder.anova$iMSMS_ID %in% rownames(weight),]
    # confounder2 = confounder2[!is.na(confounder2[,i+1]) & confounder2[,i+1] != "", ]
    dis<-data
    adonis.res[[i]] = vegan::adonis(as.formula(paste("dis ~",names[i], sep = "")), data = confounder.anova)
    #   if(method =="bray"){
    #     abun = weight[rownames(weight) %in% confounder2$iMSMS_ID, ]
    #     dis = vegdist(abun, method = "bray",diag =T, upper =T)
    #   }else{
    #     dis = weight[match(confounder2$iMSMS_ID, rownames(weight)), match(confounder2$iMSMS_ID, colnames(weight))]
    #   }
    #   adonis.res[[i]] = vegan::adonis(as.formula(paste("dis ~",names[i], sep = "")), data = confounder2)
  }
  names(adonis.res) = c(names)
  # extract the R2 and Pvalue
  result = matrix(NA, nrow = length(names), ncol =2)
  for(i in 1:(length(names))){
    result[i,1] = adonis.res[[i]]$aov.tab$R2[1]
    result[i,2] = adonis.res[[i]]$aov.tab$`Pr(>F)`[1]
  }
  rownames(result) = c(names)
  colnames(result) = c("R2", "Pvalue")
  result = data.frame(result, stringsAsFactors = F)
  result$Padjust = p.adjust(result$Pvalue, method = "fdr")
  
  result$ID = rownames(result)
  # for(i in 1:nrow(result)){
  #   result$Group[i] = termgroup[termgroup$Term == result$ID[i], "Group"]
  # }
  subresult<-result[result$Pvalue < 0.1,]
  subresult<-subresult[order(subresult$R2,decreasing = T),]
  name_<-as.character(subresult$ID)
  presult = result[result$Pvalue < 0.05,]
  padj.result = result[result$Padjust < 0.05,]
  #  anova.cols = c("Demography" = "#0571B0", "Disease"="#CCCCCC", "Life style" ="#92C5DE", "Medication"="#F4A582", "Physiology"="#CA0020")
  # p<-ggplot(subresult, aes(x = reorder(ID, R2),y=R2, fill = ID)) +
  subresult$ID<-factor(subresult$ID,levels =rev(name_))
  p<-ggplot(subresult, aes(x = ID,y=R2, fill = ID)) +
    geom_bar(stat='identity') +
    #   scale_x_discrete(limits=name_)+
    coord_flip() + ylab("Adonis R2") + xlab("") +
    #    scale_fill_brewer()+
    scale_fill_manual(values = pbrbl(nrow(subresult)))+
    #    scale_fill_manual(values = anova.cols) +
    #    scale_fill_bre
    geom_text(data = presult, aes(ID, R2),label="*", col= "black",nudge_y = 0.0005, nudge_x = -0.15)+
    geom_text(data = padj.result, aes(ID, R2),label="*", col= "red",nudge_y = 0.0005,nudge_x = -0.15)+
    theme_article()
  ggsave(p,filename = paste0(path_16s_result,"/Diversity/permo2/",j,"_barplot.pdf"),width=8,height=6)
  write.csv(result,file = paste0(path_16s_result,"/Diversity/permo2/",j,"_permoanova.csv"))
}

outcome<-c("unweighted_unifrac","weighted_unifrac","bray_curtis")
dir.create(paste0(path_16s_result,"/Diversity/permo2"))
for (j in outcome) {
  def_permoanova(j)
}



# Figure 2c-e -------------------------------------------------------------

metadata = read.csv("Data/data_imputed.csv", header = T,sep = ",")
ID_new<-paste0("trim.",metadata$number,".fq")

rownames(metadata)<-ID_new
path_meta_result<-"Result/GM_GDM/metagenomic"

package_list <- c("vegan","dplyr","readr","cluster","tidyverse","clusterSim")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# dir.create("Result/GM_GDM")
# dir.create("Result/GM_GDM/metagenomic")
# dir.create("Result/GM_GDM/metagenomic/Diversity")

# PAM

# change the metaphlan result to a composition table, select top n most abundant features
CompositionTable <- function(x,n) { 
  require(foreach)
  #  x[is.na(x)]=0
  mean_value <- data.frame(Taxa=colnames(x), Mean_abundance=colSums(x)/nrow(x))
  most <- as.character(mean_value[order(mean_value$Mean_abundance,decreasing = T),]$Taxa[1:n])
  print(paste("Most abundant taxa is",most,sep = " "))
  
  composition_table <- foreach(i=1:length(most),.combine = rbind) %do%  {
    return.string = data.frame(ID = rownames(x), Relative=x[,most[i]],Level=colnames(x[,most[i],drop=F]))
  }
  
  first <- composition_table[grep(most[1],composition_table$Level),]
  first <- first[order(first$Relative,decreasing = T),]
  level <- as.factor(first$ID)
  composition_table$ID <- factor(composition_table$ID,levels = level)
  
  return(composition_table)
}

# calculator for JSD distance
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

# partition around medoid (PAM) clustering
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
}


source_dir<-"D:/project/pes_GM_GDM/Data_review"
setwd(source_dir)
dat_meta<-read.csv("Data/Sequencing/metagenomic/Data/Classes/Species_filtered.csv", header = T,sep = ",")
#table_with_taxonomy<-table_with_taxonomy[,grep("s__",colnames(table_with_taxonomy))]
rownames(dat_meta)<-paste0("seq",seq(1,nrow(dat_meta)))

tax_mat<- dat_meta[,1:7]

#rownames(table_with_taxonomy)<-rownames(tax_mat)
otu_mat<- dat_meta[,-c(1:7)]

rownames(otu_mat)<-dat_meta$Species

#transform to relative abundance
#dag3_species <- otu_mat/100.0
dag3_species <- otu_mat

dag3_species<-as.data.frame(t(dag3_species))

colnames(dag3_species)=lapply(colnames(dag3_species),function(x){
  strsplit(x,"s__")[[1]][2]
})

dag3_species_plot=CompositionTable(dag3_species,10)

dag3_species_plot$Level=factor(dag3_species_plot$Level,levels = c("Anaerobutyricum_hallii","Faecalibacterium_prausnitzii",
                                                                  "Blautia_wexlerae",
                                                                  "Candidatus_Cibiobacter_qucibialis",
                                                                  "Escherichia_coli",
                                                                  "Anaerostipes_hadrus",
                                                                  "Bifidobacterium_pseudocatenulatum",
                                                                  "Clostridia_bacterium",
                                                                  "Adlercreutzia_equolifaciens",
                                                                  "Streptococcus_salivarius"))
dag3_species_plot[dag3_species_plot==0]=NA
dag3_species_plot=na.omit(dag3_species_plot)
dag3_species_plot$Relative=-log2(dag3_species_plot$Relative)
set.seed(10)
n <- 10
library(randomcoloR)
library(ggridges)
palette <- distinctColorPalette(n)

g <- ggplot(dag3_species_plot, aes(x = Relative, y = Level,fill=Level,color=Level)) + theme_classic() +
  geom_density_ridges(alpha=0.8,scale = 2)+scale_fill_manual(values = palette) +
  scale_color_manual(values = palette)+ylab("")+xlab("Norm. Rel. Abundance (-log2)")

print(g)
dir.create("Result/Figure2")
ggsave(plot = g, filename = 'Result/Figure2/Figure2c.png')
#ggsave(plot = g, filename = 'Result/microbiome_clustering/cluster.pdf')

# species clustering using PAM
dag3_species=as.data.frame(t(dag3_species))
dag3_species.dist=dist.JSD(dag3_species)
dag3_species.cluster=pam.clustering(dag3_species.dist, k=3)

# > perform clustering
nclusters = index.G1(t(dag3_species), dag3_species.cluster, d = dag3_species.dist, centrotypes = "medoids")
nclusters=NULL

# > calculate CH index to identify optimal number of clusters
for (k in 1:10) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    dag3_species.cluster_temp=pam.clustering(dag3_species.dist, k)
    nclusters[k]=index.G1(t(dag3_species),dag3_species.cluster_temp,  d = dag3_species.dist,
                          centrotypes = "medoids")
  }
}
png(file = "Result/Figure2/Figure2d.png")
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
dev.off()
#ggsave(p, filename = 'microbiome_clustering/optimal_cluster_n.pdf')
cluster=data.frame(row.names = colnames(dag3_species),Cluster=dag3_species.cluster)
write.table(cluster,"Result/Figure2/DMP_species_PAM.txt",row.names = T,quote = F,sep = "\t")
write.table(nclusters,"Result/Figure2/ncluster.txt",row.names = T,quote = F,sep = "\t")
dag3_species=as.data.frame(t(dag3_species))


# P.copri abundance per cluster

clusters=read.table("Result/Figure2/DMP_species_PAM.txt",sep = "\t",header = T,stringsAsFactors = F)
clusters$sampleID=rownames(clusters)

#cluster_plot=merge(clusters,dag3_species[,"Prevotella_copri",drop=F],by.x="sampleID",by.y = "row.names",all=F)
cluster_plot=merge(clusters,dag3_species[,"Prevotella_stercorea",drop=F],by.x="sampleID",by.y = "row.names",all=F)

cluster_plot$Cluster=factor(cluster_plot$Cluster,levels = c("3","2","1"))
cluster_plot[cluster_plot==0]=NA
cluster_plot=na.omit(cluster_plot)
cluster_plot$Prevotella_stercorea=-log2(cluster_plot$Prevotella_stercorea)
cluster_plot[is.na(cluster_plot)]=0

##other plot
data.cluster=dag3_species.cluster
data.dist=dag3_species.dist
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])


#remove noise
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}
data=as.data.frame(t(dag3_species))
data.denoized=noise.removal(data, percent=0.01)

library(ade4)
pdf(file = "Result/Figure2/enterotyping.pdf")
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F,sub="Between-class analysis",col=c(4,2,3))


s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(4,2,3))

obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), sub="Principal coordiante analysis",grid=F,col=c(4,2,3))

s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4))


dev.off()



# TableS3----------------------------------------------------------------
setwd("D:/Project/pes_GM_GDM/Data_review")
load("Data/microbialFeatures_standardized.RData")

# read meta data 
load("Data/metadata.RData")

load("Data/data_ers.RData")


GDM_id<-subset(metadata.meta,GDM==1)$SampleID
control_id<-subset(metadata.meta,GDM==0)$SampleID

abund_GDM<-metagTaxa.features.st[GDM_id,]
abund_GDM<-as.data.frame(t(abund_GDM))
#abund_GDM$x<-rownames(abund_GDM)


abund_control<-metagTaxa.features.st[control_id,]
abund_control<-as.data.frame(t(abund_control))
#abund_control$x<-rownames(abund_control)

library(NetMoss2)
library(Hmisc)
library(dplyr)

cr<-t(abund_control)
net_control<-rcorr(cr,type="spearman")
net_control<-net_control$r


cr<-t(abund_GDM)
net_GDM<-rcorr(cr,type="spearman")
net_GDM<-net_GDM$r

# NetMoss 

nodes_result_healthy.vs.GDM = 
  NetMoss(case_dir = abund_GDM,
          control_dir = abund_control,
          net_case_dir = net_GDM,
          net_control_dir = net_control)


NMSS_1.2 = nodes_result_healthy.vs.GDM[[1]]



dir.create("Result/TableS3")
write.table(NMSS_1.2, file = "Result/TableS3/NetMoss_metagen.txt", sep = "\t", quote = F, row.names = F)  #Table S3


# Figure4 -----------------------------------------------------------------

dir.create("Result/Figure4")
dir.create("Result/Figure4/d")

netPlot2 <- function(result,
                     num.top = 5,
                     num.score = 30,
                     e.th = 0.4,
                     nodeSize = 2,
                     nodeDize = 1.2,
                     edgeSize = 0.2,
                     edgeDize = 0.3,
                     arrowSize = 0,
                     my.layout = layout_as_star,
                     my.label = TRUE)
{
  ############################################
  my.wd = getwd()
  nodes_result = result[[1]]
  nodes_result = nodes_result[order(nodes_result$NetMoss_Score,decreasing = T),]
  top.tax = as.character(nodes_result[1:num.top,1])
  
  ##abundance fold change
  case.all.sample = result[[4]]
  control.all.sample = result[[5]]
  case.all.sample = case.all.sample[order(as.character(case.all.sample$X)),]
  control.all.sample = control.all.sample[order(as.character(control.all.sample$X)),]
  all.sample.data = data.frame()
  for (k in 1:nrow(case.all.sample))
  {
    all.sample.data[k,1] = case.all.sample[k,1]
    all.sample.data[k,2] = mean(as.numeric(as.character(case.all.sample[k,-1])))
    all.sample.data[k,3] = mean(as.numeric(as.character(control.all.sample[k,-1])))
    all.sample.data[k,4] = log2((all.sample.data[k,2]+1)/(all.sample.data[k,3]+1))
  }
  colnames(all.sample.data) = c("genus","case","control","FC")
  rownames(all.sample.data) = all.sample.data$genus
  
  ####plot NetMoss score
  nodes_plot = nodes_result[1:num.score,]
  nodes_plot$taxon_names = factor(nodes_plot$taxon_names,levels = rev(nodes_plot$taxon_names))
  nodes_plot$FC = all.sample.data[as.character(nodes_plot$taxon_names),"FC"]
  nodes_plot[which(nodes_plot$FC < 0),"type"] = "Enriched in control"
  nodes_plot[which(nodes_plot$FC >= 0),"type"] = "Enriched in case"
  
  ##point
  p1 = ggplot(nodes_plot,aes(taxon_names,NetMoss_Score, color = p.adj))+
    scale_colour_distiller(palette = "Spectral")+
    geom_bar(stat = 'identity',width = 0.03, fill = "black")+
    geom_point(size = 4)+
    coord_flip()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    labs(y = "NetMoss score", x = paste0("Top ",num.score," bacteria"))+
    theme(legend.position = "top")
  
  #fold change
  p2 = ggplot(nodes_plot,aes(taxon_names,FC,color = type))+
    geom_point(size = 4)+coord_flip()+
    theme_bw()+theme(panel.grid = element_blank())+
    geom_hline(yintercept = 0)+
    scale_color_brewer(palette = "Pastel1", name = "")+
    labs(x = "", y = "Log2(FC)")+
    theme(legend.position = "top")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  pp = ggpubr::ggarrange(p1,p2,NULL,ncol = 3,nrow = 1,widths = c(1,0.3,0.1),align = "h")
  
  
  #####highlight top taxon
  ##case
  e1.case = data.frame(result[[2]])
  e1.case[upper.tri(e1.case)] <- 0
  e1.case$genus = rownames(e1.case)
  e2.case = melt(e1.case, id.vars = "genus")
  e2.case = e2.case[which(e2.case$value != 0), ]
  ##
  e2.case$genus = as.character(e2.case$genus)
  e2.case$genus = gsub("/.", "-", e2.case$genus)
  e2.case$genus = gsub(" ", "-", e2.case$genus)
  e2.case$genus = gsub('/[', '', e2.case$genus)
  e2.case$genus = gsub('/]', '', e2.case$genus)
  e2.case$genus = gsub("/.", "-", e2.case$genus)
  e2.case$variable = as.character(e2.case$variable)
  e2.case$variable = gsub("/.", "-", e2.case$variable)
  e2.case$variable = gsub(" ", "-", e2.case$variable)
  e2.case$variable = gsub('/[', '', e2.case$variable)
  e2.case$variable = gsub('/]', '', e2.case$variable)
  e2.case$variable = gsub("/.", "-", e2.case$variable)
  
  edge.case = e2.case[which(e2.case$genus != e2.case$variable), ]
  edge.case2 = edge.case[which(abs(edge.case$value) > e.th), ]
  
  node.case = data.frame(unique(c(
    as.character(edge.case2$genus),
    as.character(edge.case2$variable)
  )))
  colnames(node.case) = "node"
  rownames(node.case) = node.case$node
  rownames(nodes_result) = nodes_result$taxon_names
  
  top.tax2 = intersect(top.tax,node.case$node)
  if (length(top.tax2) == 0)
  {
    warning("the threshold is too high to find the target taxon in the case networks!")
  }
  
  inter.genus = intersect(as.character(node.case$node),as.character(nodes_result$taxon_names))
  node.case[inter.genus,"score"] = nodes_result[inter.genus,"NetMoss_Score"]
  node.case2 = node.case
  if (ncol(node.case2[which(is.na(node.case2$score)),]) != 0)
  {
    node.case2[which(is.na(node.case2$score)),"score"] = 0
  }else
  {
    node.case2 = node.case2
  }
  
  node.case2$weight = abs(node.case2$score)
  node.case2$class = 1   ### gray
  node.case2[top.tax2, "class"] = 2  ###red
  
  edge.case2$weight = 0.1
  for (kk in 1:length(top.tax2))
  {
    for (mm in 1:nrow(edge.case2))
    {
      if (as.character(edge.case2[mm,1]) == top.tax2[kk] || as.character(edge.case2[mm,2]) == top.tax2[kk])
      {
        edge.case2[mm,"weight"] = 1
      }else
      {
        next
      }
    }
  }
  
  edge.case2$class = 1  ##positive red
  edge.case2[which(edge.case2$value < 0), "class"] = 2   #negative blue
  
  
  ###plot
  g1 <- graph.empty()
  g1 <- graph_from_data_frame(edge.case2, vertices = node.case2)
  #
  nodeSize <- nodeSize
  nodeDize <- nodeDize
  edgeSize <- edgeSize
  edgeDize <- edgeDize
  VColor <- c("#8f8f8f", "#e57265")
  EColor <- c("#e57265", "#78a0c4")
  VText <- c(0.2,1)
  V(g1)$size <-
    nodeSize + nodeDize * 10 * as.numeric(as.vector(node.case2$weight))
  V(g1)$color <- VColor[node.case2$class]
  V(g1)$label.cex <- VText[node.case2$class]
  V(g1)$frame.color <- "black"
  E(g1)$width <-
    edgeSize + (edgeDize * abs(3 * as.numeric(as.vector(
      edge.case2$weight
    ))))
  E(g1)$color <- EColor[edge.case2$class]
  E(g1)$arrow.size <- arrowSize
  
  
  ##control
  e1.control = data.frame(result[[3]])
  e1.control[upper.tri(e1.control)] <- 0
  e1.control$genus = rownames(e1.control)
  e2.control = melt(e1.control, id.vars = "genus")
  e2.control = e2.control[which(e2.control$value != 0), ]
  ##
  e2.control$genus = as.character(e2.control$genus)
  e2.control$genus = gsub("/.", "-", e2.control$genus)
  e2.control$genus = gsub(" ", "-", e2.control$genus)
  e2.control$genus = gsub('/[', '', e2.control$genus)
  e2.control$genus = gsub('/]', '', e2.control$genus)
  e2.control$genus = gsub("/.", "-", e2.control$genus)
  e2.control$variable = as.character(e2.control$variable)
  e2.control$variable = gsub("/.", "-", e2.control$variable)
  e2.control$variable = gsub(" ", "-", e2.control$variable)
  e2.control$variable = gsub('/[', '', e2.control$variable)
  e2.control$variable = gsub('/]', '', e2.control$variable)
  e2.control$variable = gsub("/.", "-", e2.control$variable)
  
  edge.control = e2.control[which(e2.control$genus != e2.control$variable), ]
  edge.control2 = edge.control[which(abs(edge.control$value) > e.th), ]
  
  node.control = data.frame(unique(c(
    as.character(edge.control2$genus),
    as.character(edge.control2$variable)
  )))
  colnames(node.control) = "node"
  rownames(node.control) = node.control$node
  rownames(nodes_result) = nodes_result$taxon_names
  
  top.tax2 = intersect(top.tax,node.control$node)
  if (length(top.tax2) == 0)
  {
    warning("the threshold is too high to find the target taxon in the control networks!")
  }
  
  inter.genus = intersect(as.character(node.control$node),as.character(nodes_result$taxon_names))
  node.control[inter.genus,"score"] = nodes_result[inter.genus,"NetMoss_Score"]
  node.control2 = node.control
  if (ncol(node.control2[which(is.na(node.control2$score)),]) != 0)
  {
    node.control2[which(is.na(node.control2$score)),"score"] = 0
  }else
  {
    node.control2 = node.control2
  }
  
  node.control2$weight = abs(node.control2$score)
  node.control2$class = 1   ### gray
  node.control2[top.tax2, "class"] = 2  ###red
  
  edge.control2$weight = 0.1
  for (kk in 1:length(top.tax2))
  {
    for (mm in 1:nrow(edge.control2))
    {
      if (as.character(edge.control2[mm,1]) == top.tax2[kk] || as.character(edge.control2[mm,2]) == top.tax2[kk])
      {
        edge.control2[mm,"weight"] = 1
      }else
      {
        next
      }
    }
  }
  
  edge.control2$class = 1  ##positive red
  edge.control2[which(edge.control2$value < 0), "class"] = 2   #negative blue
  
  
  ###plot
  g2 <- graph.empty()
  g2 <- graph_from_data_frame(edge.control2, vertices = node.control2)
  #
  nodeSize <- nodeSize
  nodeDize <- nodeDize
  edgeSize <- edgeSize
  edgeDize <- edgeDize
  VColor <- c("#8f8f8f", "#e57265")
  EColor <- c("#e57265", "#78a0c4")
  VText <- c(0.2,1)
  V(g2)$size <-
    nodeSize + nodeDize * 10 * as.numeric(as.vector(node.control2$weight))
  V(g2)$color <- VColor[node.control2$class]
  V(g2)$label.cex <- VText[node.control2$class]
  V(g2)$frame.color <- "black"
  E(g2)$width <-
    edgeSize + (edgeDize * abs(3 * as.numeric(as.vector(
      edge.control2$weight
    ))))
  E(g2)$color <- EColor[edge.control2$class]
  E(g2)$arrow.size <- arrowSize
  
  ########
  if (!my.label)
  {
    par(mfrow = c(1, 1))
    ##half circle
    p11 = ggplot() +
      geom_ellipse(aes(x0 = 0, y0 = 0, a = 10, b = 6, angle = 0),
                   fill = "#f6a5c0", alpha = 0.1, color = NA) +
      scale_y_continuous(limits = c(-6,0),expand = c(0,0))+coord_flip()+
      theme_void()
    ##half circle
    p12 = ggplot() +
      geom_ellipse(aes(x0 = 0, y0 = 0, a = 10, b = 6, angle = 0),
                   fill = "#93cf96", alpha = 0.1, color = NA) +
      scale_y_continuous(limits = c(0,6),expand = c(0,0))+coord_flip()+
      theme_void()
    
    multiplot(p11,p12,cols = 2)
    
    ##case
    par(fig=c(0,0.5,0.1,0.9),new=TRUE)
    plot(
      g1,
      layout = my.layout,
      vertex.frame.color = "NA",
      vertex.shapes = "none",
      vertex.label = ""
    )
    
    ##title
    par(fig=c(0.3,0.4,0.8,0.95),new=TRUE)
    title("case network")
    
    ##legend
    par(fig=c(0.1,0.4,0.1,0.3),new=TRUE)
    legend(
      x = -1,
      y = -1.5,
      bty = "n",
      ####edges color
      c("positive correlation", "negative correlation"),
      lty = 1,
      lwd = 2,
      col = EColor
    )
    
    ##control
    par(fig=c(0.5,1,0.1,0.9),new=TRUE)
    plot(
      g2,
      layout = my.layout,
      vertex.frame.color = "NA",
      vertex.shapes = "none",
      vertex.label = ""
    )
    
    ##title
    par(fig=c(0.5,0.7,0.8,0.95),new=TRUE)
    title("control network")
    
    ##legend
    par(fig=c(0.5,0.8,0.1,0.3),new=TRUE)
    legend(
      x = -1.5,
      y = -1.5,
      bty = "n",
      ####nodes color
      c("others", "key taxon"),
      pch = 21,
      pt.bg = VColor
    )
  }else
  {
    par(mfrow = c(1, 1))
    ##half circle
    p11 = ggplot() +
      geom_ellipse(aes(x0 = 0, y0 = 0, a = 10, b = 6, angle = 0),
                   fill = "#f6a5c0", alpha = 0.1, color = NA) +
      scale_y_continuous(limits = c(-6,0),expand = c(0,0))+coord_flip()+
      theme_void()
    ##half circle
    p12 = ggplot() +
      geom_ellipse(aes(x0 = 0, y0 = 0, a = 10, b = 6, angle = 0),
                   fill = "#93cf96", alpha = 0.1, color = NA) +
      scale_y_continuous(limits = c(0,6),expand = c(0,0))+coord_flip()+
      theme_void()
    
    #multiplot(p11,p12,cols = 2)
    multiplot(p11,p12,ncol  = 2)
    
    ##case
    par(fig=c(0,0.5,0.1,0.9),new=TRUE)
    plot(
      g1,
      layout = my.layout,
      vertex.frame.color = "NA",
      vertex.shapes = "none"
    )
    
    ##title
    par(fig=c(0.3,0.4,0.8,0.95),new=TRUE)
    title("case network")
    
    ##legend
    par(fig=c(0.1,0.4,0.1,0.3),new=TRUE)
    legend(
      x = -1,
      y = -1.5,
      bty = "n",
      ####edges color
      c("positive correlation", "negative correlation"),
      lty = 1,
      lwd = 2,
      col = EColor
    )
    
    ##control
    par(fig=c(0.5,1,0.1,0.9),new=TRUE)
    plot(
      g2,
      layout = my.layout,
      vertex.frame.color = "NA",
      vertex.shapes = "none"
    )
    
    ##title
    par(fig=c(0.5,0.7,0.8,0.95),new=TRUE)
    title("control network")
    
    ##legend
    par(fig=c(0.5,0.8,0.1,0.3),new=TRUE)
    legend(
      x = -1.5,
      y = -1.5,
      bty = "n",
      ####nodes color
      c("others", "key taxon"),
      pch = 21,
      pt.bg = VColor
    )
  }
  
  setwd(my.wd)
  
  ###score
  ggsave("NetMoss_score.pdf", pp)
  
  print (paste0("the NetMoss score saved in ", my.wd))
}
setwd("D:/Project/pes_GM_GDM/Data_review/Result/Figure4")
dir.create("d")
dir.create("a")
dir.create("b")

###figure 4d
setwd("D:/Project/pes_GM_GDM/Data_review/Result/Figure4/d")
#plot networks
netPlot2(result = nodes_result_healthy.vs.GDM,
         num.top = 5,
         num.score = 30,
         e.th = 0.5,
         my.layout = layout_as_star,
         my.label = TRUE)

#ggsave(filename = "Network_plot_metagen.pdf",height = 8,width = 12)


##Figure 4a
setwd(source_dir_root1)

rownames(metadata.meta)<-metadata.meta$SampleID
meta<-cbind(metadata.meta,metagTaxa.features.st)
meta$Disease<-ifelse(meta$GDM==1,"GDM","Healthy control")

combined.data<-meta[,c(colnames(metagTaxa.features.st),"Disease")]
Row.names<-rownames(combined.data)

combined.data<-data.frame(Row.names,combined.data)
rownames(combined.data)<-NULL
# calculate networks 
library(dplyr)
library(data.table)
library(Hmisc)
library(tidyverse)

edges_3Diseases <- NULL
for(dss in unique(meta$Disease)){
  
  dat.tmp <- combined.data %>% 
    filter(Disease == dss) %>% 
    tibble::column_to_rownames("Row.names") %>%
    dplyr::select(-Disease)
  
  Corr <- rcorr(as.matrix(dat.tmp) , type="spearman")
  occor.r <- Corr$r
  occor.p <- Corr$P
  
  # Hide lower triangle for r (so each correlation is only going to appear once)
  lower <- occor.r
  lower[lower.tri(occor.r)]<-NA
  
  # edges
  r_df.l <- lower %>% reshape2::melt(value.name = "r")
  p_df.l <- occor.p %>% reshape2::melt(value.name = "p")
  
  if(all(r_df.l$Var1 == p_df.l$Var1) & all(r_df.l$Var2 == p_df.l$Var2)){
    occor.rp_l <- cbind.data.frame(r_df.l, p=p_df.l$p, stringsAsFactors=F) # quicker
  }else{
    occor.rp_l <- merge(occor.r %>% reshape2::melt(value.name = "r"), 
                        occor.p %>% reshape2::melt(value.name = "p"),
                        by=c("Var1", "Var2"))  # too slowwwwwwwww
  }
  
  
  edges <- occor.rp_l %>% 
    filter(!is.na(r)) %>% 
    filter(abs(r) < 1) %>% 
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) 
  
  edges$edge = paste(edges$Var1, edges$Var2, sep = "|")
  edges$disease = dss
  edges$edgeType <-
    paste0(substr(edges$edge,1,1), 
           substr(sapply(strsplit(edges$edge, "|", fixed = T), "[[", 2), 1, 1))
  
  edges_3Diseases <- bind_rows(edges_3Diseases, edges)
  
}

edges_3Diseases$r <- as.numeric(edges_3Diseases$r)
edges_3Diseases$p <- as.numeric(edges_3Diseases$p)

# edges_3Diseases$source <- sapply(edges_3Diseases$Var1,
#                                  function(x) microbeAbb$V2[which(microbeAbb$abb == x)])
# edges_3Diseases$target <- sapply(edges_3Diseases$Var2,
#                                  function(x) microbeAbb$V2[which(microbeAbb$abb == x)])

edges_3Diseases$source<-edges_3Diseases$Var1
edges_3Diseases$target<-edges_3Diseases$Var2


for(dss in c( "GDM","Healthy control")){
  
  edges <- edges_3Diseases %>% 
    filter(disease == dss) %>%
    mutate(abs.r =  abs(r)) %>%
    filter(abs.r > 0.4) %>%
    filter(p < 0.05) %>%
    relocate(source, target) 
  
  write.csv(edges, file = paste0("Result/Figure4/a/",dss,"_edges.csv"),quote = F, row.names = F)
}


##GDM 

library(openxlsx)
dss = "GDM" # manually change: Health, preCOPD, COPD
edges <- fread(paste0("Result/Figure4/a/",dss,"_edges.csv"), data.table = F) 
#nodeColors <- eval(parse(text = paste0("nodesColors_",dss)))
nodeColors<-read.xlsx("Result/Figure4/nodeColors_GDM.xlsx", sheet = 1) 

#microbeAbb<-read.xlsx("metadata/Taxa_16S.xlsx",sheet = 1)
microbeAbb<-read.csv("Data/Sequencing/metagenomic/Data/Classes/Species.csv",header=T,sep=",")
microbeAbb<-microbeAbb[,1:7]
microbeAbb$phylum<-gsub("p__","",microbeAbb$Phylum)

internal.node<-unique(c(edges$source,edges$target))

extraNodes<-setdiff(nodeColors$Id,internal.node)

nodeColors<-subset(nodeColors,Id%in%internal.node)

#extraNodes <- microbeAbb$V2[!microbeAbb$V2 %in% as.character(nodeColors$Id )]

selfEdges <- cbind.data.frame(source = extraNodes, 
                              target = extraNodes,
                              edgeType = "selfLink",
                              stringsAsFactors=F) 

edges <- bind_rows(edges, selfEdges)

nodes <- 
  data.frame(table(c(edges$source, edges$target)) ) %>%
  mutate(id = Var1) %>% dplyr::select(-Var1) %>% 
  mutate(abb = sapply(id, function(x) microbeAbb$Genus[which(microbeAbb$Species == x)])) %>%
  mutate(type = sub("_/d+","", abb)) %>%
  mutate(phylum =  sapply(id, function(x) microbeAbb$phylum[which(microbeAbb$Species == x)])) %>%
  mutate(phylum6 = sapply(phylum, 
                          function(x){
                            if(x %in% c("Bacteroidetes","Firmicutes","Proteobacteria","Actinobacteria")) x else "Others"
                          } )) %>%
  relocate(id)  



# plot with igraph
library(igraph)
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
net


V(net)$module.color <- 
  sapply(names(V(net)),
         function(x){
           if(x %in% nodeColors$Id){
             nodeColors$color[which(nodeColors$Id == x)]
           }else "gray"
         }) 

df <- cbind.data.frame(
  V(net)$module.color,
  V(net)$phylum6,
  V(net)$type
)
assign(paste0(dss,"_df"), df, envir = .GlobalEnv)

V(net)$type.color <-
  sapply(1:length(V(net)),
         function(i){
           if(V(net)$type[i] == "bacteria") "white" else V(net)$module.color[i]
         })


E(net)$source.module.color  <- 
  sapply(E(net)$Var1,
         function(x) {
           if(is.na(x) | !x %in% nodeColors$Id) "gray" else nodeColors$color[which(nodeColors$Id == x)]
         })


net.simp <- net - E(net)[E(net)$edgeType=="selfLink"]  #remove self links, only keep the nodes

# Figure4a:
pdf(paste("Result/Figure4/a/", dss,".pdf",sep = ""), 
    width = 4 , height = 4) 
#储存的图片大小（对应圆的面积）跟node个数成正比，所以长宽与node个数的平方根成正比
par(mar=c(0.1,0.1,0.1,0.1)) 
set.seed(100)

plot(net.simp, vertex.label=NA, 
     vertex.size = 4,
     vertex.color = V(net)$type.color,
     vertex.frame.color= V(net)$module.color,
     edge.size = 1, 
     edge.color = E(net)$source.module.color, 
     layout=layout_with_fr(net)
)

dev.off()


## Control -----------------------------------------------------------------


dss = "Healthy control" # manually change: Health, preCOPD, COPD
edges <- fread(paste0("Result/Figure4/a/",dss,"_edges.csv"), data.table = F) 
#nodeColors <- eval(parse(text = paste0("nodesColors_",dss)))
nodeColors<-read.xlsx(paste0("Result//Figure4/","nodeColors_control.xlsx"), sheet = 1) 
microbeAbb<-microbeAbb[,1:7]
microbeAbb$phylum<-gsub("p__","",microbeAbb$Phylum)

internal.node<-unique(c(edges$source,edges$target))

extraNodes<-setdiff(nodeColors$Id,internal.node)

nodeColors<-subset(nodeColors,Id%in%internal.node)

#extraNodes <- microbeAbb$V2[!microbeAbb$V2 %in% as.character(nodeColors$Id )]

selfEdges <- cbind.data.frame(source = extraNodes, 
                              target = extraNodes,
                              edgeType = "selfLink",
                              stringsAsFactors=F) 

edges <- bind_rows(edges, selfEdges)

nodes <- 
  data.frame(table(c(edges$source, edges$target)) ) %>%
  mutate(id = Var1) %>% dplyr::select(-Var1) %>% 
  mutate(abb = sapply(id, function(x) microbeAbb$Species[which(microbeAbb$Species == x)])) %>%
  mutate(type = sub("_/d+","", abb)) %>%
  mutate(phylum =  sapply(id, function(x) microbeAbb$phylum[which(microbeAbb$Species == x)])) %>%
  mutate(phylum6 = sapply(phylum, 
                          function(x){
                            if(x %in% c("Bacteroidetes","Firmicutes","Proteobacteria","Actinobacteria")) x else "Others"
                          } )) %>%
  relocate(id)  



# plot with igraph 
library(igraph)
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
net


V(net)$module.color <- 
  sapply(names(V(net)),
         function(x){
           if(x %in% nodeColors$Id){
             nodeColors$color[which(nodeColors$Id == x)]
           }else "gray"
         }) 

df <- cbind.data.frame(
  V(net)$module.color,
  V(net)$phylum6,
  V(net)$type
)
assign(paste0(dss,"_df"), df, envir = .GlobalEnv)

V(net)$type.color <-
  sapply(1:length(V(net)),
         function(i){
           if(V(net)$type[i] == "bacteria") "white" else V(net)$module.color[i]
         })


E(net)$source.module.color  <- 
  sapply(E(net)$Var1,
         function(x) {
           if(is.na(x) | !x %in% nodeColors$Id) "gray" else nodeColors$color[which(nodeColors$Id == x)]
         })


net.simp <- net - E(net)[E(net)$edgeType=="selfLink"]  #remove self links, only keep the nodes

#dir.create("Result/GM_GDM/metagenomic/Network/NetPlots")
# Fig4a:
pdf(paste("Result/Figure4/a/", dss,".pdf",sep = ""), 
    width = 4 , height = 4) 
#储存的图片大小（对应圆的面积）跟node个数成正比，所以长宽与node个数的平方根成正比
par(mar=c(0.1,0.1,0.1,0.1)) 
set.seed(100)

plot(net.simp, vertex.label=NA, 
     vertex.size = 4,
     vertex.color = V(net)$type.color,
     vertex.frame.color= V(net)$module.color,
     edge.size = 1, 
     edge.color = E(net)$source.module.color, 
     layout=layout_with_fr(net)
)

dev.off()





# Figure 4b ---------------------------------------------------------------

tmp <- 
  rbind.data.frame(
    `Healthy control_df` %>% mutate(Disease = "Healthy control"),
    GDM_df %>% mutate(Disease = "GDM"),
    stringsAsFactors = F
  ) %>% 
  mutate(module.color=`V(net)$module.color`,
         phylum = `V(net)$phylum6`,
         type = `V(net)$type`) %>%
  mutate(module =  sapply(module.color,
                          function(x){
                            if(x %in% c("#43997A", "#47A265")){
                              "green"
                            }else if(x %in% c("#7F5477", "#91569B")){
                              "purple"
                            }else if(x == "#5D6795"){
                              "blue"
                            }else if(x == "#CB6651"){
                              "brown"
                            }else if(x == "#E41A1C"){
                              "red"
                            }else if(x == "#FFD422"){
                              "yellow"
                            }else x
                          }))

plotDat.phylumPie <- tmp %>%
  group_by(module, phylum) %>%
  dplyr::summarise(n=n()) %>%
  #calculate relative proportions
  group_by(module) %>%
  mutate(perc=n/sum(n))

library(ggplot2)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

phylumPie <- ggplot(plotDat.phylumPie %>% filter(module != "gray"), aes(x="", y=perc, fill=phylum))+
  facet_grid(.~module) + 
  geom_bar(width = 1, stat = "identity", color="white") +
  coord_polar("y", start=0) +
  blank_theme +
  theme(axis.text.x=element_blank())+
  scale_fill_manual(values = c('#84AADA','#65A484',"#D2BDDB","#DAD7D7","#E1E0AB"))

ggsave(phylumPie,file=paste("Result/Figure4/b/pieplot.pdf",sep = ""),width = 12,height = 4)


# plot module size 
plotDat.moduleSize <- tmp %>%
  filter(module != "gray") %>%
  group_by(Disease, module) %>%
  dplyr::summarise(n=n())


plotDat.moduleSize$Disease <- factor(plotDat.moduleSize$Disease, levels = c("Low","High"))
moduleSize <- ggplot(plotDat.moduleSize ) +
  geom_point(aes(x=Disease, y=module, size=n, color=module), shape = 21) +
  scale_color_manual(values = c("#5D6795","#CB6651","#43997A","pink2","#7F5477","#FFD422")) +
  scale_size(range = c(2, 12)) +
  theme_bw() + theme(panel.grid = element_blank())+
  labs(x="")
ggsave(moduleSize,file="Result/Figure4/b/moduleSize.pdf",width = 6,height = 6)



# plot module type percentage 
plotDat.moduleTypePerc <- tmp %>%
  filter(module != "gray") %>%
  group_by(Disease, module, type) %>%
  dplyr::summarise(n=n()) %>%
  group_by(Disease, module) %>%
  mutate(perc = n/sum(n)) 

#type<-rep("bacteria",times=nrow(plotDat.moduleTypePerc))

#plotDat.moduleTypePerc<-data.frame(plotDat.moduleTypePerc,type)
plotDat.moduleTypePerc$fillType <- 
  sapply(1:nrow(plotDat.moduleTypePerc),
         function(i){
           if(plotDat.moduleTypePerc$type[i] == "bacteria"){
             "white"
           }else{
             if(plotDat.moduleTypePerc$module[i] == "blue"){
               "#5D6795"
             }else if(plotDat.moduleTypePerc$module[i] == "brown"){
               "#CB6651"
             }else if(plotDat.moduleTypePerc$module[i] == "green"){
               "#43997A"
             }else if(plotDat.moduleTypePerc$module[i] == "pink2"){
               "pink2"
             }else if(plotDat.moduleTypePerc$module[i] == "purple"){
               "#7F5477"
             }else if(plotDat.moduleTypePerc$module[i] == "red"){
               "#E41A1C"
             }else if(plotDat.moduleTypePerc$module[i] == "yellow"){
               "#FFD422"
             }
           }
         })

plotDat.moduleTypePerc$Disease <- factor(plotDat.moduleTypePerc$Disease, levels = c("Low","High"))
Fig4b.moduleTypePerc <- ggplot(plotDat.moduleTypePerc, aes(x="", y=perc, fill=fillType)) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  facet_grid(module ~ Disease) +
  scale_fill_manual(values =  unique(plotDat.moduleTypePerc$fillType)[order(unique(plotDat.moduleTypePerc$fillType))])+
  #scale_color_manual(values = unique(plotDat.moduleTypePerc$module)[order(unique(plotDat.moduleTypePerc$module))]) +
  theme_void() # remove background, grid, numeric labels
ggsave(Fig4b.moduleTypePerc,file="Result/Figure4/b/moduleTypePerc.pdf",width = 6,height = 6)


# Figure 4c ---------------------------------------------------------------
setwd("D:/Project/pes_GM_GDM/Data_review")
load("Data/microbialFeatures_standardized.RData")
# read meta data 
load("Data/metadata.RData")
load("Data/data_ers.RData")

GDM_id<-subset(amplicon.metadata,GDM==1)$SampleID
control_id<-subset(amplicon.metadata,GDM==0)$SampleID

abund_GDM<-bact.features.st[GDM_id,-c(1:2)]
abund_GDM<-as.data.frame(t(abund_GDM))
#abund_GDM$x<-rownames(abund_GDM)


abund_control<-bact.features.st[control_id,-c(1:2)]
abund_control<-as.data.frame(t(abund_control))
#abund_control$x<-rownames(abund_control)

library(NetMoss2)
library(Hmisc)
library(dplyr)

cr<-t(abund_control)
net_control<-rcorr(cr,type="spearman")
net_control<-net_control$r


cr<-t(abund_GDM)
net_GDM<-rcorr(cr,type="spearman")
net_GDM<-net_GDM$r

nodes_result_healthy.vs.GDM = 
  NetMoss(case_dir = abund_GDM,
          control_dir = abund_control,
          net_case_dir = net_GDM,
          net_control_dir = net_control)


result = nodes_result_healthy.vs.GDM[[1]]

dir.create("Result/Figure4/c")
write.table(result, file = "Result/Figure4/c/NetMoss.txt", sep = "\t", quote = F, row.names = F)

detach("egg")
library(ggpubr)
#run netplot2
setwd("D:/Project/pes_GM_GDM/Data_review/Result/Figure4/c")
#plot networks
netPlot2(result = nodes_result_healthy.vs.GDM,
         num.top = 5,
         num.score = 30,
         e.th = 0.4,
         my.layout = layout_as_star,
         my.label = TRUE)


# Figure 5 interaction term -----------------------------------------------
setwd(source_dir_root1)
dir.create("Result/Figure5")
rownames(amplicon.metadata)<-amplicon.metadata$SampleID
rownames(data_ers)<-data_ers$SampleID
meta <- cbind(amplicon.metadata,bact.features.st) 
meta<-cbind(meta,data_ers[,c("SampleID","ers")])

Prevotella_median<-median(meta$g__Prevotella)
meta$Group<-ifelse(meta$g__Prevotella>Prevotella_median,"Prevotella_High","Prevotella_Low")

meta$Group<-factor(meta$Group,levels=c("Prevotella_Low","Prevotella_High"))

y<-"GDM"
ep<-"ers"
mf<-"Group"
dat<-meta
covar<-covariate
df.name<-"bact.features.st"
fml <- paste0(y, " ~ ", paste(colnames(covar), collapse = " + "), " + ", ep, " + ", mf, " + ",
              ep, " * ", mf )
m <- glm(as.formula(fml), family = "binomial",data = dat)
an1 = summary(m)
Pvalue<-tryCatch({summary(m)$coefficients[10,4]},error=function(e){"NA"})
OR<-tryCatch({exp(stats::coef(m)[10])},error=function(e){"NA"})
OR_l<-tryCatch({exp(confint(m)[10,1])},error=function(e){"NA"})
OR_H<-tryCatch({exp(confint(m)[10,2])},error=function(e){"NA"})

res_vec <- c("Glucose" = y, "Exposure" = ep,"MicrobType" = df.name, "MicrobFeature" = mf,
             "OR.value" = OR,"OR_lower"=OR_l,"OR_higher"=OR_H,"P.interaction" = Pvalue)
#Pvalue_res <- bind_rows(Pvalue_res, res_vec)

#dir.create("Result/mediation/16S/interaction")

write.csv(res_vec, file = "Result/Figure5/Interaction_ers.Prevotella_GDM.csv", quote = F, row.names = F)  ##Figure 5c

library(interactions)

p<-interact_plot(m,
                 pred = ers,  # 自变量
                 interval = TRUE,
                 #plot.points = TRUE,
                 modx = Group)    # 调节变量
ggsave(p, file = "Result/Figure5/Interaction_ers.Prevotella_plot.pdf") #Figure 5a



####stratified analysis
meta1<-subset(meta,Group=="Prevotella_Low")
meta2<-subset(meta,Group=="Prevotella_High")

m1<-glm(GDM~ers+Age+pre_BMI+Parity+Educational_level+weightgain+passsmokHis_1y,family = "binomial",data=meta1)
m2<-glm(GDM~ers+Age+pre_BMI+Parity+Educational_level+weightgain+passsmokHis_1y,family = "binomial",data=meta2)

Pvalue<-tryCatch({summary(m1)$coefficients[2,4]},error=function(e){"NA"})
OR<-tryCatch({exp(stats::coef(m1)[2])},error=function(e){"NA"})
OR_l<-tryCatch({exp(confint(m1)[2,1])},error=function(e){"NA"})
OR_H<-tryCatch({exp(confint(m1)[2,2])},error=function(e){"NA"})

dat_plot1<-c(OR,OR_l,OR_H,Pvalue)

Pvalue<-tryCatch({summary(m2)$coefficients[2,4]},error=function(e){"NA"})
OR<-tryCatch({exp(stats::coef(m2)[2])},error=function(e){"NA"})
OR_l<-tryCatch({exp(confint(m2)[2,1])},error=function(e){"NA"})
OR_H<-tryCatch({exp(confint(m2)[2,2])},error=function(e){"NA"})

dat_plot2<-c(OR,OR_l,OR_H,Pvalue)

dat_plot<-t(data.frame(dat_plot1,dat_plot2))

colnames(dat_plot)<-c("OR","Lower","Upper","Pvalue")
Group<-c("Prevotella_Low","Prevotella_High")
dat_plot<-data.frame(dat_plot,Group)

write.csv(dat_plot,file = "Result/Figure5/Interaction_ers.Prevotella_GDM.csv")
###errorbar

p<-ggplot(data = dat_plot,aes(x = Group,y = OR)) +
  geom_point(shape=18,size=3) +
  geom_errorbar(aes(ymin = Lower,ymax = Upper),width=0.2) +
  theme_bw()+
  geom_hline(yintercept = 1,color="red",linetype="dotted")+
  labs(x="",y="GDM OR (95%CI)")
ggsave(p, file = "Result/Figure5/Interaction_ers.Prevotella_errorbar.pdf",width = 4,height = 6) #Figure 5b


##pcoa plot
library(ggplot2)
library(RColorBrewer )
# pick a nice color palette
palette(brewer.pal(n = 8, name = "Set2"))
myPalette <- c("#F1DE33", "#B5D236","#90C549","#56B664","#2BA67B","#229683","#268388","#2A7189",
               "#325D87","#3A4281","#442A75","#3B2955")

median_ers<-median(meta$ers)
meta$group_ers<-ifelse(meta$ers<median_ers,"L","H")

dat_plot<-meta[,c("PCOA1","PCOA2","g__Prevotella","Group","group_ers")]

p0<- ggplot() +
  geom_point(data = dat_plot,aes(x=PCOA1, y=PCOA2, color=g__Prevotella), size=2.8, alpha=0.7) +
  stat_ellipse(data = dat_plot %>% filter(group_ers == "H"),
               aes(x=PCOA1, y=PCOA2, group = Group,linetype=Group),
               color="#A81F24", level = 0.8, lwd = 1) +
  stat_ellipse(data = dat_plot %>% filter(group_ers == "L"),
               aes(x=PCOA1, y=PCOA2, group = Group,linetype = Group),
               color="darkgray",level = 0.8, lwd = 1) +
  scale_color_gradientn(colours = myPalette) +
  theme_bw() #+ theme(panel.grid = element_blank(),
# axis.title = element_blank(),
# axis.text = element_blank(),
# axis.ticks = element_blank())


dat_plot$group_ers<-as.factor(dat_plot$group_ers)
p1 <- ggplot(data = dat_plot) +
  geom_density(aes(x = PCOA1, group = Group, linetype = Group, color=Group),lwd = 1) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     legend.position = "none",
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank()) +
  scale_color_manual(values = c("darkgray","#A81F24")) 

p2 <- ggplot(data = dat_plot) +
  geom_density(aes(x = PCOA2, group = Group, linetype = Group, color=Group),lwd = 1) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     legend.position = "none",
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank()) +
  scale_color_manual(values = c("darkgray","#A81F24")) +
  coord_flip()


library(ggpubr)

p<-ggarrange(
  ggarrange(p1, ggplot(), widths = c(0.8,0.2)),
  ggarrange(p0 + theme(legend.position =  c(0.85, 0.5)), p2, widths = c(0.8,0.2)),
  ncol = 1,  heights = c(0.2,0.8)
)
ggsave(p,file = "Result/Figure5/Interaction_ers.Prevotella_PCOA_plot.pdf",width = 8,height = 8)  ##Figure 5c




# Figure 6 correlationplot and mediation ----------------------------------

library(data.table)
library(dplyr)
library(foreach)
setwd(source_dir_root1)
dir.create("Result/Figure6")
# Figure 6a ---------------------------------------------------------------

load("Data/microbialFeatures_standardized.RData")

# read meta data 
load("Data/metadata.RData")
load("Data/data_ers.RData")
#rownames(amplicon.metadata)<-amplicon.metadata$SampleID

amplicon.metadata<-merge(amplicon.metadata,data_ers[,c("number","ers")],by="number")
meta <- cbind(amplicon.metadata,bact.features.st)
vars_toAnalyze <-c(pes_log_new,"ers","OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24")
all(vars_toAnalyze %in% colnames(meta))

# df of covariables, rowname=sample, colnames=variable 
covariables <- c("Age","pre_BMI","passsmokHis_1y","Educational_level","Parity","weightgain")

rownames(meta)<-NULL
covar <- meta %>% dplyr::select(SampleID, all_of(covariables)) %>% tibble::column_to_rownames("SampleID") 
#covar$Gender <- as.integer(sub( "F", 1, sub("M",0,covar$Gender)))

# df of phenotypes excluding covariables,  rowname=sample, colnames=phenotypes 
pheno <- meta %>% 
  dplyr::select(SampleID, all_of(vars_toAnalyze)) %>% 
  tibble::column_to_rownames("SampleID")  
# pheno[pheno == "N"] <- 0
# pheno[pheno == "Y"] <- 1
colnames(pheno) <- gsub("/W","_", colnames(pheno))

# lm to calculate association between phenotypes and microbial features 
# matrix of features：rowname=sample, colnames=features

for(df.name in c("bact.features.st")){
  # df.name = "metagMod.features.st" 
  ftrs_transformed <- eval(parse(text = df.name))
  # belowing are codes from DMP paper: https://github.com/GRONINGEN-MICROBIOME-CENTRE/DMP 
  # Run multivariate models, multi-thread implementation
  # prep parallelization
  writeLines(paste0("parallel calculation of lm started for:  ",df.name))
  # registerDoSEQ() # debug mode = single threaded
  threads = 4
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  result_ftrs = foreach(i = 1:ncol(pheno),.combine = rbind) %:%
   foreach(j = 1:ncol(ftrs_transformed),.combine = rbind) %dopar% {  # parallel implementation
    predictors = data.frame(covar[!is.na(pheno[,i]),],
                            model.matrix(
                              as.formula(paste0("~ ",colnames(pheno)[i])),data = pheno)[,-1,drop = F])
    
    cleaned_data = predictors[complete.cases(predictors),]
    rn <- rownames(cleaned_data)
    rn <- rn[rn %in% rownames(ftrs_transformed)]
    ftrs.cleaned = ftrs_transformed[rn,]
    cleaned_data = cleaned_data[rn,]
    if (nrow(cleaned_data) > 3) {
      if(T){
        # make model
        s1 = lm(
          as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)))),
          data = cleaned_data
        )
        # debug: print model
        s0 = lm(
          as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[1:ncol(covar)]))),
          data = cleaned_data
        )
        
        # compare models
        an1 = anova(s1,s0)
        output = data.frame(
          phenotype = colnames(pheno)[i],
          taxon = colnames(ftrs.cleaned)[j],
          Nsamples = nrow(cleaned_data),
          levels = if(class(pheno[,i]) == "factor") paste(collapse=":",levels(pheno[,i])) else "Not Applicable",
          levels_SampleSize = 
            if(class(pheno[,i]) == "factor" | length(table(pheno[,i]))==2) paste(collapse= ":",table(pheno[,i])) else "Not Applicable",
          effect.size =
            if(class(pheno[,i]) == "factor") {
              paste(collapse = ":",c(0,round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))])))
            } else round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))]) ,
          
          R2 = summary(s1)$r.squared - summary(s0)$r.squared,
          F.stat = an1[2,5],
          Pvalue = an1[2,6]
        )
      }
      
      
      #add covariates
      output
    }# condition: nrow(cleaned_data)>3 etc
  }# parallel inside
  writeLines("ftrs_done")
  rownames(result_ftrs) <- NULL
  on.exit(stopCluster(cl))
  # return results
  result_ftrs
  # export results
  write.table(result_ftrs, file = paste0("Result/Figure6/lmRes_",df.name,"_withCovariates.txt"), quote = F, row.names = F, sep = "\t" )
  writeLines(paste0(df.name, " has finished and results exported."))
}

#####metagen
rownames(metadata.meta)<-metadata.meta$SampleID
meta <- cbind(metadata.meta,metagMod.features.st,metagTaxa.features.st)

##combine ers and metagenomic data
dat_ers<-data_ers[,c("number","ers")]

rownames(metadata.meta)<-metadata.meta$SampleID
meta <- cbind(metadata.meta,metagMod.features.st,metagTaxa.features.st)

meta<-left_join(meta,dat_ers[,c("number","ers")],by="number")

# SampleID<-rownames(meta)
# meta<-data.frame(SampleID,meta)
vars_toAnalyze <-c(pes_log_new,"OGTT0_24","OGTT1_24","OGTT2_24","HbAlc_24","ers","GDM")
all(vars_toAnalyze %in% colnames(meta))

# df of covariables, rowname=sample, colnames=variable 
covariables <- c("Age","pre_BMI","passsmokHis_1y","Educational_level","Parity","weightgain")

rownames(meta)<-NULL
covar <- meta %>% dplyr::select(SampleID, all_of(covariables)) %>% tibble::column_to_rownames("SampleID") 
#covar$Gender <- as.integer(sub( "F", 1, sub("M",0,covar$Gender)))

# df of phenotypes excluding covariables,  rowname=sample, colnames=phenotypes 
pheno <- meta %>% 
  dplyr::select(SampleID, all_of(vars_toAnalyze)) %>% 
  tibble::column_to_rownames("SampleID")  
# pheno[pheno == "N"] <- 0
# pheno[pheno == "Y"] <- 1
colnames(pheno) <- gsub("/W","_", colnames(pheno))


for(df.name in c("metagMod.features.st", "metagTaxa.features.st")){
  # df.name = "metagMod.features.st" 
  ftrs_transformed <- eval(parse(text = df.name))
  # belowing are codes from DMP paper: https://github.com/GRONINGEN-MICROBIOME-CENTRE/DMP 
  # Run multivariate models, multi-thread implementation
  # prep parallelization
  writeLines(paste0("parallel calculation of lm started for:  ",df.name))
  # registerDoSEQ() # debug mode = single threaded
  threads = 4
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  # debug: timer
  # t1 <- Sys.time()
  # loop over all phenotypes
  result_ftrs = foreach(i = 1:ncol(pheno),.combine = rbind) %:%
  foreach(j = 1:ncol(ftrs_transformed),.combine = rbind) %dopar% {  # parallel implementation
    predictors = data.frame(covar[!is.na(pheno[,i]),],
                            model.matrix(
                              as.formula(paste0("~ ",colnames(pheno)[i])),data = pheno)[,-1,drop = F])
    
    cleaned_data = predictors[complete.cases(predictors),]
    rn <- rownames(cleaned_data)
    rn <- rn[rn %in% rownames(ftrs_transformed)]
    ftrs.cleaned = ftrs_transformed[rn,]
    cleaned_data = cleaned_data[rn,]
    if (nrow(cleaned_data) > 3) {
      # debug: print model
      #print(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data))))
      # lm as applied in DMP scripts: use lm because Y is features, always continuous variable with normal distribution
      # if covariates are considered
      if(T){
        # make model
        s1 = lm(
          as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)))),
          data = cleaned_data
        )
        # debug: print model
        #print(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[1:ncol(covar)])))
        
        # make model with extra covariates
        s0 = lm(
          as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[1:ncol(covar)]))),
          data = cleaned_data
        )
        
        # compare models
        an1 = anova(s1,s0)
        output = data.frame(
          phenotype = colnames(pheno)[i],
          taxon = colnames(ftrs.cleaned)[j],
          Nsamples = nrow(cleaned_data),
          levels = if(class(pheno[,i]) == "factor") paste(collapse=":",levels(pheno[,i])) else "Not Applicable",
          levels_SampleSize = 
            if(class(pheno[,i]) == "factor" | length(table(pheno[,i]))==2) paste(collapse= ":",table(pheno[,i])) else "Not Applicable",
          effect.size =
            if(class(pheno[,i]) == "factor") {
              paste(collapse = ":",c(0,round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))])))
            } else round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))]) ,
          
          R2 = summary(s1)$r.squared - summary(s0)$r.squared,
          F.stat = an1[2,5],
          Pvalue = an1[2,6]
        )
      }
      #add covariates
      output
    }# condition: nrow(cleaned_data)>3 etc
  }# parallel inside
  # debug
  #t2 <-  Sys.time()
  writeLines("ftrs_done")
  # debug
  #print(t2-t1)
  rownames(result_ftrs) <- NULL
  on.exit(stopCluster(cl))
  # return results
  result_ftrs
  # export results
  write.table(result_ftrs, file = paste0("Result/Figure6/lmRes_",df.name,"_withCovariates.txt"), quote = F, row.names = F, sep = "\t" )
  writeLines(paste0(df.name, " has finished and results exported."))
}

# loop through data types

# padj for each results data frame 
res <-fread("Result/Figure6/lmRes_bact.features.st_withCovariates.txt")
test<-
  bind_rows( res %>%
               filter(!taxon %in% c("PCOA1","PCOA2")) %>%
               group_by(phenotype) %>%
               mutate(padj = p.adjust(Pvalue)),
             res %>% filter(taxon %in% c("PCOA1","PCOA2"))) %>%
  dplyr::select(-levels,-levels_SampleSize) %>%
  arrange(Pvalue)
write.table(test, file = "Result/Figure6/lmRes_bact.features.st_withCovariates_padj.txt", quote = F, row.names = F, sep = "\t")

res <-fread("Result/Figure6/lmRes_metagMod.features.st_withCovariates.txt")
test<-
  bind_rows( res %>%
               #filter(!taxon %in% c("ko.alpha.st","MetagMod.pco1.st")) %>%
               group_by(phenotype) %>%
               mutate(padj = p.adjust(Pvalue)),
             res #%>% filter(taxon %in%  c("ko.alpha.st","MetagMod.pco1.st"))
  ) %>%
  dplyr::select(-levels,-levels_SampleSize) %>%
  arrange(Pvalue)
write.table(test, file = "Result/Figure6/lmRes_metagMod.features.st_withCovariates_padj.txt", quote = F, row.names = F, sep = "\t")

res <-fread("Result/Figure6/lmRes_metagTaxa.features.st_withCovariates.txt")
test<-
  bind_rows( res %>%
               #  filter(!taxon %in% c("metagTaxa.alpha.st","metagTaxa.pco1")) %>%
               group_by(phenotype) %>%
               mutate(padj = p.adjust(Pvalue)),
             res# %>% filter(taxon %in% c("metagTaxa.alpha.st","metagTaxa.pco1"))
  ) %>%
  dplyr::select(-levels,-levels_SampleSize) %>%
  arrange(Pvalue)
write.table(test, file = "Result/Figure6/lmRes_metagTaxa.features.st_withCovariates_padj.txt", quote = F, row.names = F, sep = "\t")


# circos plot 

dat1<-fread("Result/Figure6/lmRes_metagTaxa.features.st_withCovariates_padj.txt")
dat2<-fread("Result/Figure6/lmRes_metagMod.features.st_withCovariates_padj.txt")
dat<-na.omit(rbind(dat1,dat2))

# 创建一个空的列用于存放group信息
dat$group <- ""

# 根据条件填充group列
dat$group[grep("KEGG_", dat$taxon)] <- "KEGG"
dat$group[grep("s__", dat$taxon)] <- "Bacteria"
dat$group[grep("MetaCyc_", dat$taxon)] <- "MetaCyc"
dat$group[dat$group == ""] <- "KOs"

dat<-dat[,c("phenotype","group","taxon","Pvalue")]

dat$log10P<- -log10(dat$Pvalue)
colnames(dat) <- c("Exposure","MicroType","Microbiome","Pval","log10P")

# adjust the log10P value to a range 
library(scales)
dat$flow <-  rescale(dat$log10P, to = c(1.5, 5))

dat$Exposure<-gsub("_log","",dat$Exposure)
dat$Exposure<-gsub("_24","",dat$Exposure)
# dat<-subset(dat,Pval<0.001)
#dat<-subset(dat,Exposure%in%c("Dimethoate","Chlorpyrifos","Atrazine"))
dat<-subset(dat,Exposure%in%c("Dimethoate","Dimethenamid","ers","GDM","OGTT0","OGTT1","OGTT2","HbAlc"))
dat<-subset(dat,Pval<0.01)
dat <- dat %>%
  mutate(Exposure_MicroType = paste(Exposure, MicroType,sep = "_"))

dat.d <- dat %>% 
  mutate(orig = Exposure_MicroType,
         #   flow=log10P,
         dest = Microbiome ) %>%
  dplyr::select(orig, dest, flow, Pval)


# define the colors  

config_type <- cbind.data.frame(
  order = 1:11,
  names = c(unique(dat$Exposure), unique(dat$MicroType)),
  labels = c(unique(dat$Exposure), unique(dat$MicroType)),
  col = c("#FF0000", "#D2960C", "#7DAF00", "#7500FF", "#A0007D", # seven colors for exposures
          "#95C6C6", "#49B1DD", "#7EABCA", "#E98DAF", "#E9CF91",
          "#2F7E77"), #"#9C6625"),  # four colors for Microtype
  stringsAsFactors=F
)

config_element <- cbind.data.frame(
  order = 1:length(c(unique(dat$Exposure_MicroType), unique(dat$Microbiome))),
  names = c(unique(dat$Exposure_MicroType), unique(dat$Microbiome)),
  labels = c(rep(NA, length(unique(dat$Exposure_MicroType))), unique(dat$Microbiome)), 
  stringsAsFactors=F
) %>% 
  mutate(type = sapply(names, function(x){
    if(x %in% dat$Exposure_MicroType) {
      strsplit(x,"_",fixed = T)[[1]][1]
    } else {
      dat$MicroType[which(dat$Microbiome == x)[1]]
    }
  }))%>%
  mutate(typeLabel = sapply(type,
                            function(x) if(x %in% c("ko","arg","vf")) toupper(x) else x )) %>%
  mutate(col = sapply(type, function(x) config_type$col[which(config_type$names == x)])) %>%
  mutate(type_end = sapply(names, function(x){
    if(x %in% dat$Exposure_MicroType) {
      strsplit(x,"_",fixed = T)[[1]][2]
    } else {
      dat$MicroType[which(dat$Microbiome == x)[1]]
    }
  })) %>%
  mutate(col_end = sapply(type_end, function(x) config_type$col[which(config_type$names == x)]))


config_element$type <- factor(config_element$type, levels = c(config_type$labels))

config_element <- config_element %>%
  arrange(names) %>% arrange(type)

# find appropriate label position for bend2
unique(config_element$type)
i <- which(config_element$type == "KEGG")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""


i <- which(config_element$type == "Bacteria")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] <- ""


i <- which(config_element$type == "MetaCyc")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] <- ""


i <- which(config_element$type == "Dimethoate")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] <- ""


i <- which(config_element$type == "Dimethenamid")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""

i <- which(config_element$type == "ers")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""

i <- which(config_element$type == "GDM")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""


i <- which(config_element$type == "HbAlc")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""


i <- which(config_element$type == "OGTT0")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""

i <- which(config_element$type == "OGTT1")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""

i <- which(config_element$type == "OGTT2")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""


# label with KEGG and metacyc names 
path_func<-"Data/Sequencing/metagenomic/Data"
dat_metacyc_name<-read.csv(paste0(path_func,"/metacyc_n.csv"),header = T,sep = ",",row.names = 1)
dat_KEGG_name<-read.csv(paste0(path_func,"/KEGG_n.csv"),header = T,sep = ",",row.names = 1)
# 
config_element$labels_full <-
  sapply(config_element$labels,
         function(x){
           if(is.na(x)){
             NA
           }else if(startsWith(x,"KEGG")){
             dat_KEGG_name$KEGG_name[which(dat_KEGG_name$KEGG_name_new == x)]
           }else if(startsWith(x,"MetaCyc")){
             dat_metacyc_name$pathway[which(dat_metacyc_name$metacyc_name_new == x)]
           }else x
           
         })


# plot with the names of element types 
library(migest)
library(tidyverse)

config_element$type<-factor(config_element$type,levels = c("ers","Dimethoate","Dimethenamid","GDM","OGTT0","OGTT1","OGTT2","HbAlc","Bacteria","KEGG","MetaCyc"))
config_element<-config_element[order(config_element$type),]
pdf("Result/Figure6/circular_typeNames.pdf")
mig_chord(x = dat.d, 
          order = config_element$names,
          grid.col = config_element %>% 
            dplyr::select(names, col_end) %>%
            deframe(), 
          #transparency = dat.d$Pval,
          
          #lab = config_element %>% dplyr::select(names, labels_full) %>% deframe(), 
          lab_bend2 = config_element %>% dplyr::select(names, typeLabel) %>% deframe(), 
          
          gap.degree = 1, 
          label_size = 0.8,
          no_axis = T
          #axis_breaks = 100
)
dev.off()

# plot with the names of elements 
pdf("Result/Figure6/circular_elementFullnames.pdf")
mig_chord(x = dat.d,
          order = config_element$names,
          grid.col = config_element %>%
            dplyr::select(names, col_end) %>%
            deframe(),
          #transparency = dat.d$Pval,
          
          lab = config_element %>% dplyr::select(names, labels_full) %>% deframe(),
          # lab_bend2 = config_element %>% dplyr::select(names, typeLabel) %>% deframe(),
          
          gap.degree = 1,
          label_size = 0.5,
          no_axis = T
          #axis_breaks = 100
)
dev.off()


# Table S2 ----------------------------------------------------------------
library(data.table)
setwd(source_dir_root1)
dir.create("Result/TableS2")
dat1<-fread("Result/mediation/lmRes_bact.features.st_withCovariates_padj.txt")
dat1$Group<-rep("Bacteria 16S",times=nrow(dat1))

dat2<-fread("Result/mediation/lmRes_metagTaxa.features.st_withCovariates_padj.txt")
dat2$Group<-rep("Bacteria metagenome",times=nrow(dat2))


dat3<-fread("Result/mediation/lmRes_metagMod.features.st_withCovariates_padj.txt")

dat3$Group <- ""

# 根据条件填充Group列
dat3$Group[grep("KEGG_", dat3$taxon)] <- "KEGG Module"
dat3$Group[grep("MetaCyc_", dat3$taxon)] <- "MetaCyc Pathway"
dat3$Group[dat3$Group == ""] <- "KOs"

dat<-na.omit(rbind(dat1,dat2))
dat<-rbind(dat,dat3)


path_func<-"Data/Sequencing/metagenomic/Data"
dat_metacyc_name<-read.csv(paste0(path_func,"/metacyc_n.csv"),header = T,sep = ",",row.names = 1)
names(dat_metacyc_name)<-c("ID","Pathway","Description")
dat_KEGG_name<-read.csv(paste0(path_func,"/KEGG_n.csv"),header = T,sep = ",",row.names = 1)
dat_KEGG_name<-dat_KEGG_name[,c(1,3,2)]
names(dat_KEGG_name)<-c("ID","Pathway","Description")
dat_annotation<-rbind(dat_metacyc_name,dat_KEGG_name)


dat_full<-left_join(dat,dat_annotation,by=c("taxon"="ID"))

dat_full2<-subset(dat_full,Pvalue<0.05)
write.csv(dat_full,"Result/TableS2/Result_combine.csv")
write.csv(dat_full2,"Result/TableS2/Result_combine2.csv")


# Figure 6b ---------------------------------------------------------------

library(ggtree)
library(dplyr)
library(miMediation)
# data_total_16S<-read.csv(paste0(path_16s,"data_taxonomy_transformed.csv"),header = T,sep = ",")
# dat_ers<-read.csv("Result/Figure3/data_ers.csv",header = T,sep = ",")
# 
# dat_ers<-merge(dat_ers,data_total_16S[,c("X","number")])
path_meta<-"Data/Sequencing/metagenomic/Data/"

#dir.create("Result/mediation")
data_total_meta<-read.csv(paste0(path_meta,"data_taxonomy_transformed.csv"),header = T,sep = ",")
data_total_meta<-left_join(data_total_meta,dat_ers[,c("number","ers")])
# 
# ers<-dat_ers$ers
#tree<-read.tree("Sequencing/16S/1OTU/decontam/pruned_tree.nwk")

####full taxa

####Species level
Species<-read.csv(paste0(path_meta,"Classes/Species_filtered.csv"),header=T,sep=",")
Species_name<-Species$Species
tree<-Species[,c(1:7)]
#tree<-read.xlsx("metadata/Taxa_16S.xlsx",sheet = 1)

tree$Kingdom<-sub("^k__","",tree$Kingdom)
tree$Phylum<-sub("^p__","",tree$Phylum)
tree$Class<-sub("^c__","",tree$Class)
tree$Order<-sub("^o__","",tree$Order)
tree$Family<-sub("^f__","",tree$Family)
tree$Genus<-sub("^g__","",tree$Genus)
#tree$Species<-sub("^s__","",tree$Species)
# tree<-subset(tree,Family!="Unassigned")
# genus_name_new<-as.character(tree$Genus)
rownames(tree)<-tree$Species
tree<-as.matrix(tree)

rownames(Species)<-Species$Species

M<-Species[,-c(1:7)]

M <- mutate_all(M, function(x) round(x * 100000))
M<-t(M)

Trt<-as.numeric(as.character(data_total_meta$ers))
Y <- as.numeric(as.character(data_total_meta$GDM))

demo.rsltlst <- phyloMed(Trt, M, Y, tree = tree, fdr.alpha =0.2,
                         graph = "circular")

ggsave(filename = "Result/Figure6/Mediation_plot.pdf")

demo.physeq <- demo.rsltlst$clean.data
demo.physeq

result<-demo.rsltlst$rslt$PhyloMed.A
p_node<-result$node.pval
p_clade<-result$sig.clade
p_global<-result$global.pval

write.table(p_node,file = "Result/Figure6/p_node.txt")
write.table(p_global,file = "Result/Figure6/p_global.txt")
write.table(p_clade,file = "Result/Figure6/p_clade.txt")

covariate_dat<-data_total_meta[covariate]
covariate_dat<-as.matrix(covariate_dat)
demo.rsltlst <- phyloMed(Trt, M, Y, tree = tree, confounders = covariate_dat,fdr.alpha =0.2,
                         graph = "circular")

ggsave(filename = "Result/Figure6/Figure6b.pdf")



# TableS6 -----------------------------------------------------------------

dir.create("Result/TableS6")
demo.physeq <- demo.rsltlst$clean.data
demo.physeq

result<-demo.rsltlst$rslt$PhyloMed.A
p_node<-result$node.pval
p_clade<-result$sig.clade
p_global<-result$global.pval

#dir.create(paste0(path_med_metagen,"/phloymed"))
write.table(p_node,file = "Result/TableS6/p_node_adjusted.txt")
write.table(p_global,file = "Result/TableS6/p_global_adjusted.txt")
write.table(p_clade,file = "Result/TableS6/p_clade_adjusted.txt")


save(tree,M,Trt,Y,data_total_meta, file = "Result/Figure6/data_phloymed.RData")


# Figure S8. Phylogenetic tree of branches and nodes in mediating the association of dimethoate and GDM----------------------------------------------------------

dir.create("Result/FigureS8")
library(ggtree)
library(dplyr)
library(miMediation)
# data_total_16S<-read.csv(paste0(path_16s,"data_taxonomy_transformed.csv"),header = T,sep = ",")

####Species level
Species<-read.csv(paste0(path_meta,"Classes/Species_filtered.csv"),header=T,sep=",")
Species_name<-Species$Species
tree<-Species[,c(1:7)]
#tree<-read.xlsx("metadata/Taxa_16S.xlsx",sheet = 1)

tree$Kingdom<-sub("^k__","",tree$Kingdom)
tree$Phylum<-sub("^p__","",tree$Phylum)
tree$Class<-sub("^c__","",tree$Class)
tree$Order<-sub("^o__","",tree$Order)
tree$Family<-sub("^f__","",tree$Family)
tree$Genus<-sub("^g__","",tree$Genus)
#tree$Species<-sub("^s__","",tree$Species)
# tree<-subset(tree,Family!="Unassigned")
# genus_name_new<-as.character(tree$Genus)
rownames(tree)<-tree$Species
tree<-as.matrix(tree)

rownames(Species)<-Species$Species

M<-Species[,-c(1:7)]

M <- mutate_all(M, function(x) round(x * 100000))
M<-t(M)

Trt<-as.numeric(as.character(data_total_meta$Dimethoate_log))
Y <- as.numeric(as.character(data_total_meta$GDM))

demo.rsltlst <- phyloMed(Trt, M, Y, tree = tree, fdr.alpha =0.2,
                         graph = "circular")

ggsave(filename = "Result/FigureS8/Mediation_plot_dimethoate.pdf")

demo.physeq <- demo.rsltlst$clean.data
demo.physeq

result<-demo.rsltlst$rslt$PhyloMed.A
p_node<-result$node.pval
p_clade<-result$sig.clade
p_global<-result$global.pval

#dir.create(paste0(path_med_metagen,"/phloymed"))
write.table(p_node,file = "Result/FigureS8/p_node_dimethoate.txt")
write.table(p_global,file = "Result/FigureS8/p_global_dimethoate.txt")
write.table(p_clade,file = "Result/FigureS8/p_clade_dimethoate.txt")


covariate_dat<-data_total_meta[covariate]
covariate_dat<-as.matrix(covariate_dat)
demo.rsltlst <- phyloMed(Trt, M, Y, tree = tree, confounders = covariate_dat,fdr.alpha =0.2,
                         graph = "circular")

ggsave(filename = "Result/FigureS8/Mediation_plot_adjusted_dimethoate.pdf")

demo.physeq <- demo.rsltlst$clean.data
demo.physeq

result<-demo.rsltlst$rslt$PhyloMed.A
p_node<-result$node.pval
p_clade<-result$sig.clade
p_global<-result$global.pval

dir.create(paste0(path_med_metagen,"/phloymed"))
write.table(p_node,file = "Result/FigureS8/p_node_adjusted_dimethoate.txt")
write.table(p_global,file = "Result/FigureS8/p_global_adjusted_dimethoate.txt")
write.table(p_clade,file = "Result/FigureS8/p_clade_adjusted_dimethoate.txt")


save(tree,M,Trt,Y,data_total_meta, file = "Result/FigureS8/data_phloymed.RData")



# Figure S1. The overview of bacterial taxonomic profiles ---------------------------------------------------------------
setwd("D:/Project/pes_GM_GDM/Data_review")
load("Data/microbialFeatures_unstandardized.RData")

# read meta data 
load("Data/metadata.RData")

load("Data/data_ers.RData")

# 1.bacterial taxa plot-
dir.create("Result/FigureS1")

rownames(amplicon.metadata)<-amplicon.metadata$SampleID
rownames(data_ers)<-data_ers$X
meta <- cbind(amplicon.metadata,bact.features) 
meta<-cbind(meta,data_ers[,c("SampleID","ers")])

dag3_species <- bact.features[,-c(1:2)]

#dag3_species<-as.data.frame(t(dag3_species))

colnames(dag3_species)=lapply(colnames(dag3_species),function(x){
  strsplit(x,"g__")[[1]][2]
})


CompositionTable <- function(x,n) {
  require(foreach)
  #  x[is.na(x)]=0
  mean_value <- data.frame(Taxa=colnames(x), Mean_abundance=colSums(x)/nrow(x))
  most <- as.character(mean_value[order(mean_value$Mean_abundance,decreasing = T),]$Taxa[1:n])
  print(paste("Most abundant taxa is",most,sep = " "))
  
  composition_table <- foreach(i=1:length(most),.combine = rbind) %do%  {
    return.string = data.frame(ID = rownames(x), Relative=x[,most[i]],Level=colnames(x[,most[i],drop=F]))
  }
  
  first <- composition_table[grep(most[1],composition_table$Level),]
  first <- first[order(first$Relative,decreasing = T),]
  level <- as.factor(first$ID)
  composition_table$ID <- factor(composition_table$ID,levels = level)
  
  return(composition_table)
}

dag3_species_plot=CompositionTable(dag3_species,14)
top14_taxa<-unique(dag3_species_plot$Level)

taxa_top<-dag3_species[top14_taxa]
taxa_top <- taxa_top %>%
  mutate(Others = 100 - rowSums(.))

taxa_top_new=taxa_top[order(-taxa_top[,top14_taxa[1]]),]

ID_order1<-rownames(taxa_top_new)


# one box per variety
library(ggsci)
library(reshape2)
library(ggplot2)
taxa_top_new$ID<-rownames(taxa_top_new)
plot_dat<-melt(taxa_top_new,id.vars="ID")
Group<-rep("amplicon",times=nrow(plot_dat))
plot_dat<-data.frame(plot_dat,Group)

####metagen
dag3_species <- metagTaxa.features

#dag3_species<-as.data.frame(t(dag3_species))

colnames(dag3_species)=lapply(colnames(dag3_species),function(x){
  strsplit(x,"s__")[[1]][2]
})

dag3_species_plot=CompositionTable(dag3_species,14)
top14_taxa<-unique(dag3_species_plot$Level)

taxa_top<-dag3_species[top14_taxa]
taxa_top <- taxa_top %>%
  mutate(Others = 100 - rowSums(.))

taxa_top_new=taxa_top[order(-taxa_top[,top14_taxa[1]]),]

ID_order2<-rownames(taxa_top_new)


# one box per variety
library(ggsci)

library(ggplot2)
taxa_top_new$ID<-rownames(taxa_top_new)
plot_dat2<-melt(taxa_top_new,id.vars="ID")
Group<-rep("Metagenomic",times=nrow(plot_dat2))
plot_dat2<-data.frame(plot_dat2,Group)


plot_data_full<-rbind(plot_dat,plot_dat2)
# path_16s_result<-"Result/GM_GDM/16S"
# 
# path_meta_result<-"Result/GM_GDM/metagenomic"


pdf(file="Result/FigureS1/taxplot_detail.pdf",width = 8,height = 6)
ggplot(plot_dat,aes(x=ID,y=value,fill=variable))+geom_bar(stat="identity")+
  labs(x="",y="Relative abundance (%)")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_x_discrete(limits=ID_order1)+
  scale_fill_hue(h=c(10,240),c=30)+
  theme(axis.text.x = element_blank())

dev.off()

pdf(file="Result/FigureS1/taxplot_detail_metagen.pdf",width = 8,height = 6)

ggplot(plot_dat2,aes(x=ID,y=value,fill=variable))+geom_bar(stat="identity")+
  labs(x="",y="Relative abundance (%)")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_x_discrete(limits=ID_order2)+
  scale_fill_hue(h=c(10,240),c=30)+
  theme(axis.text.x = element_blank())

dev.off()






# Figure S5 GDM-associated multi-kingdom microbiome signature---------------------------------------------------------------
#1. LOAD PACKAGES 
library(dplyr)        #v1.0.8   
library(tidyr)        #v1.2.0
library(phyloseq)     #v1.38.0
library(ggplot2)      #v3.3.5
library(ape)          #v5.6.2
library(Maaslin2)     #v1.10.0
library(tableone)     #v0.13.2
library(knitr)        #v1.39
library(RColorBrewer) #v1.1.3
library(forcats)      #v0.5.1
library(scales)       #v1.1.1
library(vegan)        #v2.5.7
library(ggExtra)      #v0.10.0
library(Hmisc)        #v4.7.0
library(corrplot)     #v0.92
library(corpcor)      #1.6.10
library(Boruta)       #v7.0.0
library(caret)        #v6.0.86
library(VIM)          #v6.1.0  
library(stringr)      #v1.4.0
library(egg)          #v0.4.5
library(purrr)        #v0.3.4
library(broom)        #v1.0.1

#2. LOAD DATA
setwd("D:/Project/pes_GM_GDM/Data_review")
load(paste0(path_func,"/metagenomic.RData"))
path_func<-"Data/Sequencing/metagenomic/Data"

dir.create("Result/FigureS5")
#### meta 
path_func_result<-"Result/FigureS5"

## metaphlan3.raw: format for phyloseq

metaphlan3.raw<-read.table(paste0(path_func,"/taxonomy.tsv"),header=T)
metaphlan3.raw<-metaphlan3.raw[,c("clade_name",as.character(meta$SampleID_new))]

humann3.raw<-read.table(paste0(path_func,"/pathabundance_relab_unstratified.tsv"),header = T,sep="\t")
humann3.raw<-humann3.raw[,c("Pathway",as.character(meta$SampleID_new))]
## m3.tree: trim "GCA--" leading string from tip labels---_---------------------
#m3.tree<-read.tree(paste0(path_func,"/mpa_vOct22_CHOCOPhlAnSGB_202212.nwk"))
m3.tree<-read.tree(paste0(path_func,"/mpa_v30_CHOCOPhlAn_201901.nwk"))
genus_colors.df<-read.csv(paste0(path_func,"/genus_color.csv"),header = T,sep = ",")
#save.image(paste0(path_func,"/metagenomic.RData"))

rownames(metaphlan3.raw) <- metaphlan3.raw$clade_name
#metaphlan3.raw[, c('clade_name', 'NCBI_tax_id')] <- list(NULL)
metaphlan3.raw[, c('clade_name')] <- list(NULL)

#take only rows that include s__ (species level annotation)
metaphlan3.species <- metaphlan3.raw[grepl('s__', rownames(metaphlan3.raw)), ]

#exclude rows that include t__
metaphlan3.species<-subset(metaphlan3.species,!grepl("t__",rownames(metaphlan3.species)))

#m3.tree: trim "GCA--" leading string from tip labels
m3.tree$tip.label <- gsub('GCA_[0-9]+/|', '', m3.tree$tip.label)

## Create phyloseq taxa table
species <- data.frame(Names = rownames(metaphlan3.species))
species <- data.frame(do.call('rbind', strsplit(as.character(species$Names),'|',fixed=TRUE)))
rownames(species) <- rownames(metaphlan3.species)
colnames(species) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

species <- as.matrix(species)

## Make phyloseq object-
all(rownames(metaphlan3.species) == rownames(species))
all(colnames(metaphlan3.species) == unlist(meta[,'SampleID_new']))
metaphlan3.species<-as.data.frame(metaphlan3.species)
meta<-as.data.frame(meta)

ps <- phyloseq(otu_table(metaphlan3.species, taxa_are_rows = TRUE),
               sample_data(meta),
               tax_table(species),
               phy_tree(m3.tree))
ps

library(ggtree)
library(ape)
trimed_tree<-phy_tree(ps)
write.tree(trimed_tree,file = paste0(path_func_result,"/trimed_tree.nwk"))

## Filter low abundance taxa
ps2.filt <- filter_taxa(ps, function(x) mean(x) > 0.01, TRUE)
ps2.filt
# This filters out taxa with mean abundance <= 0.1% (Mphln3 relative abundance
# values are out of 100).
# ps2.filt should have 115 taxa across 164 samples, with 64 metadata variables.
detectionRate=10

detect<-function(x){
  length(which(x > 0))*100/length(x)
}

idx3= apply(otu_table(ps2.filt), 1, detect) >detectionRate
remain_name<-taxa_names(ps2.filt)[idx3]

ps2.filt<-prune_taxa(remain_name,ps2.filt)
ps2.filt

#This filters out taxa with detection rate >= 10%

#TAXONOMIC ALPHA DIVERISTY COMPARISONS 
set.seed(123)
ps.int <- transform_sample_counts(ps, function(x) trunc(x*100000))
adiv <- estimate_richness(ps.int, measures=c('Observed', 'Shannon'))

# These alpha div vars are already included in 'meta' dataframe, and therefore in 
# the phyloseq object. Going forward appending of re-derived variables to the 
# phyloseq objects will be commented out:

sample_data(ps2.filt)$Richness <- adiv$Observed
sample_data(ps2.filt)$Shannon <- adiv$Shannon


## Plot alpha diversity by 

# Access metadata and gather to long format by alpha diversity metric
ps2.filt.df <- data.frame(sample_data(ps2.filt))
ps2.filt.df.a <- gather(ps2.filt.df, Alpha_Measure, Value, Richness:Shannon, factor_key=TRUE)


p_alpha<- ggboxplot(ps2.filt.df.a, x="GDM", y="Value",color = "GDM",
                    add = "mean",outlier.shape = NA)+
  #   stat_compare_means(aes(group=group))+
  stat_compare_means(method = "t.test",
                     aes(group=GDM,label = paste0("p = ", ..p.format..)))+
  labs(x="",y="Taxa (MetaPhlAn3)")+
  #scale_color_discrete(labels=c("Low", "Medium", "High"))+
  geom_jitter(aes(color=GDM),width=0.2)+
  scale_color_brewer(palette='Dark2')+
  facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
  #scale_fill_manual(values = rev(cbPalette))+
  scale_x_discrete(limits=c("0","1"),label=c("GDM-", "GDM+"))+
  theme_article()
ggsave(p_alpha,file=paste0(path_func_result,"/FigureS5a.pdf"),width = 8,height = 6)


# 
# p_alpha<- ggboxplot(ps2.filt.df.a, x="Batches", y="Value",color = "Batches",
#                     add = "mean",outlier.shape = NA)+
#   #   stat_compare_means(aes(group=group))+
#   stat_compare_means(method = "t.test",
#                      aes(group=Batches,label = paste0("p = ", ..p.format..)))+
#   labs(x="",y="Taxa (MetaPhlAn3)")+
#   #scale_color_discrete(labels=c("Low", "Medium", "High"))+
#   geom_jitter(aes(color=Batches),width=0.2)+
#   scale_color_brewer(palette='Dark2')+
#   facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
#   #scale_fill_manual(values = rev(cbPalette))+
#   scale_x_discrete(limits=c("1","2"),label=c("Batch 1", "Batch 2"))+
#   theme_article()
# 
# 
# ggsave(p_alpha,file=paste0(path_func_result,"/p_alpha_batches.pdf"),width = 8,height = 6)

## Test for differences in alpha diversity by GDM status
t.obs <- t.test(Richness ~ GDM, data = ps2.filt.df)
t.sha <- t.test(Shannon ~ GDM, data = ps2.filt.df)
t.obs
t.sha


# #5. FIRMICUTES/BACTEROIDES RATIO 
# 
# ## Agglomerate at the phylum level and calculate ratio by individual
# ps2.filt.phylum <- tax_glom(ps2.filt, taxrank = 'Phylum')
# ps2.phylum.df <- psmelt(ps2.filt.phylum)
# 
# ps.phylum <- subset(ps2.phylum.df, 
#                     select=c('Sample', 'Abundance', 'GDM', 'Phylum'))
# 
# # Spread table wide by phylum and calculate ratio
# ps.phylum.wide <- ps.phylum %>% spread(Phylum, Abundance)
# ps.phylum.wide$FBratio <- ps.phylum.wide$p__Firmicutes/ps.phylum.wide$p__Bacteroidetes
# 
# ## Summarize FBratio by group and test for differences
# # Remove 3 samples with Bacteroidetes == 0 or near 0 to avoid 'Inf' or nonsensical
# # ratios. All three samples were from healthy (non preclinical AD) individuals.
# ps.phylum.wide2 <- subset(ps.phylum.wide, !is.infinite(FBratio) & FBratio < 1000) 
# 
# ps.phylum.FBsumm <- ps.phylum.wide2 %>% group_by(GDM) %>%
#   summarise(ci = list(mean_cl_normal(FBratio) %>% rename(mean=y, lwr=ymin, upr=ymax))) %>%
#   unnest(cols = c(ci))
# 
# t.FBratio <- t.test(FBratio ~ GDM, ps.phylum.wide2, 
#                     na.action = na.omit) 
# 
# 
# t.FBratio
# # Retrieve FBratio vector:
# FBratio <- subset(ps.phylum.wide, select=c('Sample', 'FBratio'))
# rownames(FBratio) <- FBratio$Sample
# FBratio <- FBratio %>% rename(Participant = Sample)
# 
# # Add FBratio to the phyloseq object:
# FBratio.ph <- sample_data(FBratio)
# ps2.filt <- merge_phyloseq(ps2.filt, FBratio.ph)
# 



#6. STACKED TAXONOMIC BARPLOTS 

## Agglomerate at the genus level and color/order by phylogeny -----------------
## select 
ps2.filt2 <- filter_taxa(ps2.filt, function(x) mean(x) > 0.2, TRUE)

ps.genus <- tax_glom(ps2.filt2, taxrank = 'Genus')
ps.genus
ps.genus.df <- psmelt(ps.genus)

# genus_selected<-tax_table(ps.genus)
# genus_selected<-as.data.frame(genus_selected)
# write.csv(genus_selected,file = paste0(path_func,"/genus_selected.csv"))
genus_colors.df<-read.csv(paste0(path_func,"/genus_color.csv"),header = T,sep = ",")
# prepare custom color palette for genera
genus_colors <- genus_colors.df$Color
names(genus_colors) <- genus_colors.df$Genus
genusColScale <- scale_fill_manual(values = genus_colors)

# reorder genus levels according to phylogeny (conveniently, this order has  
# been integrated into the custom color palette)
ps.genus.df$Genus <- factor(ps.genus.df$Genus, levels = genus_colors.df$Genus)

## Plot by GDM status (abundances summarized by group)

p_bar_genus.st <- ggplot(ps.genus.df, aes(x=GDM, y=Abundance, fill=Genus))+
  genusColScale+
  geom_bar(stat='identity', position='fill')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))+
  scale_x_discrete(labels=c('GDM-', 'GDM+'))

ggsave(p_bar_genus.st,file=paste0(path_func_result,"/FigureS5c.pdf"),width = 12,height = 8)


#7. TAXONOMIC ORDINATION ANALYSES 
##should seek tree file
## Set up ordination wrapper function

# Required packages: phyloseq, vegan, ggplot2, ggExtra
# Arguments: 
#  phyloseq.obj           = phyloseq object
#  method (string)        = ordination method supported by phyloseq (calls vegan)
#  distance (string)      = between sample distance metric supported by phyloseq 
#                           (see vegdist) 
#  adonisformula (string) = formula to pass to adonis2 (PERMANOVA);
#                           LHS should be 'ps.dist'
#  group (obj)            = variable by to which color sample points
#  seed (int)             = seed for adonis test
#  Rpalette (string)      = name of RColorBrewer palette for plotting
#  markershape (int)      = to recreate, enter 19 for metaphlan, 17 for humann
#  saveplot (logical)     = TRUE if figure should be saved to drive (will overwrite)
#                           NOTE! if TRUE, update filepath within function.
#
# Return:                 
#  list(ordination = ps.ordination, dist = ps.dist, plot = ordination.plot2, 
#       axis1test = axis1.test, axis2test = axis2.test, adonis = ps.adonis)


ordinate2.mphln.AF <- function(phyloseq.obj, method, distance, adonisformula, 
                               group, seed, Rpalette, markershape, saveplot) {
  #ARG HOUSEKEEPING
  metadata <- data.frame(phyloseq::sample_data(phyloseq.obj))
  
  # group.name <- deparse(substitute(group))
  group.name <- group
  group.idx <- grep(group.name, colnames(metadata))
  
  #ORDINATE
  ps.ordination <- phyloseq::ordinate(phyloseq.obj, method=method, distance=distance)
  
  #PLOT
  ordination.plot <- phyloseq::plot_ordination(phyloseq.obj, ps.ordination,
                                               type="samples", color=group.name)+
    geom_point(size=3, shape=markershape)+
    stat_ellipse()+
    scale_color_brewer(palette=Rpalette)+ 
    theme_classic()+
    theme(axis.text.y = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          legend.text = element_text(size=20),
          legend.position = 'left',
          plot.title = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", size=3))
  
  ordination.plot2 <- ggMarginal(ordination.plot, type='boxplot', groupFill=TRUE, 
                                 size=10)
  
  #SAVE PLOT
  phyloseq.obj.name <- deparse(substitute(phyloseq.obj))
  filepath <- paste(paste0(path_func_result,'/'), 
                    Sys.Date(),"_", 
                    phyloseq.obj.name,"_",
                    group.name,"_",
                    method,"_",
                    distance, ".pdf", sep="")
  
  if (saveplot == TRUE){
    ggsave(filepath, 
           plot = ordination.plot2, 
           device ='pdf', 
           width=20, height=13, units='cm')
  }
  
  #SAMPLE COORDINATES: MARGINAL t (or anova)-tests 
  coordinates <- data.frame(ps.ordination$vectors)
  print('Vector lengths: Coord axis 1, Coord axis 2, metadata$group')
  print(c(length(coordinates$Axis.1), length(coordinates$Axis.2), length(metadata[, group.idx])))
  print('Variable compared')
  print(group.name)
  
  df <- data.frame(Coord1 = coordinates$Axis.1,
                   Coord2 = coordinates$Axis.2,
                   Group = metadata[, group.idx])
  
  axis1.test <- aov(Coord1 ~ Group, data=df)
  axis2.test <- aov(Coord2 ~ Group, data=df)
  
  #ADONIS2 (PERMANOVA) TESTs
  ps.dist <- phyloseq::distance(phyloseq.obj, method=distance) 
  
  ps.adonis <- vegan::adonis2(formula(adonisformula),  
                              data=metadata,
                              na=na.omit,
                              permutations = 10000,
                              subset=complete.cases(ps.ordination))  
  
  #RETURN 
  list(ordination = ps.ordination, dist = ps.dist, plot = ordination.plot2, 
       axis1test = axis1.test, axis2test = axis2.test, 
       adonis = ps.adonis)
  
}



## Execute PCoA on taxonomic abundances, with group comparisons

adonis.formula <- 'ps.dist ~ Age + pre_BMI + Educational_level + pre_BMI + weightgain + Parity + passsmokHis_1y+Batches'

ordination3b <- ordinate2.mphln.AF(ps2.filt, method = 'PCoA', distance = 'unifrac', 
                                   adonisformula = adonis.formula,
                                   group = "GDM", seed = 2023,
                                   Rpalette = 'Dark2', markershape = 19, 
                                   saveplot = TRUE)

adonis.formula2 <- 'ps.dist ~ Age + pre_BMI + Educational_level + pre_BMI + weightgain + Parity + passsmokHis_1y'

ordination3b2 <- ordinate2.mphln.AF(ps2.filt, method = 'PCoA', distance = 'unifrac', 
                                    adonisformula = adonis.formula,
                                    group = "Batches", seed = 2023,
                                    Rpalette = 'Dark2', markershape = 19, 
                                    saveplot = TRUE)

# Add the PCOA1 and PCOA2 coordinates to the phyloseq obj metadata
sample_data(ps2.filt)$Tax.PCoA1 <- ordination3b$ordination$vectors[,1]
sample_data(ps2.filt)$Tax.PCoA2 <- ordination3b$ordination$vectors[,2]

## CAP analysis of taxonomic abundances

# CAP (see Anderson and Willis, Ecology 2003) has a number of implementations in
# R. The phyloseq implementation calls vegan::capscale, and is equivalent to 
# carrying out vegan::rda on PCoA coordinate vectors. The same 'distance' matrix
# that was calculated for the original PCoA is used, and in this case we fit the
# same model as in the PERMANOVA (adonis2) carried out in the above section.
# One idiosyncrasy is that rather than use the original PCoA coordinate vectors
# (generated in phyloseq by calling ape::pcoa), the phyloseq implementation of 
# CAP re-performs PCoA with vegan::cmdscale, and uses these as input. This has
# negligible impact on results, however.


# Access distance matrix from PCoA ordination and run CAP 
u.dist <- ordination3b$dist
ps2.filt.cap <- ordinate(ps2.filt,
                         method = 'CAP',
                         distance = u.dist,
                         formula = ~Age + pre_BMI + Educational_level + pre_BMI + weightgain + Parity + passsmokHis_1y+Batches,
                         na.action = na.omit)


# Biplot (out of curiosity)
cap.biplot <- ps2.filt.cap$CCA$biplot
cap.biplot
# Plot. Note! phyloseq::plot_ordination tries to add values for % variance
# explained on each of the axes-- but these values are not relevant to CAP, 
# and indeed CAP plots are not shown with any '% variance explained' information
# on the axes. These auto-populated values should be ignored. 
p_ps2.filt.cap <- plot_ordination(ps2.filt, ps2.filt.cap, color = 'GDM')+
  geom_point(size=3)+
  stat_ellipse()+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = 'left',
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))
p_ps2.filt.cap
p_ps2.filt.cap2 <- ggMarginal(p_ps2.filt.cap, type='boxplot', groupFill=TRUE, size=10)
p_ps2.filt.cap2
ggsave(p_ps2.filt.cap2,file=paste0(path_func_result,"/FigureS5b.pdf"),height = 6,width = 8)
# Test for the significance of the constraints (model terms)
cap.anova <- anova(ps2.filt.cap, by='term', permutations=10000)


# Access CAP axis coordinates and test for significant differences on CAP1 and CAP2 
ps2.filt.cap.df <- p_ps2.filt.cap$data

cap1.t <- t.test(CAP1 ~ GDM, data=ps2.filt.cap.df) 
cap2.t <- t.test(CAP2 ~ GDM, data=ps2.filt.cap.df)

cap1.t ## p=0.008696
cap2.t ## p=0.2742



###batches
u.dist2 <- ordination3b$dist
ps2.filt.cap2 <- ordinate(ps2.filt,
                          method = 'CAP',
                          distance = u.dist,
                          formula = ~Age + pre_BMI + Educational_level + pre_BMI + weightgain + Parity + passsmokHis_1y,
                          na.action = na.omit)


# Biplot (out of curiosity)
cap.biplot2 <- ps2.filt.cap2$CCA$biplot
cap.biplot2
# Plot. Note! phyloseq::plot_ordination tries to add values for % variance
# explained on each of the axes-- but these values are not relevant to CAP, 
# and indeed CAP plots are not shown with any '% variance explained' information
# on the axes. These auto-populated values should be ignored. 
dir.create("Result/FigureS10")
###Figure S10b
p_ps2.filt.cap2 <- plot_ordination(ps2.filt, ps2.filt.cap2, color = 'Batches')+
  geom_point(size=3)+
  stat_ellipse()+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = 'left',
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))
p_ps2.filt.cap
p_ps2.filt.cap3 <- ggMarginal(p_ps2.filt.cap2, type='boxplot', groupFill=TRUE, size=10)
p_ps2.filt.cap3
ggsave(p_ps2.filt.cap2,file="Result/FigureS10/FigureS10b.pdf",height = 6,width = 8)

# Access CAP axis coordinates and test for significant differences on CAP1 and CAP2 
ps2.filt.cap.df <- p_ps2.filt.cap2$data

cap1.t <- t.test(CAP1 ~ GDM, data=ps2.filt.cap.df) 
cap2.t <- t.test(CAP2 ~ GDM, data=ps2.filt.cap.df)

cap1.t ## p=0.1309
cap2.t ## p=0.5082
#8. CREATE PATHWAYS PHYLOSEQ OBJECT 

## HOUSEKEEPING ##

## meta: source metadata from current phyloseq object (ps2.filt)
meta.psh <- data.frame(sample_data(ps2.filt))

## humann3.raw: format for phyloseq
rownames(humann3.raw) <- humann3.raw$Pathway
humann3.raw$Pathway <- NULL

# Omit rows that break down the species associations (contain '|' delimiter)
# humann3 <- humann3.raw[!grepl('/|', rownames(humann3.raw)), ]
# humann3<-na.omit(humann3)
humann3<-na.omit(humann3.raw)
# Check pathway relative abundances sum to 1 (unlike taxa, which sum to 100)
humann3.sums <- colSums(humann3)
humann3.sums
samples_humann <- colnames(humann3) 
#N = 456, creation of the phyloseq object will subset to those in meta.psh

## Create phyloseq 'taxa' table, adapted for pathways
pathways <- data.frame(Pathway=rownames(humann3))
rownames(pathways) <- pathways$Pathway
pathways <- as.matrix(pathways)

## Create pathways phyloseq object
psh <- phyloseq(otu_table(humann3, taxa_are_rows = TRUE),
                sample_data(meta.psh),
                tax_table(pathways))
psh

# 457 samples by 381 pathways
## Prune unintegrated/unmapped pathways and renormalize
psh <- subset_taxa(psh, Pathway!="UNMAPPED"& Pathway!="UNINTEGRATED") #379 pathways
psh  <- transform_sample_counts(psh, function(x) (x / sum(x))*100 )

## Filter low abundance pathways
psh.filt <- filter_taxa(psh, function(x) mean(x) > 0.01, TRUE) 
psh.filt
# This filters out pathways with mean abundance <= 0.01% (Humann3 relative abundance
# values were re-normalized to 100 in the previous section).
# psh.filt should have 233 pathways across 457 samples, with 212 metadata variables.

#9. PATHWAY ALPHA DIVERSITY COMPARISONS 

# Alpha diversity metrics should be calculated using the unfiltered abundance 
# data, and require count data rather than fractional relative abundance: so,
# transform unfiltered relative abundances to integer counts, preserving 
# precision to 0.00001%. Function estimate_richness() will throw warning about 
# singletons that may be ignored (more relevant to ASV data).

psh.alpha.int <- transform_sample_counts(psh, function(x) trunc(x*100000))
h.adiv <- estimate_richness(psh.alpha.int, measures=c('Observed', 'Shannon'))

# These alpha div vars were already included in 'meta.psh' dataframe, and therefore in 
# the phyloseq object. Appending of re-derived variables to the phyloseq objects 
# is therefore commented out:
# 
# rownames(h.adiv) <- gsub("X", "", rownames(h.adiv))
# 
# # Add pathway alpha diversities to humann phyloseq object (no need to reorder)
sample_data(psh.filt)$Fnl.Richness <- h.adiv$Observed
sample_data(psh.filt)$Fnl.Shannon <- h.adiv$Shannon
# 
# # Add pathway alpha diversities to mphln phyloseq object 
# # Important! Sample order is inverted. This is handled by converting dataframe
# # to phyloseq 'sample_data' object before addition of new metadata. Actually,
# # this is the safe way to add any additional metadata to a phyloseq object.
# 
# h.adiv.ph <- sample_data(h.adiv)
# sample_data(ps2.filt)$Fnl.Richness <- h.adiv.ph$Fnl.Richness
# sample_data(ps2.filt)$Fnl.Shannon <- h.adiv.ph$Fnl.Shannon

## Plot pathway alpha diversity by GDM status

# Access metadata and gather to long format by alpha diversity metric
psh.filt.df <- data.frame(sample_data(psh.filt))
psh.filt.df.a <- gather(psh.filt.df, Alpha_Measure, Value, Fnl.Richness:Fnl.Shannon, factor_key=TRUE)



p_alpha_h<- ggboxplot(psh.filt.df.a, x="GDM", y="Value",color = "GDM",
                      add = "mean",outlier.shape = NA)+
  #   stat_compare_means(aes(group=group))+
  stat_compare_means(method = "t.test",
                     aes(group=GDM,label = paste0("p = ", ..p.format..)))+
  labs(x="",y="Pathways (HUMAnN3))")+
  #scale_color_discrete(labels=c("Low", "Medium", "High"))+
  geom_jitter(aes(color=GDM),width=0.2)+
  scale_color_brewer(palette='Dark2')+
  facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
  #scale_fill_manual(values = rev(cbPalette))+
  scale_x_discrete(limits=c("0","1"),label=c("GDM-", "GDM+"))+
  theme_article()
ggsave(p_alpha_h,file=paste0(path_func_result,"/p_alpha_pathway.pdf"),width = 8,height = 6)


# FigureS10 Batch effect --------------------------------------------------


dir.create("Result/FigureS10")

psh.filt.df.a.plot<-subset(psh.filt.df.a,Alpha_Measure=="Fnl.Shannon")
p_alpha<- ggboxplot(psh.filt.df.a.plot, x="Batches", y="Value",color = "Batches",
                    add = "mean",outlier.shape = NA)+
  #   stat_compare_means(aes(group=group))+
  stat_compare_means(method = "t.test",
                     aes(group=Batches,label = paste0("p = ", ..p.format..)))+
  labs(x="",y="Taxa (MetaPhlAn3)")+
  #scale_color_discrete(labels=c("Low", "Medium", "High"))+
  geom_jitter(aes(color=Batches),width=0.2)+
  scale_color_brewer(palette='Dark2')+
  #facet_wrap(~Alpha_Measure, nrow=1, scales = 'free_y')+
  #scale_fill_manual(values = rev(cbPalette))+
  scale_x_discrete(limits=c("1","2"),label=c("Batch 1", "Batch 2"))+
  theme_article()
# theme(axis.text.x=element_text(size=8,vjust = 1,hjust = 1))+
# theme(panel.background = element_rect(color = "black"),
#       legend.position = "none")#+

ggsave(p_alpha,file="Result/FigureS10/FigureS10a.pdf",width = 5,height = 6)




# Figure S5b --------------------------------------------------------------


#10. PATHWAY ORDINATION ANALYSES 

## Set up ordination wrapper function-------------------------------------------
# Note: this function is identical to ordinate2.mphln.AF with the addition of the
# 'binary' argument to enable use of binary Bray-Curtis index.
# An intrepid scholar may one day combine these into one streamlined function. 
#
# Required packages: phyloseq, vegan, ggplot2, ggExtra
# Arguments: 
#  phyloseq.obj           = phyloseq object
#  method (string)        = ordination method supported by phyloseq (calls vegan)
#  distance (string)      = between sample distance metric supported by phyloseq 
#                           (see vegdist) 
#  adonisformula (string) = formula to pass to adonis2 (PERMANOVA);
#                           LHS should be 'ps.dist'
#  group (obj)            = variable by to which color sample points
#  seed (int)             = seed for adonis test
#  Rpalette (string)      = name of RColorBrewer palette for plotting
#  markershape (int)      = to recreate, enter 19 for metaphlan, 17 for humann
#  binary (logical)       = TRUE if the vegdist distance metric should be the 
#                           binary variant.
#  saveplot (logical)     = TRUE if figure should be saved to drive (will overwrite)
#                           NOTE! if TRUE, update filepath within function.
#
# Return:                 
#  list(ordination = ps.ordination, dist = ps.dist, plot = ordination.plot2, 
#       axis1test = axis1.test, axis2test = axis2.test, adonis = ps.adonis)


ordinate2.humann.AF <- function(phyloseq.obj, method, distance, adonisformula, 
                                group, seed, Rpalette, markershape, binary, 
                                saveplot) { 
  #ARG HOUSEKEEPING
  metadata <- data.frame(phyloseq::sample_data(phyloseq.obj))
  
  # group.name <- deparse(substitute(group))
  group.name <- group
  group.idx <- grep(group.name, colnames(metadata))
  
  
  #ORDINATE
  ps.ordination <- phyloseq::ordinate(phyloseq.obj, method=method, 
                                      distance=distance, binary = binary)
  
  #PLOT
  ordination.plot <- phyloseq::plot_ordination(phyloseq.obj, ps.ordination, 
                                               type="samples", color=group.name)+
    geom_point(size=3, shape=markershape)+
    stat_ellipse()+
    scale_color_brewer(palette=Rpalette)+ 
    theme_classic()+
    theme(axis.text.y = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          legend.text = element_text(size=20),
          legend.position = 'left',
          plot.title = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", size=3))
  
  ordination.plot2 <- ggMarginal(ordination.plot, type='boxplot', groupFill=TRUE, 
                                 size=10)
  
  #SAVE PLOT
  phyloseq.obj.name <- deparse(substitute(phyloseq.obj))
  filepath <- paste(paste0(path_func_result,'/'), 
                    Sys.Date(),"_", 
                    phyloseq.obj.name,"_",
                    group.name,"_",
                    method,"_",
                    distance, ".pdf", sep="")
  
  if (saveplot == TRUE){
    ggsave(filepath, 
           plot = ordination.plot2, 
           device ='pdf', 
           width=20, height=13, units='cm')
  }
  
  #SAMPLE COORDINATES: MARGINAL t (or anova)-tests (print to check)
  coordinates <- data.frame(ps.ordination$vectors)
  print('Vector lengths: Coord axis 1, Coord axis 2, metadata$group')
  print(c(length(coordinates$Axis.1), length(coordinates$Axis.2), length(metadata[, group.idx])))
  print('Variable compared')
  print(group.name)
  
  df <- data.frame(Coord1 = coordinates$Axis.1,
                   Coord2 = coordinates$Axis.2,
                   Group = metadata[, group.idx])
  
  axis1.test <- aov(Coord1 ~ Group, data=df)
  axis2.test <- aov(Coord2 ~ Group, data=df)
  
  #ADONIS TESTs 
  ps.dist <- phyloseq::distance(phyloseq.obj, method=distance, binary=binary)
  
  ps.adonis <- vegan::adonis2(formula(adonisformula),  
                              data=metadata,
                              na=na.omit,
                              permutations = 10000,
                              subset=complete.cases(ps.ordination))  
  
  #RETURN 
  list(ordination = ps.ordination, dist = ps.dist, plot = ordination.plot2, 
       axis1test = axis1.test, axis2test = axis2.test, 
       adonis = ps.adonis)
  
}



## Execute PCoA on pathway abundances, with group comparisons-------------------

# Using same model formula as before:
#adonis.formula <- 'ps.dist ~ Age + APOE4.status + Diabetes + BMI + Hypertension + Interval_Days_Amyloid + GDM'

adonis.formula <- 'ps.dist ~ Age + pre_BMI + Educational_level + pre_BMI + weightgain + Parity + passsmokHis_1y+Batches'

h.ordination4b <- ordinate2.humann.AF(psh.filt, method = 'PCoA', distance = 'bray', binary = TRUE, 
                                      adonisformula = adonis.formula,
                                      group = "GDM", seed = 13,
                                      Rpalette = 'Dark2', markershape=17, saveplot = TRUE)
# 
# ordination4b <- ordinate2.mphln.AF(psh.filt, method = 'PCoA', distance = 'unifrac', 
#                                    adonisformula = adonis.formula,
#                                    group = "GDM", seed = 2023,
#                                    Rpalette = 'Dark2', markershape = 17,
#                                    saveplot = TRUE)
# # Add the pathway PCOA1 and PCOA2 coordinates to the phyloseq objects
hmn.pcoa.vectors <- data.frame(h.ordination4b$ordination$vectors[,1:2])
colnames(hmn.pcoa.vectors) <- c('Fnl.PCoA1', 'Fnl.PCoA2')
hmn.pcoa.vectors.ph <- sample_data(hmn.pcoa.vectors)

ps2.filt <- merge_phyloseq(ps2.filt, hmn.pcoa.vectors.ph)
psh.filt <- merge_phyloseq(psh.filt, hmn.pcoa.vectors.ph)



## CAP analysis of pathway abundances-------------------------------------------

# Access distance matrix from PCoA ordination and run CAP 
bb.dist <- h.ordination4b$dist
psh.filt.cap <- ordinate(psh.filt,
                         method = 'CAP',
                         distance = bb.dist,
                         formula = ~Age + pre_BMI + Educational_level + pre_BMI + weightgain + Parity + passsmokHis_1y+Batches+GDM,
                         na.action = na.omit)

# Biplot 
h.cap.biplot <- psh.filt.cap$CCA$biplot

h.cap.biplot
# Plot. Note! phyloseq::plot_ordination tries to add values for % variance
# explained on each of the axes-- but these values are not relevant to CAP, 
# and indeed CAP plots are not shown with any '% variance explained' information
# on the axes. These auto-populated values should be ignored. 
p_psh.filt.cap <- plot_ordination(psh.filt, psh.filt.cap, color = 'GDM')+
  geom_point(size=3, shape=17)+
  stat_ellipse()+
  scale_color_brewer(palette='Dark2')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = 'left',
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))

p_psh.filt.cap2 <- ggMarginal(p_psh.filt.cap, type='boxplot', groupFill=TRUE, size=10)
p_psh.filt.cap2
ggsave(p_psh.filt.cap2,file=paste0(path_func_result,"/p_PCOA_pathway.pdf"),width = 8,height = 6)

# Test for the significance of the constraints (model terms) 
h.cap.anova <- anova(psh.filt.cap, by='term', permutations=10000)

# Access CAP axis coordinates and test for significant differences on CAP1 and CAP2 
psh.filt.cap.df <- p_psh.filt.cap$data

h.cap1.t <- t.test(CAP1 ~ GDM, data=psh.filt.cap.df)
h.cap2.t <- t.test(CAP2 ~ GDM, data=psh.filt.cap.df)
h.cap1.t #0.6939
h.cap2.t #0.0355

#11. PAIRWISE BIOMARKER CORRELATION ANALYSES 

## needing tree information
## Prepare biomarker data

# Access sample metadata, including GM metrics of alpha diversity and PCoA 
# coordinates 
metadata.corr <- data.frame(sample_data(ps2.filt))

# Subset to biomarkers and GM metrics of interest
corr.vars <- c('Tax.PCoA1', 'Tax.PCoA2', 'Fnl.PCoA1', 'Fnl.PCoA2',
               'OGTT0_24', 'OGTT1_24', 
               'OGTT2_24', 'HbAlc_24',
               'pre_BMI', 'Age') 

metadata.corr2 <- metadata.corr[, corr.vars]

# # Numerically encode APOE e4 carrier status.
# metadata.corr2 <- metadata.corr2 %>% 
#   mutate(APOE4.status = if_else(APOE4.status == 'e4-', 0, 1))


## Determine pairwise spearman correlations and plot correlogram
corr.spearman <- rcorr(as.matrix(metadata.corr2, type='spearman'))

# Convert symmetric p-value matrix to non-redundant vector, adjust p-values, and 
# return to named symmetric matrix. Convert diagonal NAs to 0 for visualization. 
corr.spearman.p <- sm2vec(corr.spearman$P, diag = FALSE)
idx <- which(corr.spearman.p < 0.05)
corr.spearman.p[idx] <- p.adjust(corr.spearman.p[idx], method='BH')
corr.spearman.p.adj <- vec2sm(corr.spearman.p, diag = FALSE)
colnames(corr.spearman.p.adj) <- colnames(corr.spearman$P)
rownames(corr.spearman.p.adj) <- rownames(corr.spearman$P)
corr.spearman.p.adj[is.na(corr.spearman.p.adj)] <- 0


##FigureS5d 手动保存
p_corr.spearman <- corrplot(corr.spearman$r, order = 'original', p.mat=corr.spearman.p.adj,
                            sig.level=0.05, addCoef.col='black', insig='blank')
ggsave(p_corr.spearman,file="Result/FigureS5/FigureS5d.pdf",width = 6,height = 6)

#ggsave(p_corr.spearman,file=paste0(path_func_result,"p_corr.spearman.pdf"),width = 6,height = 6)
#should mannually save the figures



# Figure S6 Association of exposure factors with microbiome features---------------------------------------------------------------
dir.create("Result/FigureS6")
# pesticide and gut microbiota--manhatoon plot 

dat_plot<-read.csv("Result/FigureS6/metagenomic/reg/linear.csv",header = T,sep = ",")
dat_plot2<-read.csv("Result/FigureS6/16S/reg/linear.csv",header = T,sep = ",")

library(qqman)
library(tidyverse)
library(RcmdrMisc)
library(correlation)

#manhattan(dat_cor, chr="X1",  snp="Query", p="p" )
library(hudson)
dat_plot<-data.frame(dat_plot,str_split_fixed(dat_plot$species, "_log", 2))
shape<-rep("log",times=nrow(dat_plot))
dat_plot<-data.frame(dat_plot,shape)

ewas.d1<-dat_plot[,c(2,10,11,13)]

dat_plot2<-data.frame(dat_plot2,str_split_fixed(dat_plot2$genus, "_log", 2))
shape<-rep("log",times=nrow(dat_plot2))
dat_plot2<-data.frame(dat_plot2,shape)

ewas.d1<-dat_plot[,c(2,10,11,13)]
ewas.d2<-dat_plot2[,c(2,10,11,13)]
names(ewas.d1)<-c("Variable","pvalue","Group","Shape")
names(ewas.d2)<-c("Variable","pvalue","Group","Shape")
# emirror(top=ewas.d1, bottom=ewas.d1, annotate_p = c(0.05,0.05), 
#         highlight_p=c(0.05,0.05), highlighter="green", 
#         toptitle = "Metals and Metabolites Association", 
#         bottomtitle = "EWAS Comparison",
#         color1 = "#00A087",
#         color2="#E64B35",
#         type = "pdf",
#         file = "Result/pes_GM/metagenomic/reg/metal_matabolites_new",hgt = 16,wi=12)

ewas.d1<-subset(ewas.d1,Variable%in%Species_name)


#ewas.d2<-subset(ewas.d2,Variable%in%Genus_name)
emirror(top=ewas.d1, bottom=ewas.d2, annotate_p = c(0.001,0.001), 
        highlight_p= c(0.001,0.001), highlighter="#E64B35", 
        toptitle = "Pesticides and microbial species Association", 
        bottomtitle = "Pesticides and microbial genus Association",
        annotate_var="Null",
        color1 = "#8491B4",
        color2="#91D1C2",
        rotatelabels=TRUE,
        labelangle = 45,
        type = "pdf",
        hgtratio=0.8,
        file = "Result/FigureS6/FigureS6",hgt = 16,wi=12)

# emirror(top=ewas.d1, bottom=ewas.d2, annotate_p = c(0.0001,0.0001), 
#         highlight_p= c(0.0001,0.0001), highlighter="#E64B35", 
#         toptitle = "Pesticides and microbial species Association", 
#         bottomtitle = "Pesticides and microbial genus Association",
#         annotate_var="Null",
#         color1 = "#8491B4",
#         color2="#91D1C2",
#         rotatelabels=TRUE,
#         labelangle = 45,
#         type = "pdf",
#         hgtratio=0.8,
#         file = "Result/pes_GM/metagenomic/reg/Manhatton_adjusted_species_0.0001",hgt = 16,wi=12)


# Figure S7 Microbial interaction networks of bacterial genera among individuals with high and low ers index---------------------------------------------------------------

setwd(source_dir_root1)
load("Data/microbialFeatures_standardized.RData")
dir.create("Result/FigureS7")
# read meta data 
load("Data/metadata.RData")

load("Data/data_ers.RData")

# # amplicon ----------------------------------------------------------------
# 
# 
# High_id<-subset(data_ers,ers>=median(data_ers$ers))$SampleID
# Low_id<-subset(data_ers,ers<median(data_ers$ers))$SampleID
# 
# abund_High<-bact.features.st[High_id,-c(1:2)]
# abund_High<-as.data.frame(t(abund_High))
# #abund_High$x<-rownames(abund_High)
# 
# 
# abund_Low<-bact.features.st[Low_id,-c(1:2)]
# abund_Low<-as.data.frame(t(abund_Low))
# #abund_Low$x<-rownames(abund_Low)
# 
# library(NetMoss2)
# library(Hmisc)
# library(dplyr)
# 
# cr<-t(abund_Low)
# net_Low<-rcorr(cr,type="spearman")
# net_Low<-net_Low$r
# 
# 
# cr<-t(abund_High)
# net_High<-rcorr(cr,type="spearman")
# net_High<-net_High$r
# 
# nodes_result_low.vs.High = 
#   NetMoss(case_dir = abund_High,
#           control_dir = abund_Low,
#           net_case_dir = net_High,
#           net_control_dir = net_Low)
# 
# result = nodes_result_healthy.vs.High[[1]]
# 
# dir.create("Result/GM_High/16S/NetMoss")
# write.table(result, file = "Result/GM_High/16S/NetMoss/NetMoss_ers.txt", sep = "\t", quote = F, row.names = F)
# 
# detach("egg")
# library(ggpubr)
# #run netplot2
# setwd("D:/Project/pes_GM_GDM/Result/GM_GDM/16S/NetMoss")
# #plot networks
# netPlot2(result = nodes_result_healthy.vs.High,
#          num.top = 5,
#          num.score = 30,
#          e.th = 0.4,
#          my.layout = layout_as_star,
#          my.label = TRUE)
# 
# ggsave(filename = "Network_plot_ers.pdf",height = 8,width = 12)
# 
# #plot roc 
# #trim markers
# marker = data.frame(result[which(result$p.adj < 0.05),])
# marker = data.frame(marker[which(marker$NetMoss_Score > 0.3),])   ####marker selection
# rownames(marker) = marker$taxon_names
# 
# #construct metadata    ######if file exists, skip
# #metadata
# case = nodes_result_healthy.vs.High[[4]]
# Low = nodes_result_healthy.vs.High[[5]]
# metadata = data.frame(sample_id = c(colnames(case[,-1]),
#                                     colnames(Low[,-1])),
#                       type = c(rep("High",length(colnames(case[,-1]))),
#                                rep("Low",length(colnames(Low[,-1])))))
# 
# #metadata <- metadata[!metadata$sample_id %in% c("x"), ]
# 
# metadata$sample_id = as.character(metadata$sample_id)
# metadata$type = as.factor(metadata$type)
# rownames(metadata) = metadata$sample_id
# myROC = netROC(case_dir = abund_High,
#                control_dir = abund_Low,
#                marker = marker,
#                metadata = metadata,
#                plot.roc = TRUE, 
#                train.num = 20)    ####image saved


# metagenomic -------------------------------------------------------------

High_id<-subset(data_ers,ers>=median(data_ers$ers))$number
Low_id<-subset(data_ers,ers<median(data_ers$ers))$number

High_id_meta<-subset(metadata.meta,number%in%High_id)$SampleID
Low_id_meta<-subset(metadata.meta,number%in%Low_id)$SampleID
# 
# High_id<-subset(metadata.meta,High==1)$SampleID
# Low_id<-subset(metadata.meta,High==0)$SampleID

abund_High<-metagTaxa.features.st[High_id_meta,]
abund_High<-as.data.frame(t(abund_High))
#abund_High$x<-rownames(abund_High)


abund_Low<-metagTaxa.features.st[Low_id_meta,]
abund_Low<-as.data.frame(t(abund_Low))
#abund_Low$x<-rownames(abund_Low)

library(NetMoss2)
library(Hmisc)
library(dplyr)

cr<-t(abund_Low)
net_Low<-rcorr(cr,type="spearman")
net_Low<-net_Low$r


cr<-t(abund_High)
net_High<-rcorr(cr,type="spearman")
net_High<-net_High$r

# NetMoss 
nodes_result_healthy.vs.High = 
  NetMoss(case_dir = abund_High,
          control_dir = abund_Low,
          net_case_dir = net_High,
          net_control_dir = net_Low)



NMSS_1.2 = nodes_result_healthy.vs.High[[1]]
dir.create("Result/FigureS7")
write.table(NMSS_1.2, file = "Result/FigureS7/NetMoss_ers.txt", sep = "\t", quote = F, row.names = F)

setwd("D:/Project/pes_GM_GDM/Data_review/Result/FigureS7")
#plot networks
library(NetMoss2)
netPlot(result = nodes_result_healthy.vs.High,
        num.top = 5,
        num.score = 30,
        e.th = 0.5,
        my.layout = layout_as_star,
        my.label = TRUE)

ggsave(filename = "Network_plot_ers.pdf",height = 8,width = 12)

#plot roc 
#trim markers
result<-NMSS_1.2
marker = data.frame(result[which(result$p.adj < 0.05),])
marker = data.frame(marker[which(marker$NetMoss_Score > 0.3),])   ####marker selection
rownames(marker) = marker$taxon_names

#construct metadata    ######if file exists, skip
#metadata
case = nodes_result_healthy.vs.High[[4]]
Low = nodes_result_healthy.vs.High[[5]]
metadata = data.frame(sample_id = c(colnames(case[,-1]),
                                    colnames(Low[,-1])),
                      type = c(rep("High",length(colnames(case[,-1]))),
                               rep("Low",length(colnames(Low[,-1])))))

# calculate and plot the networks in three datasets
setwd("D:/Project/pes_GM_GDM/Data_review")
load("Data/microbialFeatures_standardized.RData")

# read meta data 
load("Data/metadata.RData")

load("Data/data_ers.RData")


#2. metagenomic 

rownames(metadata.meta)<-metadata.meta$SampleID
meta<-cbind(metadata.meta,metagTaxa.features.st)


High_id<-subset(data_ers,ers>=median(data_ers$ers))$number
Low_id<-subset(data_ers,ers<median(data_ers$ers))$number

High_id_meta<-subset(metadata.meta,number%in%High_id)$SampleID
Low_id_meta<-subset(metadata.meta,number%in%Low_id)$SampleID

meta$Disease<-ifelse(meta$SampleID%in%High_id_meta,"High","Low")

combined.data<-meta[,c(colnames(metagTaxa.features.st),"Disease")]
Row.names<-rownames(combined.data)

combined.data<-data.frame(Row.names,combined.data)
rownames(combined.data)<-NULL
# calculate networks ===================
library(dplyr)
library(data.table)
library(Hmisc)


edges_3Diseases <- NULL
for(dss in unique(meta$Disease)){
  
  dat.tmp <- combined.data %>% 
    filter(Disease == dss) %>% 
    tibble::column_to_rownames("Row.names") %>%
    dplyr::select(-Disease)
  
  Corr <- rcorr(as.matrix(dat.tmp) , type="spearman")
  occor.r <- Corr$r
  occor.p <- Corr$P
  
  # Hide lower triangle for r (so each correlation is only going to appear once)
  lower <- occor.r
  lower[lower.tri(occor.r)]<-NA
  
  # edges
  r_df.l <- lower %>% reshape2::melt(value.name = "r")
  p_df.l <- occor.p %>% reshape2::melt(value.name = "p")
  
  if(all(r_df.l$Var1 == p_df.l$Var1) & all(r_df.l$Var2 == p_df.l$Var2)){
    occor.rp_l <- cbind.data.frame(r_df.l, p=p_df.l$p, stringsAsFactors=F) # quicker
  }else{
    occor.rp_l <- merge(occor.r %>% reshape2::melt(value.name = "r"), 
                        occor.p %>% reshape2::melt(value.name = "p"),
                        by=c("Var1", "Var2"))  # too slowwwwwwwww
  }
  
  
  edges <- occor.rp_l %>% 
    filter(!is.na(r)) %>% 
    filter(abs(r) < 1) %>% 
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) 
  
  edges$edge = paste(edges$Var1, edges$Var2, sep = "|")
  edges$disease = dss
  edges$edgeType <-
    paste0(substr(edges$edge,1,1), 
           substr(sapply(strsplit(edges$edge, "|", fixed = T), "[[", 2), 1, 1))
  
  edges_3Diseases <- bind_rows(edges_3Diseases, edges)
  
}

edges_3Diseases$r <- as.numeric(edges_3Diseases$r)
edges_3Diseases$p <- as.numeric(edges_3Diseases$p)

# edges_3Diseases$source <- sapply(edges_3Diseases$Var1,
#                                  function(x) microbeAbb$V2[which(microbeAbb$abb == x)])
# edges_3Diseases$target <- sapply(edges_3Diseases$Var2,
#                                  function(x) microbeAbb$V2[which(microbeAbb$abb == x)])

edges_3Diseases$source<-edges_3Diseases$Var1
edges_3Diseases$target<-edges_3Diseases$Var2


for(dss in c( "High","Low")){
  
  edges <- edges_3Diseases %>% 
    filter(disease == dss) %>%
    mutate(abs.r =  abs(r)) %>%
    filter(abs.r > 0.4) %>%
    filter(p < 0.05) %>%
    relocate(source, target) 
  
  write.csv(edges, file = paste0("Result/FigureS7/",dss,"_edges.csv"),quote = F, row.names = F)
}

# for each edge file, calculate net modules using Gephi
# Manually identify overlapped nodes in modules between healthy and preCOPD nets, and between preCOPD and COPD nets
# 
# nodes and corresponding module colors were saved in : Module_colors.RData

# remove(list = ls())
# load("NodeModuleColors.RData")
# load("microbeAbb.RData")

#2.1 High ---------------------------------------------------------------------

library(openxlsx)
dss = "High" # manually change: Health, preCOPD, COPD
edges <- fread(paste0("Result/FigureS7/",dss,"_edges.csv"), data.table = F) 
#nodeColors <- eval(parse(text = paste0("nodesColors_",dss)))
nodeColors<-read.xlsx(paste0("Result/FigureS7/","nodeColors_GDM.xlsx"), sheet = 1) 

#microbeAbb<-read.xlsx("metadata/Taxa_16S.xlsx",sheet = 1)
microbeAbb<-read.csv("Data/Sequencing/metagenomic/Data/Classes/Species.csv",header=T,sep=",")
microbeAbb<-microbeAbb[,1:7]
microbeAbb$phylum<-gsub("p__","",microbeAbb$Phylum)

internal.node<-unique(c(edges$source,edges$target))

extraNodes<-setdiff(nodeColors$Id,internal.node)

nodeColors<-subset(nodeColors,Id%in%internal.node)

#extraNodes <- microbeAbb$V2[!microbeAbb$V2 %in% as.character(nodeColors$Id )]

selfEdges <- cbind.data.frame(source = extraNodes, 
                              target = extraNodes,
                              edgeType = "selfLink",
                              stringsAsFactors=F) 

edges <- bind_rows(edges, selfEdges)

nodes <- 
  data.frame(table(c(edges$source, edges$target)) ) %>%
  mutate(id = Var1) %>% dplyr::select(-Var1) %>% 
  mutate(abb = sapply(id, function(x) microbeAbb$Genus[which(microbeAbb$Species == x)])) %>%
  mutate(type = sub("_/d+","", abb)) %>%
  mutate(phylum =  sapply(id, function(x) microbeAbb$phylum[which(microbeAbb$Species == x)])) %>%
  mutate(phylum6 = sapply(phylum, 
                          function(x){
                            if(x %in% c("Bacteroidetes","Firmicutes","Proteobacteria","Actinobacteria")) x else "Others"
                          } )) %>%
  relocate(id)  



# plot with igraph ====================
library(igraph)
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
net


V(net)$module.color <- 
  sapply(names(V(net)),
         function(x){
           if(x %in% nodeColors$Id){
             nodeColors$color[which(nodeColors$Id == x)]
           }else "gray"
         }) 

df <- cbind.data.frame(
  V(net)$module.color,
  V(net)$phylum6,
  V(net)$type
)
assign(paste0(dss,"_df"), df, envir = .GlobalEnv)

V(net)$type.color <-
  sapply(1:length(V(net)),
         function(i){
           if(V(net)$type[i] == "bacteria") "white" else V(net)$module.color[i]
         })


E(net)$source.module.color  <- 
  sapply(E(net)$Var1,
         function(x) {
           if(is.na(x) | !x %in% nodeColors$Id) "gray" else nodeColors$color[which(nodeColors$Id == x)]
         })


net.simp <- net - E(net)[E(net)$edgeType=="selfLink"]  #remove self links, only keep the nodes

dir.create("Result/FigureS7")
# Fig5a:
pdf(paste("Result/FigureS7/", dss,".pdf",sep = ""), 
    width = 4 , height = 4) 
#储存的图片大小（对应圆的面积）跟node个数成正比，所以长宽与node个数的平方根成正比
par(mar=c(0.1,0.1,0.1,0.1)) 
set.seed(100)

plot(net.simp, vertex.label=NA, 
     vertex.size = 4,
     vertex.color = V(net)$type.color,
     vertex.frame.color= V(net)$module.color,
     edge.size = 1, 
     edge.color = E(net)$source.module.color, 
     layout=layout_with_fr(net)
)

dev.off()


#2.2 Low -----------------------------------------------------------------


dss = "Low" # manually change: Health, preCOPD, COPD
edges <- fread(paste0("Result/FigureS7/",dss,"_edges.csv"), data.table = F) 
#nodeColors <- eval(parse(text = paste0("nodesColors_",dss)))
nodeColors<-read.xlsx(paste0("Result/FigureS7/","nodeColors_control.xlsx"), sheet = 1) 
microbeAbb<-microbeAbb[,1:7]
microbeAbb$phylum<-gsub("p__","",microbeAbb$Phylum)

internal.node<-unique(c(edges$source,edges$target))

extraNodes<-setdiff(nodeColors$Id,internal.node)

nodeColors<-subset(nodeColors,Id%in%internal.node)

#extraNodes <- microbeAbb$V2[!microbeAbb$V2 %in% as.character(nodeColors$Id )]

selfEdges <- cbind.data.frame(source = extraNodes, 
                              target = extraNodes,
                              edgeType = "selfLink",
                              stringsAsFactors=F) 

edges <- bind_rows(edges, selfEdges)

nodes <- 
  data.frame(table(c(edges$source, edges$target)) ) %>%
  mutate(id = Var1) %>% dplyr::select(-Var1) %>% 
  mutate(abb = sapply(id, function(x) microbeAbb$Species[which(microbeAbb$Species == x)])) %>%
  mutate(type = sub("_/d+","", abb)) %>%
  mutate(phylum =  sapply(id, function(x) microbeAbb$phylum[which(microbeAbb$Species == x)])) %>%
  mutate(phylum6 = sapply(phylum, 
                          function(x){
                            if(x %in% c("Bacteroidetes","Firmicutes","Proteobacteria","Actinobacteria")) x else "Others"
                          } )) %>%
  relocate(id)  



# plot with igraph ====================
library(igraph)
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
net


V(net)$module.color <- 
  sapply(names(V(net)),
         function(x){
           if(x %in% nodeColors$Id){
             nodeColors$color[which(nodeColors$Id == x)]
           }else "gray"
         }) 

df <- cbind.data.frame(
  V(net)$module.color,
  V(net)$phylum6,
  V(net)$type
)
assign(paste0(dss,"_df"), df, envir = .GlobalEnv)

V(net)$type.color <-
  sapply(1:length(V(net)),
         function(i){
           if(V(net)$type[i] == "bacteria") "white" else V(net)$module.color[i]
         })


E(net)$source.module.color  <- 
  sapply(E(net)$Var1,
         function(x) {
           if(is.na(x) | !x %in% nodeColors$Id) "gray" else nodeColors$color[which(nodeColors$Id == x)]
         })


net.simp <- net - E(net)[E(net)$edgeType=="selfLink"]  #remove self links, only keep the nodes

#dir.create("Result/GM_GDM/metagenomic/Network/NetPlots")
# Fig5a:
pdf(paste("Result/FigureS7/", dss,".pdf",sep = ""), 
    width = 4 , height = 4) 
#储存的图片大小（对应圆的面积）跟node个数成正比，所以长宽与node个数的平方根成正比
par(mar=c(0.1,0.1,0.1,0.1)) 
set.seed(100)

plot(net.simp, vertex.label=NA, 
     vertex.size = 4,
     vertex.color = V(net)$type.color,
     vertex.frame.color= V(net)$module.color,
     edge.size = 1, 
     edge.color = E(net)$source.module.color, 
     layout=layout_with_fr(net)
)

dev.off()






