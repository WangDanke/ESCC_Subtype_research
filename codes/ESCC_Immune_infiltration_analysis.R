rm(list=ls())

setwd("D:\\ESCCmRNA-from rawdata\\combined-array")

##loading packages
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(dplyr)

##
load("combat_data_combined.Rdata")
load("data_pheno_merge.Rdata")
##
####
####select the tumor tissue data
a <- data_pheno_merge[which(data_pheno_merge$sample_type=="tumor"),]
##
a <- arrange(a,groups)
##
data_1 <- combat_data_combined[,rownames(a)]
### load the marker genes of each immune cell constitute the background gene set for immune infiltration analysis 
gene_set<- read.csv('gene_reference.csv',
                    header = T)
gene_set<-gene_set[, 1:2]
head(gene_set)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
##ssGSEA
data_1_gsva<- gsva(as.matrix(data_1), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

##pheatmap
library(pheatmap)
data_1_gsva1<- t(scale(t(data_1_gsva)))#
data_1_gsva1[data_1_gsva1< -2] <- -2
data_1_gsva1[data_1_gsva1>2] <- 2
anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
anti<- gsub('^ ','',rownames(data_1_gsva1))%in%anti_tumor
pro<- gsub('^ ','',rownames(data_1_gsva1))%in%pro_tumor
non <- !(anti|pro)##
data_1_gsva1<- rbind(data_1_gsva1[anti,],data_1_gsva1[pro,],data_1_gsva1[non,])#
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}# set the normalization function
nor_data_1_gsva1 <- normalization(data_1_gsva1)
annotation_col = data.frame(group=a$groups)#
rownames(annotation_col)<-colnames(data_1)#
bk = unique(c(seq(0,1, length=100)))# set the pheatmap parameters
pheatmap(nor_data_1_gsva1,
         show_colnames = F,
         cluster_rows = F,cluster_cols = F,
         annotation_col = annotation_col,
         breaks=bk,cellwidth=1,cellheight=2,
         fontsize=2,gaps_row = c(12,20),
         filename = 'ssgsea_all.pdf',width = 8)#
save(data_1_gsva,data_1_gsva1,file = 'immune_infiltration_score.Rdata')




































