####Validata these three Subtypes with RNA-seq data from ESCC based on the TCGA database
dir = choose.dir(default = "",caption = "")
setwd(dir)
getwd()
####Use the 80 ESCC tumor samples to get three Subtypes  
mydata <- t(TCGA_ESCC_RNA_seq_80_TPM[rownames(combat_data_combined),])
d <- dist(mydata, method = "euclidean") 
fit2 <- hclust(d, method="ward.D")

plot(fit2)
hcd = as.dendrogram(fit2)
##
groups <- cutree(fit2, k=3)
###plot the clustering tree
# load code of A2R function
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
# colored dendrogram
op = par(bg = "#FFFFFF")
A2Rplot(fit2, k = 3, boxes = FALSE, col.up = "gray50", col.down = c("#999999", "#E69F00", "#56B4E9"),
        main = "clustering of TCGA_80 ESCC samples")


####Correlation analysis was performed between the three TCGA_Subtypes obtained based on RNA-seq data and the three Subtypes obtained based on RNA-microarray data
##################################################################################
GSE_pheno <- subset(data_pheno_merge,data_pheno_merge$sample_type=="tumor")
GSE_data <- combat_data_combined[,rownames(GSE_pheno)]

TCGA_data <- log2(TCGA_ESCC_RNA_seq_80_TPM)
TCGA_data <- TCGA_data[!is.infinite(rowSums(TCGA_data)),]
AA <- intersect(rownames(GSE_data),rownames(TCGA_data))


library(Hmisc)
ESCC_221_cor <- cor(GSE_data[AA,],TCGA_data[AA,],method = "spearman")
TCGA_subtype1 <- subset(TCGA_ESCC_survive_pheno_80,TCGA_ESCC_survive_pheno_80$subtypes=="TCGA_Subtype1")
ESCC_cor_subtype1_TCGA <- ESCC_221_cor[,rownames(TCGA_subtype1)]
ESCC_cor_subtype1_TCGA$mean_cor <- apply(ESCC_cor_subtype1_TCGA,1,mean)
ESCC_cor_subtype1_TCGA$group <- GSE_pheno$SUBTYPE

save(ESCC_cor_subtype1_TCGA,ESCC_cor_subtype2_TCGA,ESCC_cor_subtype3_TCGA,file = "TCGA_Subtypes_cor_GEO_Subtypes.RData")

##Figure 8 B,C,D
##violin plot shows the  distribution of mean correlation values among subtypes
library(ggpubr)
my_comparisons <- list(c('subtype1','subtype2'),
                       c('subtype2','subtype3'),c('subtype1','subtype3'))
ggviolin(ESCC_cor_subtype3_TCGA,x="group",y="mean_cor",fill='group', 
         palette = c("#999999", "#E69F00", "#56B4E9"),
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = my_comparisons)


#### Survival analysis between TCGA_Subtypes
library(survival)
library(survminer)
library(preprocessCore)
sfit <- survfit(Surv(time, event)~subtypes, data=TCGA_ESCC_survive_pheno_80)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)



















