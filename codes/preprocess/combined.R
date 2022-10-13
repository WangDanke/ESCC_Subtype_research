options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("limma")) BiocManager::install("limma",update = F,ask = F)
if(!require("bladderbatch")) BiocManager::install("bladderbatch",update = F,ask = F)
if(!require("sva")) BiocManager::install("sva",update = F,ask = F)
if(!require("preprocessCore")) BiocManager::install("preprocessCore",update = F,ask = F)
if(!require("tidyr")) install.packages("tidyr",update = F,ask = F)
if(!require("dplyr")) install.packages("dplyr",update = F,ask = F)
if(!require("tibble")) install.packages("tibble",update = F,ask = F)
if(!require("factoextra")) install.packages("factoextra",update = F,ask = F)
if(!require("FactoMineR")) install.packages("FactoMineR",update = F,ask = F)
if(!require("ggplot2")) install.packages("ggplot2",update = F,ask = F)
if(!require("ggrepel")) install.packages("ggrepel",update = F,ask = F)
if(!require("ggpubr")) install.packages("ggpubr",update = F,ask = F)
if(!require("dendextend")) install.packages("dendextend",update = F,ask = F)
if(!require("RobustRankAggreg")) install.packages("RobustRankAggreg",update = F,ask = F)
#load the packages that we needed

dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\combined-array")
setwd(dir)
getwd()

###add batch in the pheno data
GSE100942_pheno$batch <- "GSE100942"
GSE161533_pheno$batch <- "GSE161533"
GSE17351_pheno$batch <- "GSE17351"
GSE20347_pheno$batch <- "GSE20347"
GSE23400_pheno$batch <- "GSE23400"                             
GSE38129_pheno$batch <- "GSE38129"
GSE77861_pheno$batch <- "GSE77861"
##combined the pheno data of different batches
data_pheno <- rbind(GSE100942_pheno,GSE161533_pheno,GSE17351_pheno,GSE20347_pheno,
                    GSE23400_pheno,GSE38129_pheno,GSE77861_pheno)
##find the shared genes in diffenent batches
shared_gene<-Reduce(intersect,list(rownames(GSE100942_anno),
                          rownames(GSE161533_anno),
                          rownames(GSE17351_anno),
                          rownames(GSE20347_anno),
                          rownames(GSE23400_anno),
                          rownames(GSE38129_anno),
                          rownames(GSE77861_anno)))
###combined the different batches by the shared genes
GSE100942_anno_shar <- GSE100942_anno[shared_gene,]
GSE161533_anno_shar <- GSE161533_anno[shared_gene,]
GSE17351_anno_shar <- GSE17351_anno[shared_gene,]
GSE20347_anno_shar <- GSE20347_anno[shared_gene,]
GSE23400_anno_shar <- GSE23400_anno[shared_gene,]
GSE38129_anno_shar <- GSE38129_anno[shared_gene,]
GSE77861_anno_shar <- GSE77861_anno[shared_gene,]

data_combined <- cbind(GSE100942_anno_shar,GSE161533_anno_shar,GSE17351_anno_shar,GSE20347_anno_shar,
                       GSE23400_anno_shar,GSE38129_anno_shar,GSE77861_anno_shar)

#save the data
save(data_pheno,file = "data_pheno.RData",quote = F, row.names = T)
write.csv(data_pheno,file = "data_pheno.csv",quote = F, row.names = T)
save(data_combined,file = "data_combined.RData",quote = F, row.names = T)
write.csv(data_combined,file = "data_combined.csv",quote = F, row.names = T)

#####batch correct
#show tne distribution of the data
data_pheno$batch <- as.factor(data_pheno$batch)
#boxplot
boxplot(data_combined,col = data_pheno$batch, main="Expression before_batchcor")
#PCA  
ddb.pca <- PCA(t(data_combined), graph = FALSE)
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = data_pheno$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
###batch correction use the sva package
library(sva)
model <- model.matrix(~as.factor(data_pheno$sample_type))
combat_data_combined <- ComBat(dat = data_combined, batch = data_pheno$batch, mod = model)
#save the data
save(combat_data_combined,file = "combat_data_combined.RData",quote = F, row.names = T)
write.csv(combat_data_combined,file = "combat_data_combined.csv",quote = F, row.names = T)

#show tne distribution of the data after batch correction
data_pheno_merge$groups <- as.factor(data_pheno_merge$groups)
#boxplot
boxplot(combat_data_combined,col = data_pheno$batch, main="Expression aft_batchcorrected")
pdf("boxplot-aft-batchcor.pdf")
dev.off()
#PCA  
ddb.pca <- PCA(t(combat_data_combined_1), graph = FALSE)
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = a$groups, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "groups",
             range=0
)



