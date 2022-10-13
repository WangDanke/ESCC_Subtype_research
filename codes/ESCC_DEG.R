dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\combined-array")
setwd(dir)
getwd()
##load packages
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

#check the structure of the data
str(combat_data_combined) 
head(combat_data_combined) 
dim(combat_data_combined) 

###use the limma package to find the DEGs
data_pheno_merge <- read.csv("data_pheno_merge.csv",row.names=1)
save(data_pheno_merge,file = "data_pheno_merge.Rdata",quote = F, row.names = T)
subtype_1 <- data_pheno_merge[which(data_pheno_merge$SUBTYPE=="subtype1"),]
subtype_2 <- data_pheno_merge[which(data_pheno_merge$SUBTYPE=="subtype2"),]
subtype_3 <- data_pheno_merge[which(data_pheno_merge$SUBTYPE=="subtype3"),]
#### Differential expression analysis of individual subtype and paired normal tissues
data<-combat_data_combined[,rownames(subtype_3)]
group <- as.character(subtype_3$sample_type)
group <- factor(group,levels = c("normal","tumor"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_subtype3=topTable(fit2,adjust='fdr',coef=2,number=Inf)


save(allDiff_subtype1,file = "allDiff_subtype1.Rdata",quote = F, row.names = T)
save(allDiff_subtype2,file = "allDiff_subtype2.Rdata",quote = F, row.names = T)
save(allDiff_subtype3,file = "allDiff_subtype3.Rdata",quote = F, row.names = T)

####
####Fig.S2
####plot the volcano
allDiff_subtype3$result = as.factor(ifelse(allDiff_subtype3$adj.P.Val < 0.05 & abs(allDiff_subtype3$logFC) > 1,
                                         ifelse(allDiff_subtype3$logFC > 1 ,'UP','DOWN'),'NOT')
)

this_tile <- paste0( 'subtype3','\nThe number of up gene is ',nrow(allDiff_subtype3[allDiff_subtype3$result =='UP',]) ,
                     '\nThe number of down gene is ',nrow(allDiff_subtype3[allDiff_subtype3$result =='DOWN',])
)

ggplot(data=allDiff_subtype3, aes(x=logFC, y=-log10(adj.P.Val), color=result)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))+
  theme_bw()













