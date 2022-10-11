#####GSE77861
#####14 samples

###load packages
BiocManager::install("biomaRt")
BiocManager::install("FactoMineR")
library(WGCNA)
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(GEOquery)
library(ggplot2)
library(FactoMineR)

##use function choose.dir() to filter working directory
dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE77861_RAW")##dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE77861_RAW")
setwd(dir)
getwd()
# Select the file with the suffix cel.gz
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
# check the names of files
basename(cel.files)
GSE77861_raw<- ReadAffy(filenames = cel.files)

##The default sample name to read into the chip is the file name, which is
##viewed and modified using the sampleNames function 
sampleNames(GSE77861_raw)
sampleNames(GSE77861_raw)<-stri_sub(sampleNames(GSE77861_raw),1,10)#the length of name maybe 8,9 or 10
sampleNames(GSE77861_raw)
###make the metadata
pData(GSE77861_raw)
GSE77861_pheno=pData(GSE77861_raw)
GSE77861_pheno$sample=rownames(GSE77861_pheno)
GSE77861_pheno$sample_type=rep(c("normal","tumor"),time=7)
GSE77861_pheno$Group=as.numeric(rep(c("0","1"),time=7))
GSE77861_pheno
##save the metadata
save(GSE77861_pheno, file = "GSE77861_pheno.RData", quote = F, row.names = F)

###use RMA function to preprocess the rawdata
###:Background processing, LOG2 transformation,
###quantile scaling and probe expression measurement were performed successively
GSE77861_raw_rma <- rma(GSE77861_raw)
GSE77861 <- exprs(GSE77861_raw_rma)


##clean the metadata
GSE77861_pheno$Group <- as.numeric(GSE77861_pheno$Group)
GSE77861_pheno$sample_type <- as.factor(GSE77861_pheno$sample_type)

####regress  unknown covariates 
mod = model.matrix(~as.factor(Group),data = GSE77861_pheno)
mod0 = model.matrix(~1,data = GSE77861_pheno)
n.sv = num.sv(GSE77861, mod, method="be")
svobj = sva(GSE77861,mod, mod0,n.sv=n.sv)
cor = cor(GSE77861_pheno$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 14)##check the significance of the correlation
### Regress all Technical and Non-DX biological covariates
X = svobj$sv
Y = GSE77861
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))

GSE77861 = GSE77861-t(to_regress)



###anno the gene symbol
GSE77861_anno <- as.data.frame(GSE77861)
GSE77861_anno$ID<-rownames(GSE77861_anno)
gpl570<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl570)[c('ID','Gene Symbol')])
GSE77861_anno<-merge(x=GSE77861_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE77861_anno <- GSE77861_anno[,-1]
###The expression levels of the same gene names were averaged
GSE77861_anno<- aggregate(GSE77861_anno,by = list(GSE77861_anno$`Gene Symbol`),FUN = mean)
head(GSE77861_anno)

GSE77861_anno <- GSE77861_anno[-1,]##blank gene name was dropped
rownames(GSE77861_anno) <- GSE77861_anno$Group.1
GSE77861_anno <- GSE77861_anno[,-c(1,16)]



###plots of PCA and boxplot
PCA(GSE77861_anno)
boxplot(GSE77861_anno,main="Boxplot after preprocess")

save(GSE77861_anno,file = "GSE77861_anno.RData",quote = F, row.names = T)
write.csv(GSE77861_anno,file = "GSE77861_anno.csv",quote = F, row.names = T)





