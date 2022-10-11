#####GSE38129
#####60 samples

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
dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE38129_RAW")##dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE38129_RAW")
setwd(dir)
getwd()
# Select the file with the suffix cel.gz
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
# check the names of files
basename(cel.files)
GSE38129_raw<- ReadAffy(filenames = cel.files)

##The default sample name to read into the chip is the file name, which is
##viewed and modified using the sampleNames function 
sampleNames(GSE38129_raw)
sampleNames(GSE38129_raw)<-stri_sub(sampleNames(GSE38129_raw),1,9)#the length of name maybe 8,9 or 10
sampleNames(GSE38129_raw)
###make the metadata
pData(GSE38129_raw)
GSE38129_pheno=pData(GSE38129_raw)
GSE38129_pheno$sample=rownames(GSE38129_pheno)
GSE38129_pheno$sample_type=rep(c("normal","tumor"),time=30)
GSE38129_pheno$Group=as.numeric(rep(c("0","1"),time=30))
GSE38129_pheno
##save the metadata
save(GSE38129_pheno, file = "GSE38129_pheno.RData", quote = F, row.names = F)

###use RMA function to preprocess the rawdata
###:Background processing, LOG2 transformation,
###quantile scaling and probe expression measurement were performed successively
GSE38129_raw_rma <- rma(GSE38129_raw)
GSE38129 <- exprs(GSE38129_raw_rma)


##clean the metadata
GSE38129_pheno$Group <- as.numeric(GSE38129_pheno$Group)
GSE38129_pheno$sample_type <- as.factor(GSE38129_pheno$sample_type)

####regress  unknown covariates 
mod = model.matrix(~as.factor(Group),data = GSE38129_pheno)
mod0 = model.matrix(~1,data = GSE38129_pheno)
n.sv = num.sv(GSE38129, mod, method="be")
svobj = sva(GSE38129,mod, mod0,n.sv=n.sv)
cor = cor(GSE38129_pheno$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 60)##check the significance of the correlation
### Regress all Technical and Non-DX biological covariates
X = svobj$sv
Y = GSE38129
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))

GSE38129 = GSE38129-t(to_regress)



###anno the gene symbol
GSE38129_anno <- as.data.frame(GSE38129)
GSE38129_anno$ID<-rownames(GSE38129_anno)
gpl571<-getGEO("GPL571",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE38129_anno<-merge(x=GSE38129_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE38129_anno <- GSE38129_anno[,-1]
###The expression levels of the same gene names were averaged
GSE38129_anno<- aggregate(GSE38129_anno,by = list(GSE38129_anno$`Gene Symbol`),FUN = mean)
head(GSE38129_anno)

GSE38129_anno <- GSE38129_anno[-1,]##blank gene name was dropped
rownames(GSE38129_anno) <- GSE38129_anno$Group.1
GSE38129_anno <- GSE38129_anno[,-c(1,62)]



###plots of PCA and boxplot
PCA(GSE38129_anno)
boxplot(GSE38129_anno,main="Boxplot after preprocess")

save(GSE38129_anno,file = "GSE38129_anno.RData",quote = F, row.names = T)
write.csv(GSE38129_anno,file = "GSE38129_anno.csv",quote = F, row.names = T)








