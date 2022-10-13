#####GSE100942
#####8 samples

###load packages
BiocManager::install("biomaRt")
BiocManager::install("FactoMineR")
BiocManager::install("WGCNA")
library(WGCNA)
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(GEOquery)
library(ggplot2)
library(FactoMineR)
library(nlme)

##use function choose.dir() to filter working directory
dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE100942_RAW")##dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE100942_RAW")
setwd(dir)
getwd()
# Select the file with the suffix cel.gz
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
# check the names of files
basename(cel.files)
GSE100942_raw<- ReadAffy(filenames = cel.files)

##The default sample name to read into the chip is the file name, which is
##viewed and modified using the sampleNames function 
sampleNames(GSE100942_raw)
sampleNames(GSE100942_raw)<-stri_sub(sampleNames(GSE100942_raw),1,10)#the length of name maybe 8,9 or 10
sampleNames(GSE100942_raw)
###make the metadata
pData(GSE100942_raw)
GSE100942_pheno=pData(GSE100942_raw)
GSE100942_pheno$sample=rownames(GSE100942_pheno)
GSE100942_pheno$sample_type=rep(c("tumor","normal"),time=4)
GSE100942_pheno$Group=as.numeric(rep(c("1","0"),time=4))
GSE100942_pheno
##save the metadata
save(GSE100942_pheno, file = "GSE100942_pheno.RData", quote = F, row.names = F)

###use RMA function to preprocess the rawdata
###:Background processing, LOG2 transformation,
###quantile scaling and probe expression measurement were performed successively
GSE100942_raw_rma <- rma(GSE100942_raw)
GSE100942 <- exprs(GSE100942_raw_rma)


##clean the metadata
GSE100942_pheno$Group <- as.numeric(GSE100942_pheno$Group)
GSE100942_pheno$sample_type <- as.factor(GSE100942_pheno$sample_type)

####regress  unknown covariates
###special states, GSE100942 with no regress  unknown covariates
mod = model.matrix(~as.factor(Group),data = GSE100942_pheno)
mod0 = model.matrix(~1,data = GSE100942_pheno)
n.sv = num.sv(GSE100942, mod, method="be")
svobj = sva(GSE100942,mod, mod0,n.sv=n.sv)
cor = cor(GSE100942_pheno$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples =8)##check the significance of the correlation
### Regress all Technical and Non-DX biological covariates
X = svobj$sv
Y = GSE100942
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress =  as.matrix(X[,2]) %*% (as.matrix(beta[2,]))

GSE100942 = GSE100942-t(to_regress)



###anno the gene symbol
GSE100942_anno <- as.data.frame(GSE100942)
GSE100942_anno$ID<-rownames(GSE100942_anno)
gpl570<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl570)[c('ID','Gene Symbol')])
GSE100942_anno<-merge(x=GSE100942_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE100942_anno <- GSE100942_anno[,-1]
###The expression levels of the same gene names were averaged
GSE100942_anno<- aggregate(GSE100942_anno,by = list(GSE100942_anno$`Gene Symbol`),FUN = mean)
head(GSE100942_anno)

GSE100942_anno <- GSE100942_anno[-1,]##blank gene name was dropped
rownames(GSE100942_anno) <- GSE100942_anno$Group.1
GSE100942_anno <- GSE100942_anno[,-c(1,10)]



###plots of PCA and boxplot
PCA(GSE100942_anno)
boxplot(GSE100942_anno,main="Boxplot after preprocess")

save(GSE100942_anno,file = "GSE100942_anno.RData",quote = F, row.names = T)
write.csv(GSE100942_anno,file = "GSE100942_anno.csv",quote = F, row.names = T)
