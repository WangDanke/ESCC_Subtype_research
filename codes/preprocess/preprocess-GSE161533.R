#####GSE161533
#####56 samples

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
dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE161533_RAW")##dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE161533_RAW")
setwd(dir)
getwd()
# Select the file with the suffix cel.gz
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
# check the names of files
basename(cel.files)
GSE161533_raw<- ReadAffy(filenames = cel.files)

##The default sample name to read into the chip is the file name, which is
##viewed and modified using the sampleNames function 
sampleNames(GSE161533_raw)
sampleNames(GSE161533_raw)<-stri_sub(sampleNames(GSE161533_raw),1,10)#the length of name maybe 8,9 or 10
sampleNames(GSE161533_raw)
###make the metadata
pData(GSE161533_raw)
GSE161533_pheno=pData(GSE161533_raw)
GSE161533_pheno$sample=rownames(GSE161533_pheno)
GSE161533_pheno$sample_type=c(rep(c("normal"),time=28),rep(c("tumor"),time=28))
GSE161533_pheno$Group=as.numeric(c(rep(c("0"),time=28),rep(c("1"),time=28)))
GSE161533_pheno
##save the metadata
save(GSE161533_pheno, file = "GSE161533_pheno.RData", quote = F, row.names = F)

###use RMA function to preprocess the rawdata
###:Background processing, LOG2 transformation,
###quantile scaling and probe expression measurement were performed successively
GSE161533_raw_rma <- rma(GSE161533_raw)
GSE161533 <- exprs(GSE161533_raw_rma)

##clean the metadata
GSE161533_pheno$Group <- as.numeric(GSE161533_pheno$Group)
GSE161533_pheno$sample_type <- as.factor(GSE161533_pheno$sample_type)

####regress  unknown covariates 
mod = model.matrix(~as.factor(Group),data = GSE161533_pheno)
mod0 = model.matrix(~1,data = GSE161533_pheno)
n.sv = num.sv(GSE161533, mod, method="be")
svobj = sva(GSE161533,mod, mod0,n.sv=n.sv)
cor = cor(GSE161533_pheno$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples =56)##check the significance of the correlation
### Regress all Technical and Non-DX biological covariates
X = svobj$sv
Y = GSE161533
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,-c(7)]) %*% (as.matrix(beta[-c(7),])))

GSE161533 = GSE161533-t(to_regress)



###anno the gene symbol
GSE161533_anno <- as.data.frame(GSE161533)
GSE161533_anno$ID<-rownames(GSE161533_anno)
gpl570<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl570)[c('ID','Gene Symbol')])
GSE161533_anno<-merge(x=GSE161533_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE161533_anno <- GSE161533_anno[,-1]
###The expression levels of the same gene names were averaged
GSE161533_anno<- aggregate(GSE161533_anno,by = list(GSE161533_anno$`Gene Symbol`),FUN = mean)
head(GSE161533_anno)

GSE161533_anno <- GSE161533_anno[-1,]##blank gene name was dropped
rownames(GSE161533_anno) <- GSE161533_anno$Group.1
GSE161533_anno <- GSE161533_anno[,-c(1,58)]



###plots of PCA and boxplot
PCA(GSE161533_anno)
boxplot(GSE161533_anno,main="Boxplot after preprocess")

save(GSE161533_anno,file = "GSE161533_anno.RData",quote = F, row.names = T)
write.csv(GSE161533_anno,file = "GSE161533_anno.csv",quote = F, row.names = T)
