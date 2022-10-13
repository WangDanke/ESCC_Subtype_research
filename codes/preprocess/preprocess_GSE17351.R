
#####GSE17351
#####sample size：10

###loading packages
BiocManager::install("biomaRt")
BiocManager::install("FactoMineR")
library(WGCNA)
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(GEOquery)
library(PCA)
library(ggplot2)
library(FactoMineR)

##function "choose.dir（）" to choose folder
dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE17351_RAW")##dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE17351_RAW")
setwd(dir)
getwd()
# list the ".CEL" files
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)
GSE17351_raw<- ReadAffy(filenames = cel.files)#读入文件
sampleNames(GSE17351_raw)

##rename samples
sampleNames(GSE17351_raw)<-stri_sub(sampleNames(GSE17351_raw),1,9)##8或9或10
sampleNames(GSE17351_raw)
###metadata
pData(GSE17351_raw)
GSE17351_pheno=pData(GSE17351_raw)
GSE17351_pheno$sample=rownames(GSE17351_pheno)
GSE17351_pheno$sample_type=rep(c("normal","tumor"),time=5)
GSE17351_pheno$Group=as.numeric(rep(c("0","1"),time=5))
GSE17351_pheno
###save the files
save(GSE17351_pheno,file = "GSE17351_pheno.RData", quote = F, row.names = F)

####function "RNA" to preprocess
GSE17351_raw_rma <- rma(GSE17351_raw)
GSE17351 <- exprs(GSE17351_raw_rma)


##clean the metadata 
GSE17351_pheno$Group <- as.numeric(GSE17351_pheno$Group)
GSE17351_pheno$sample_type <- as.factor(GSE17351_pheno$sample_type)

####regress  unknown covariates 
mod = model.matrix(~as.factor(Group),data = GSE17351_pheno)
mod0 = model.matrix(~1,data = GSE17351_pheno)
n.sv = num.sv(GSE17351, mod, method="be")
svobj = sva(GSE17351,mod, mod0,n.sv=n.sv)
cor = cor(GSE17351_pheno$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 10)

X = svobj$sv
Y = GSE17351
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))

GSE17351 = GSE17351-t(to_regress)



###gene annotation
GSE17351_anno <- as.data.frame(GSE17351)
GSE17351_anno$ID<-rownames(GSE17351_anno)
gpl570<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl570)[c('ID','Gene Symbol')])
GSE17351_anno<-merge(x=GSE17351_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE17351_anno <- GSE17351_anno[,-1]
###The expression levels of genes with the same gene name were averaged
GSE17351_anno<- aggregate(GSE17351_anno,by = list(GSE17351_anno$`Gene Symbol`),FUN = mean)
head(GSE17351_anno)
#Group.1 GSM433786 GSM433787 GSM433788 GSM433789 GSM433790 GSM433791 GSM433792 GSM433793 GSM433794
#1           5.057084  5.015055  4.923569  5.103173  4.908568  4.823771  4.857703  4.858489  4.884393
#2     A1BG  5.969950  5.974494  5.924251  7.201177  6.272978  5.803666  5.890439  6.503881  6.157236
#3 A1BG-AS1  5.304607  5.007162  4.770723  5.618097  5.301334  4.663768  5.073829  5.102902  4.967380
#4     A1CF  4.760489  4.571400  4.829012  4.745067  4.801011  4.558535  4.837814  4.739964  5.018903
#5      A2M  8.270933  7.301832  8.776535  8.291396  8.742367  8.263867  8.628180  8.722209  8.686100
#6  A2M-AS1  5.266746  4.252222  7.447606  6.074372  5.916305  5.997484  5.759067  5.351163  5.499208
GSE17351_anno <- GSE17351_anno[-1,]
rownames(GSE17351_anno) <- GSE17351_anno$Group.1
GSE17351_anno <- GSE17351_anno[,-c(1,12)]



###PCA  and boxplot
PCA(GSE17351_anno)
boxplot(GSE17351_anno,main="Boxplot after preprocess")

save(GSE17351_anno,file = "GSE17351_anno.RData",quote = F, row.names = T)
write.csv(GSE17351_anno,file = "GSE17351_anno.csv",quote = F, row.names = T)




