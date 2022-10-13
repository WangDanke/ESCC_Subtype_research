#####GSE20347
#####34samples

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
dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE20347_RAW")##dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE20347_RAW")
setwd(dir)
getwd()
# Select the file with the suffix cel.gz
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
# check the names of files
basename(cel.files)
GSE20347_raw<- ReadAffy(filenames = cel.files)

##The default sample name to read into the chip is the file name, which is
##viewed and modified using the sampleNames function 
sampleNames(GSE20347_raw)
sampleNames(GSE20347_raw)<-stri_sub(sampleNames(GSE20347_raw),1,9)#the length of name maybe 8,9 or 10
sampleNames(GSE20347_raw)
###make the metadata
pData(GSE20347_raw)
GSE20347_pheno=pData(GSE20347_raw)
GSE20347_pheno$sample=rownames(GSE20347_pheno)
GSE20347_pheno$sample_type=c(rep(c("normal"),time=17),rep(c("tumor"),time=17))
GSE20347_pheno$Group=as.numeric(c(rep(c("0"),time=17),rep(c("1"),time=17)))
GSE20347_pheno
##save the metadata
save(GSE20347_pheno, file = "GSE20347_pheno.RData", quote = F, row.names = F)

###use RMA function to preprocess the rawdata
###:Background processing, LOG2 transformation,
###quantile scaling and probe expression measurement were performed successively
GSE20347_raw_rma <- rma(GSE20347_raw)
GSE20347 <- exprs(GSE20347_raw_rma)


##clean the metadata
GSE20347_pheno$Group <- as.numeric(GSE20347_pheno$Group)
GSE20347_pheno$sample_type <- as.factor(GSE20347_pheno$sample_type)

####regress  unknown covariates 
mod = model.matrix(~as.factor(Group),data = GSE20347_pheno)
mod0 = model.matrix(~1,data = GSE20347_pheno)
n.sv = num.sv(GSE20347, mod, method="be")
svobj = sva(GSE20347,mod, mod0,n.sv=n.sv)
cor = cor(GSE20347_pheno$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 34)##check the significance of the correlation
### Regress all Technical and Non-DX biological covariates
X = svobj$sv
Y = GSE20347
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))

GSE20347 = GSE20347-t(to_regress)



###anno the gene symbol
GSE20347_anno <- as.data.frame(GSE20347)
GSE20347_anno$ID<-rownames(GSE20347_anno)
gpl571<-getGEO("GPL571",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE20347_anno<-merge(x=GSE20347_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE20347_anno <- GSE20347_anno[,-1]
###The expression levels of the same gene names were averaged
GSE20347_anno<- aggregate(GSE20347_anno,by = list(GSE20347_anno$`Gene Symbol`),FUN = mean)
head(GSE20347_anno)
#Group.1 GSM433786 GSM433787 GSM433788 GSM433789 GSM433790 GSM433791 GSM433792 GSM433793 GSM433794
#1           5.057084  5.015055  4.923569  5.103173  4.908568  4.823771  4.857703  4.858489  4.884393
#2     A1BG  5.969950  5.974494  5.924251  7.201177  6.272978  5.803666  5.890439  6.503881  6.157236
#3 A1BG-AS1  5.304607  5.007162  4.770723  5.618097  5.301334  4.663768  5.073829  5.102902  4.967380
#4     A1CF  4.760489  4.571400  4.829012  4.745067  4.801011  4.558535  4.837814  4.739964  5.018903
#5      A2M  8.270933  7.301832  8.776535  8.291396  8.742367  8.263867  8.628180  8.722209  8.686100
#6  A2M-AS1  5.266746  4.252222  7.447606  6.074372  5.916305  5.997484  5.759067  5.351163  5.499208
GSE20347_anno <- GSE20347_anno[-1,]##blank gene name was dropped
rownames(GSE20347_anno) <- GSE20347_anno$Group.1
GSE20347_anno <- GSE20347_anno[,-c(1,36)]



###plots of PCA and boxplot
PCA(GSE20347_anno)
boxplot(GSE20347_anno,main="Boxplot after preprocess")

save(GSE20347_anno,file = "GSE20347_anno.RData",quote = F, row.names = T)
write.csv(GSE20347_anno,file = "GSE20347_anno.csv",quote = F, row.names = T)

















