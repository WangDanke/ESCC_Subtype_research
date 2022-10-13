#####GSE23400
#####100 samples

###load packages
BiocManager::install("biomaRt")
BiocManager::install("FactoMineR")
BiocManager::install("WGCNA")
BiocManager::install("oligo")

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
dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE23400_RAW")##dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE23400_RAW")
setwd(dir)
getwd()
# Select the file with the suffix cel.gz
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
# check the names of files
basename(cel.files)
GSE23400_raw<- ReadAffy(filenames = cel.files)


##The default sample name to read into the chip is the file name, which is
##viewed and modified using the sampleNames function 
sampleNames(GSE23400_raw)
sampleNames(GSE23400_raw)<-stri_sub(sampleNames(GSE23400_raw),1,9)#the length of name maybe 8,9 or 10
sampleNames(GSE23400_raw)
###make the metadata
pData(GSE23400_raw)
GSE23400_pheno=pData(GSE23400_raw)
GSE23400_pheno$sample=rownames(GSE23400_pheno)
GSE23400_pheno$sample_type=c(rep(c("normal"),time=53),rep(c("tumor"),time=53))
GSE23400_pheno$Group=as.numeric(c(rep(c("0"),time=53),rep(c("1"),time=53)))
GSE23400_pheno
##save the metadata
save(GSE23400_pheno, file = "GSE23400_pheno.RData", quote = F, row.names = F)

###use RMA function to preprocess the rawdata
###:Background processing, LOG2 transformation,
###quantile scaling and probe expression measurement were performed successively
GSE23400_raw_rma <- rma(GSE23400_raw)
GSE23400 <- exprs(GSE23400_raw_rma)


###density plot
i=1;
plot(density((GSE23400[,1]),na.rm=T), col=as.numeric(GSE23400_pheno)[1],
     main="hist of log2 EXp", xlab="log2 exp", ylim=c(0,1))
for (i in 2:dim(GSE23400)[2]) 
  lines(density((GSE23400[,i]), na.rm=T),col=as.numeric(GSE23400_pheno)[i])
i = 1; 
plot(density((datExpr.prenorm[,i]), na.rm=T), col = as.numeric(datMeta$A.C)[i], main="Hist of Log2 Exp", xlab = "log2 exp", ylim=c(0,0.8))
for(i in 2:dim(datExpr.prenorm)[2])
  lines(density((datExpr.prenorm[,i]), na.rm=T), col = as.numeric(datMeta$A.C)[i])

##clean the metadata
GSE23400_pheno$Group <- as.numeric(GSE23400_pheno$Group)
GSE23400_pheno$sample_type <- as.factor(GSE23400_pheno$sample_type)

####regress  unknown covariates 
mod = model.matrix(~as.factor(Group),data = GSE23400_pheno)
mod0 = model.matrix(~1,data = GSE23400_pheno)
n.sv = num.sv(GSE23400, mod, method="be")
svobj = sva(GSE23400,mod, mod0,n.sv=n.sv)
cor = cor(GSE23400_pheno$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples =106)##check the significance of the correlation
### Regress all Technical and Non-DX biological covariates
X = svobj$sv
Y = GSE23400
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))

GSE23400 = GSE23400-t(to_regress)



###anno the gene symbol
GSE23400_anno <- as.data.frame(GSE23400)
GSE23400_anno$ID<-rownames(GSE23400_anno)
gpl96<-getGEO("GPL96",destdir=".")
iddata<-as.data.frame(Table(gpl96)[c('ID','Gene Symbol')])
GSE23400_anno<-merge(x=GSE23400_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE23400_anno <- GSE23400_anno[,-1]
###The expression levels of the same gene names were averaged
GSE23400_anno<- aggregate(GSE23400_anno,by = list(GSE23400_anno$`Gene Symbol`),FUN = mean)
head(GSE23400_anno)

GSE23400_anno <- GSE23400_anno[-1,]##blank gene name was dropped
rownames(GSE23400_anno) <- GSE23400_anno$Group.1
GSE23400_anno <- GSE23400_anno[,-c(1,108)]



###plots of PCA and boxplot
PCA(GSE23400_anno)
ddb.pca <- PCA(t(GSE23400_anno), graph = FALSE)
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = GSE23400_pheno$sample_type, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)
##remove the outliers
##remove the normal type samples: GSM573851,GSM573889, GSM573852 and the match three  tumor samples
GSE23400_anno <- GSE23400_anno[,-c(6,7,44,59,60,97)]
GSE23400_pheno <- GSE23400_pheno[-c(6,7,44,59,60,97),]

boxplot(GSE23400_anno,main="Boxplot after preprocess")

save(GSE23400_pheno, file = "GSE23400_pheno.RData", quote = F, row.names = F)
save(GSE23400_anno,file = "GSE23400_anno.RData",quote = F, row.names = T)
write.csv(GSE23400_anno,file = "GSE23400_anno.csv",quote = F, row.names = T)





















