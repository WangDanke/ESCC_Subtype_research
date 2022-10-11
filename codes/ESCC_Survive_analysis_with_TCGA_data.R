
####single gene survival analysis
#######
######## 
dir <- choose.dir(default = "", caption = "")
setwd(dir)
getwd()
BiocManager::install("survival")
BiocManager::install("survminer")
library(survival)
library(survminer)
library(preprocessCore)
##vital_status,death equals to "1".
TCGA_ESCC_survive_pheno_80$days_to_death[is.na(TCGA_ESCC_survive_pheno_80$days_to_death)] <- 0   
TCGA_ESCC_survive_pheno_80$days_to_last_follow_up[is.na(TCGA_ESCC_survive_pheno_80$days_to_last_follow_up)] <- 0
TCGA_ESCC_survive_pheno_80$time=as.numeric(TCGA_ESCC_survive_pheno_80[,1])+as.numeric(TCGA_ESCC_survive_pheno_80[,3])
TCGA_ESCC_survive_pheno_80$event=ifelse(TCGA_ESCC_survive_pheno_80$vital_status=='Alive',0,1)
save(TCGA_ESCC_survive_pheno_80,file = "TCGA_ESCC_survive_80.RData")

##log transform
data_1 <- TCGA_ESCC_RNA_seq_80_TPM
data_1 <- log2(data_1+1)

###survival analysis of genes
MYOC<- data_1[c('MYOC'),]
TCGA_ESCC_survive_pheno_80 <- cbind(TCGA_ESCC_survive_pheno_80,as.data.frame(t(MYOC)))
TCGA_ESCC_survive_pheno_80$Group_MYOC<- ifelse(TCGA_ESCC_survive_pheno_80$MYOC>median(TCGA_ESCC_survive_pheno_80$MYOC),"High","Low")
##Figure 6
sfit <- survfit(Surv(time, event)~Group_MYOC, data=TCGA_ESCC_survive_pheno_80)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)














