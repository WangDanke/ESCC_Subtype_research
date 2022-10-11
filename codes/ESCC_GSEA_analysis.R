##set the work enviroment
dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\combined-array")
setwd(dir)
getwd()
##loading packages
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("clusterProfiler")) BiocManager::install("clusterProfiler",update = F,ask = F)
if(!require("enrichplot")) BiocManager::install("enrichplot",update = F,ask = F)
if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db",update = F,ask = F)

###three subtypes are subjective to GSEA
dif <- allDiff_subtype3[,1]
dif <- as.data.frame(dif)
dif$SYMBOL <- rownames(allDiff_subtype3)
colnames(dif)[1] <- "logFC"

dif_ID <- bitr(dif$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db")

dif_all <- merge(dif,dif_ID,all=F)
write.csv(dif_all,file="dif_all.csv",quote = F, row.names = T)
dim(dif_all)
dif_all_order <- dif_all[order(dif_all$logFC, decreasing = T),]
gene_logFC <- dif_all_order$logFC
head(gene_logFC)
names(gene_logFC) <- dif_all_order$ENTREZID
head(gene_logFC)


###one step GSEA
GSEA_GO_subtype3 <- gseGO(gene_logFC,
                 ont = "ALL",
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
##save the results
save(GSEA_GO_subtype1,file = "GSEA_GO_subtype1.RData",quote = T,row.names = T)
save(GSEA_GO_subtype2,file = "GSEA_GO_subtype2.RData",quote = T,row.names = T)
save(GSEA_GO_subtype3,file = "GSEA_GO_subtype3.RData",quote = T,row.names = T)

###plot (figure3)
###
BiocManager::install("forcats")
BiocManager::install("ggstance")
BiocManager::install("ggupset")
library(forcats)
library(ggstance)
library(ggupset)
library(pheatmap)

###find the common enriched pathways
a <- rownames(GSEA_GO_subtype1@result)
b <- rownames(GSEA_GO_subtype2@result)
c <- intersect(a,b)
d <- rownames(GSEA_GO_subtype3@result)
com <- intersect(c,d)
GSEA_GO_subtype1_com <- GSEA_GO_subtype1@result[com,]
GSEA_GO_subtype2_com <- GSEA_GO_subtype2@result[com,]
GSEA_GO_subtype3_com <- GSEA_GO_subtype3@result[com,]
GSEA_GO_com <- cbind(GSEA_GO_subtype1_com[,c(5,6)],GSEA_GO_subtype2_com[,c(5,6)],GSEA_GO_subtype3_com[,c(5,6)] )

colnames(GSEA_GO_com) <- c("suntype1_enrichmentScore","suntype1_NES","suntype2_enrichmentScore","suntype2_NES",
                           "suntype3_enrichmentScore","suntype3_NES")
GSEA_GO_com[,c(7,8,9)] <- GSEA_GO_group1_com[,c(1,2,3)]

write.csv(GSEA_GO_com,file = "GSEA_GO_com.csv")
save(GSEA_GO_com,file = "GSEA_GO_com.RData")


GSEA_GO_com_BP <-subset(GSEA_GO_com,GSEA_GO_com$ONTOLOGY =="BP")
GSEA_GO_com_CC <-subset(GSEA_GO_com,GSEA_GO_com_up$ONTOLOGY =="CC")
GSEA_GO_com_MF <-subset(GSEA_GO_com,GSEA_GO_com_up$ONTOLOGY =="MF")
rownames(GSEA_GO_com_BP) <- GSEA_GO_com_BP$Description
##select 20 pathways to demonstrate
GSEA_GO_com_20 <- GSEA_GO_com_BP[c('mitotic DNA replication','hemidesmosome assembly','endodermal cell differentiation',
                                   'collagen catabolic process','collagen fibril organization','extracellular matrix disassembly',
                                   'type I interferon signaling pathway','response to type I interferon','positive regulation of cytokinesis',
                                   'regulation of cell cycle phase transition','very long-chain fatty acid metabolic process',
                                   'epoxygenase P450 pathway','multivesicular body assembly','neurotransmitter metabolic process',
                                   'drug metabolic process','fatty acid beta-oxidation','long-chain fatty acid metabolic process',
                                   'fatty acid oxidation','unsaturated fatty acid metabolic process','fatty acid derivative metabolic process',
                                   'monocarboxylic acid catabolic process'),]
#
rownames(GSEA_GO_com_20) <- paste(GSEA_GO_com_20$ONTOLOGY,GSEA_GO_com_20$Description,sep = ":")
GSEA_GO_com_20 <- GSEA_GO_com_20[,c(1,3,5)]

###figure 3A
BiocManager::install("circlize")
BiocManager::install("ComplexHeatmap")
BiocManager::install("RColorBrewer")
##colorRamp2
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
##
GSEA_GO_com_20 <-GSEA_GO_com_20[order(GSEA_GO_com_20$suntype1_enrichmentScore,decreasing=T),]
##pheatmap
col=colorRamp2(c(1,0.5, 0, -0.5, -1), brewer.pal(n=5, name='RdBu'))
pheatmap(GSEA_GO_com_20,
         color = col,
         show_colnames = T,
         cluster_rows = F,cluster_cols = F)

###Three subtype-specific BP pathways were found

h <- rownames(GSEA_GO_group3@result)
g <- rownames(GSEA_GO_com)
g_special_group3 <- setdiff(h,g)
GSEA_special_group3 <- GSEA_GO_group3@result[g_special_group3,]
GSEA_special_group3_BP <- subset(GSEA_special_group3, GSEA_special_group3$ONTOLOGY =="BP")


###figure 3B-D
library(pheatmap)
library(dotplot)
library(clusterProfiler)

GSEA_special_group1_BP <- GSEA_special_group1_BP[-c(9),]
GSEA_special_group2_BP <- GSEA_special_group2_BP[-c(58),]
GSEA_GO_group1@result <- GSEA_special_group1_BP
GSEA_GO_group2@result <- GSEA_special_group2_BP
GSEA_GO_group3@result <- GSEA_special_group3_BP
##paopao
dotplot(GSEA_GO_group3,split=".sign",showCategory=15)+facet_grid(~.sign)


write.csv(GSEA_special_group3_BP,file = "GSEA_special_group3_BP.csv")
write.csv(GSEA_special_group2_BP,file = "GSEA_special_group2_BP.csv")
write.csv(GSEA_special_group1_BP,file = "GSEA_special_group1_BP.csv")




















