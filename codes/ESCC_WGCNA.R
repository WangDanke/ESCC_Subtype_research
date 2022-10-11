dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\combined-array")
setwd(dir)
getwd()
###packages
install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)
library('WGCNA')



####select suitable sft
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
data <- as.data.frame(t(combat_data_combined))
sft <- pickSoftThreshold(data, powerVector = powers)
str(sft)

# Plot the results:
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
####construct the co-expression network,the power= str(sft)
net <- blockwiseModules(data, power = 7, maxBlockSize = 5000,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "FPKM-TOM",
                        verbose = 3)
moduleLabelsAutomatic <- net$colors
moduleColorsAutomatic <- labels2colors(moduleLabelsAutomatic)

table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, data, net, file="net_all.RData")

#save the massage of the genes in each module
table<-net$colors

combat_data_combined <- as.data.frame(combat_data_combined)
net_information<-data.frame(gene=rownames(combat_data_combined),module=table,color=labels2colors(net$colors))
dim(net_information)
write.csv(net_information,file="net_information.csv",quote = F, row.names = T)
save(net_information, file = "net_information.RData", quote = F, row.names = T)





###find the relationship between modules 
net
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
str(net)
head(net)

###plot 
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)


##find the correlation between module and tumor
match_trait <- match(rownames(data), data_pheno_merge$sample)

Traits <- data_pheno_merge[match_trait, -c(1,2,4,6)]
modTraitCor <- cor(MEs, Traits, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, 282)
textMatrix <- paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")",sep = "")
dim(textMatrix) <- dim(modTraitCor)
save(modTraitCor, file = "modTraitCor.RData", quote = F, row.names = T)
save(modTraitP, file = "modTraitP.RData", quote = F, row.names = T)
###figure 4A
par(mar = c(5, 4, 4, 4))
labeledHeatmap(
  Matrix = modTraitCor, xLabels = names(Traits), yLabels = names(MEs),
  ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50),
  textMatrix = textMatrix, setStdMargins = FALSE,
  cex.text = 0.7, zlim = c(-1, 1), main = paste("Module-trait relationships"))

###find the membership between genes and modules
nSamples = nrow(data)
modNames = substring(names(MEs),3)
geneModuleMembership =as.data.frame(cor(data,MEs,use='p'))
write.csv(geneModuleMembership,"geneModuleMembership.csv")
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership)= paste('MM', modNames, sep='')
names(MMPvalue)= paste('p.MM', modNames, sep='')


###GO analysis in each module 
library(org.Hs.eg.db)
library(clusterProfiler)

module0 <- subset(net_information, module == "0" )
module1 <- subset(net_information, module == "1" )
module2 <- subset(net_information, module == "2" )
module3 <- subset(net_information, module == "3" )
module4 <- subset(net_information, module == "4" )
module5 <- subset(net_information, module == "5" )
module6 <- subset(net_information, module == "6" )
module7 <- subset(net_information, module == "7" )
module8 <- subset(net_information, module == "8" )
module9 <- subset(net_information, module == "9" )
module10 <- subset(net_information, module == "10" )
module11 <- subset(net_information, module == "11" )
module12 <- subset(net_information, module == "12" )
module13 <- subset(net_information, module == "13" )


###GO analysis
GO_module5<-enrichGO(gene=rownames(module5),OrgDb = "org.Hs.eg.db",ont="ALL",
                            pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")

write.csv(GO_module12, file = "GO_module12.csv", quote = F, row.names = T)

##dotplot and barplot
barplot(GO_module13, showCategory=30)
dotplot(GO_module13, showCategory=30)








