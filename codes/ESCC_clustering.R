dir <- choose.dir(default = "", caption = "")
setwd(dir)
getwd()
###Use the tumor samples to clustering
pheno_tumor<- data_pheno_merge[which(data_pheno_merge$sample_type=="tumor"),]
mydata <- t(combat_data_combined[,rownames(pheno_tumor)])
##hierarchical clustering
d <- dist(mydata, method = "euclidean") 
fit2 <- hclust(d, method="ward.D")

plot(fit2)
hcd = as.dendrogram(fit2)
##cut in to three groups
groups <- cutree(fit2, k=3)

###Figure 2A
# load code of A2R function
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
# colored dendrogram
op = par(bg = "#FFFFFF")
A2Rplot(fit2, k = 3, boxes = FALSE, col.up = "gray50", col.down = c("#999999", "#E69F00", "#56B4E9"),
        main = "clustering of 141 ESCC samples")

colnames(groups)[1]<-"clstering_groups"
save(groups,file = "clustering_groups.RData", quote = F, row.names = T)
write.csv(groups,file = "clustering_groups.csv", quote = F, row.names = T)


####Supplementary figure 2A
####hierarchical clustering tree colored by batch
temp_col <- c("yellow","green","purple","blue","red","pink","grey")[as.factor(pheno_tumor$batch)]
temp_col <- temp_col[order.dendrogram(hcd)]
temp_col <- factor(temp_col, unique(temp_col))

library(dendextend)
hcd %>% color_branches(clusters = as.numeric(temp_col), col = levels(temp_col)) %>% 
  set("labels_colors", as.character(temp_col)) %>% 
  plot


#Use another clustering method(K-means) to validate that our Subtypes are not sensitive to the clustering method
ESCC_141_kmeans<-kmeans(mydata,centers=3,nstart = 25)
names(ESCC_141_kmeans)
View(ESCC_141_kmeans)
Kmean_groups <- ESCC_141_kmeans[["cluster"]]
Kmean_groups <- as.data.frame(Kmean_groups)
pheno_tumor$Kmeans_Subtype <- Kmean_groups$Kmean_groups

#Supplementary figure 2B
####Hierarchical clustering tree graph, color labels are divided according to the subtypes clustered by K-means method
temp_col <- c("#56B4E9", "#E69F00","#999999" )[as.factor(pheno_tumor$Kmeans_Subtype)]
temp_col <- temp_col[order.dendrogram(hcd)]
temp_col <- factor(temp_col, unique(temp_col))
hcd %>% color_branches(clusters = as.numeric(temp_col), col = levels(temp_col)) %>% 
  set("labels_colors", as.character(temp_col)) %>% 
  plot

















