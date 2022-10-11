dir <- choose.dir()##dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE53625_RAW")
setwd(dir)
getwd()
library(ggplot2)
library(ggsignif)
####boxplot
e <- subset(data_pheno_merge,sample_type == "tumor")
data <- combat_data_combined[,rownames(e)]

data <- as.data.frame(data)
data <- as.data.frame(t(data))
data$SUBTYPE<-e$SUBTYPE
data$SUBTYPE<-as.factor(data$SUBTYPE)

compaired <- list(c("subtype1", "subtype2"),c("subtype2","subtype3"),c("subtype1","subtype3"))

##Plot a boxplot of a single gene
ggplot(data, aes(x=SUBTYPE, y=MYC, color=SUBTYPE)) + 
  geom_boxplot()+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)+
  theme(panel.grid=element_blank(), 
        panel.background=element_rect(color="black", fill="transparent"))




















