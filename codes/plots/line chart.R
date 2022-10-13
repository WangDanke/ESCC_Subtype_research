######figure2B-C
data <- read.csv("Amplification genes.csv",row.names=1)
data <- read.csv("Deletion or insertion genes.csv",row.names=1)

colnames(data) <- c("subtype1","subtype2","subtype3")
data <- as.data.frame(t(data))
data$SUBTYPE <- rownames(data)
data$SUBTYPE <- as.factor(data$SUBTYPE)
library(reshape2)
data <- melt(data, id="SUBTYPE")
#
ggplot(data,aes(x=SUBTYPE, y=value, group=variable, color=variable))+
  geom_line()+
  geom_point(size=3.5, shape=20)+
  labs(x="SUBTYPE", y="logFC")+
  theme(panel.grid=element_blank(), 
        panel.background=element_rect(color="black", fill="transparent"))