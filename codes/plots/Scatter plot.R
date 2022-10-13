
##Fig2D
subtype2 <- allDiff_group2[rownames(allDiff_group1),]
subtype3 <- allDiff_group3[rownames(allDiff_group1),]
dat= data.frame(sub1=allDiff_group1$`subtyp1(logFC)`, sub2=subtype2$`subtyp2(logFC)`, sub3=subtype3$`subtyp3(logFC)`)
dat <- as.data.frame(t(dat))
dat$SUBTYPE <- rownames(dat)
dat$SUBTYPE <- as.factor(dat$SUBTYPE)
dat2= melt(dat,id=1)

pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  return(list(slope,intercept))
}

fit.subtype1_2 = pcreg(alldata$`subtyp1(logFC)`, alldata$`subtyp2(logFC)`)
fit.subtype1_3 = pcreg(alldata$`subtyp1(logFC)`, alldata$`subtyp3(logFC)`)

dat2$variable = as.character(dat2$variable)
dat2$variable = gsub("subtyp2(logFC)", paste("subtyp2(logFC), slope=", signif(fit.subtype1_2[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("subtyp3(logFC)", paste("subtyp3(logFC), slope=", signif(fit.subtype1_3[[1]],2), sep=""), dat2$variable)


library(ggplot2)

subtypes_DEGs=ggplot(dat2,aes(x=sub1,y=value,color=variable)) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"))+ 
  geom_point(alpha=.7, size = 0.5) + 
  geom_abline(slope=1, lty=2) + xlim(-4,4) + ylim(-4,4) + 
  geom_abline(slope=fit.subtype1_2[[1]], intercept = fit.subtype1_2[[2]], color="#E69F00") + 
  geom_abline(slope=fit.subtype1_3[[1]], intercept = fit.subtype1_3[[2]], color="#56B4E9") +
  xlab("subtype1 (log2FC)") + ylab("subtypes (log2FC)") +
  coord_fixed(ratio=1)+
  theme(panel.grid=element_blank(), 
        panel.background=element_rect(color="black", fill="transparent"))##theme




















