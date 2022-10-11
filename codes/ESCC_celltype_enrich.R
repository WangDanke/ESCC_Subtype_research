dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\combined-array")
setwd(dir)
getwd()

####load the net imformation file
options(stringsAsFactors = TRUE)###如果不运行这个,后面的module_level <- levels(escc$module)结果是NULL
escc <- net_information

escc$module<- as.factor(escc$module)
module_level <- levels(escc$module)
module_number<- length(module_level)
length(rownames(escc))
##load the marker gene list
marker.list<-read.table("PanglaoDB_markers_27_Mar_2020.tsv",header=T,sep="\t")
marker.level<-levels(marker.list$cell.type)

head(marker.list)
length(marker.level)
#
colnames(escc)[1]<- 'official.gene.symbol'
escc_cell_type<-merge(x=escc,y=marker.list,by='official.gene.symbol',all.x=T,all.y=F)
write.csv(escc_cell_type,file = "escc_cell_type.csv")


####make a list of modules
cluster.markers<-list()
for (i in 1:length(levels(escc_cell_type$module))) {
  cluster.markers[[i]]<-subset(escc_cell_type, module==i-1)
}

####Fisher.test
ptable<-data.frame()#pvalue of fisher test
otable<-data.frame()#odds ratio
ctable<-data.frame()#count of marker genes in each module
for (j in 1:length(marker.level)){
  for(i in 1:length(module_level)){
    N <- length(rownames(escc_cell_type))
    m <- length(rownames(subset(marker.list,marker.list$cell.type==marker.level[j])))
    n <- length(rownames(cluster.markers[[i]]))
    k <- length(intersect(cluster.markers[[i]][,1],
                          subset(marker.list,marker.list$cell.type==marker.level[j])[,2]))
    ptable[i,j]<-fisher.test(rbind(c(m-k,N-m-n+k),c(k,n-k)))$p.value
    ctable[i,j]<-k
    otable[i,j]<-fisher.test(rbind(c(m-k,N-m-n+k),c(k,n-k)))$estimate 
  }
}
###
colnames(otable)<-colnames(ptable)<-colnames(ctable)<-marker.level
rownames(otable)<-rownames(ptable)<-rownames(ctable)<-module_level
#save the results data
rownames(ptable)<-paste(0:13,"module",sep="")
rownames(ctable)<-paste(0:13,"module",sep="")
rownames(otable)<-paste(0:13,"module",sep="")


write.csv(ptable,"MGEnrich_ptable.csv")
write.csv(ctable,"MGEnrich_ctable.csv")
write.csv(otable,"MGEnrich_otable.csv")
####find the Pvalue<0.05 cell
A<-data.frame()
for (i in 1:length(rownames(ptable))){
  for (j in 1:length(colnames(ptable))){
    A[i,j]<-as.data.frame(ifelse(ptable[i,j] < 0.05, ptable[i,j], NA))
  }
}


rownames(A)<-rownames(ptable)
rownames(A)<-paste(0:13,"module",sep="")
colnames(A)<-colnames(ptable)



write.csv(A,file='cell_type_enrich_all.csv')


















