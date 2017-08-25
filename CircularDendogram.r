.libPaths("/home/vimal/R_Library/")
library(GMD)
library(gplots)
library(ape)
library(cluster) 


FPKM <-read.table("Zscore_FPKM_RNAseq.txt", header=T,row.names=1)
FPKM_mat <-data.matrix(FPKM)
dist.obj<- dist(t(FPKM_mat))
hclust.obj<- hclust(dist.obj)


css.obj <- css.hclust(dist.obj,hclust.obj,k=50)
plot(css.obj[,1],css.obj[,2],type="b")
elbow.obj <- elbow.batch(css.obj)
print(elbow.obj)


heatmap.3(FPKM_mat,dendrogram = "both", kr=46)

mypal = rainbow(31)
clus31 = cutree(hclust.obj, 31)
par(mar =c(0.5, 0.4, 0.4, 0.2) )
plot(as.phylo(hclust.obj),type="fan",cex=0.4,tip.color = mypal[clus31])
