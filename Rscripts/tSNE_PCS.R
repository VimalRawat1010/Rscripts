#### PCA & tSNE

library(caret)  
library(Rtsne)
setwd("/home/vimal/Software/Github/Latest/Rscripts/")
######################################################################
## The WHOLE post is in: https://github.com/pablo14/post_cluster_tsne
######################################################################

download.file("https://github.com/pablo14/post_cluster_tsne/blob/master/data_1.txt", "data_1.txt")
## Download data from: https://github.com/pablo14/post_cluster_tsne/blob/master/data_1.txt (url path inside the gitrepo.)
data_tsne=read.delim("data_1.txt", header = T, stringsAsFactors = F, sep = "\t")

## Rtsne function may take some minutes to complete...
set.seed(9)  
tsne_model_1 = Rtsne(as.matrix(data_tsne), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)  