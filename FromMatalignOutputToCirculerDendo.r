library(GMD)
setwd("")
All_motifs <-read.table("Combined_pair_list_filtre.out.mat",header=T,row.names=1)
All_motifs_mat <-data.matrix(All_motifs)
dist.obj <- dist(All_motifs_mat)
hclust.obj <-hclust(dist.obj)
css.obj <- css.hclust(dist.obj,hclust.obj,k=50)
elbow.obj <- elbow.batch(css.obj)
print(elbow.obj)
#### Elbow Point is 44
cutree.obj <- cutree(hclust.obj,k=44)


colorCodes <- c(T="red", A ="blue")

## function to set label color
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label")
    code <- substr(label, 2, 2)
    ## use the following line to reset the label to one letter code
    # attr(x, "label") <- code
    attr(x, "nodePar") <- list(lab.col=colorCodes[code])
  }
  return(x)
}

##### OR
#  clus5 = cutree(hclust.obj, 5)
#  mypal = c("red","blue","green" ,"grey" ,"cyan")
#  plot(as.phylo(hclust.obj), type = "fan", tip.color = mypal[clus5], cex=0.3)
#########



library("ape")
plot(as.phylo(hclust.obj), tip.color=colorCodes[substr(rownames(All_motifs), 2, 2)], type="fan",cex=0.4)
