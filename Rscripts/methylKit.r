
library(methylKit)
setwd("/home/ubuntu/data/5TB/Aphid/MD_BEDGRAPH/CG/")


file.list=list( paste("A-01_LE"),
                paste("A-03_LE"),
                paste("A-11_LE"),
                paste("A-12_LE"),
                paste("A-14_LE"),
                paste("A-16_LE"),
                paste("A-20_LE"),
                paste("A-22_LE"),
                paste("A-23_LE"),
                paste("B-29_LE"),
                paste("B-33_LE"),
                paste("B-33_LE"),
                paste("B-36_LE"),
                paste("A-02_BB"),
                paste("A-04_BB"),
                paste("A-05_BB"),
                paste("A-06_BB"),
                paste("A-09_BB"),
                paste("A-13_BB"),
                paste("A-17_BB"),
                paste("A-18_BB"),
                paste("A-19_BB"),
                paste("A-24_BB"),
                paste("B-25_BB"),
                paste("B-26_BB"),
                paste("B-27_BB"),
                paste("B-28_BB"),
                paste("B-31_BB"),
                paste("B-32_BB"))


myobj.LE.BB=methRead(file.list,
                      sample.id=list(
                        "A-01_LE",
                        "A-03_LE",
                        "A-11_LE",
                        "A-12_LE",
                        "A-14_LE",
                        "A-16_LE",
                        "A-20_LE",
                        "A-22_LE",
                        "A-23_LE",
                        "B-29_LE",
                        "B-33_LE",
                        "B-33_LE",
                        "B-36_LE",
                        "A-02_BB",
                        "A-04_BB",
                        "A-05_BB",
                        "A-06_BB",
                        "A-09_BB",
                        "A-13_BB",
                        "A-17_BB",
                        "A-18_BB",
                        "A-19_BB",
                        "A-24_BB",
                        "B-25_BB",
                        "B-26_BB",
                        "B-27_BB",
                        "B-28_BB",
                        "B-31_BB",
                        "B-32_BB"
                        ),
                      assembly="tair10",
                      treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                      context="CpG"
)






#getMethylationStats(myobj[[1]],plot=FALSE,both.strands=FALSE)
#getMethylationStats(myobj.ANS.LE[[1]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.LE.BB[[1]],plot=FALSE,both.strands=FALSE)


### Plot 1
pdf("LE.MethPlot.pdf")
op <- par(mfrow=c(4,3), mar=c(2,2,3,1))

getMethylationStats(myobj.LE.BB[[1]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[2]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[3]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[4]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[5]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[6]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[7]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[9]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[10]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[11]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[12]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.LE.BB[[13]],plot=TRUE,ylim=c(0,100000),both.strands=FALSE,labels = FALSE,breaks=10)

dev.off()
### Plot 2
pdf("BB.MethPlot.pdf")
op <- par(mfrow=c(4,3), mar=c(2,2,3,1))

getCoverageStats(myobj.LE.BB[[14]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[15]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[16]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[17]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[18]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[19]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[20]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[21]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[22]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[23]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)
getCoverageStats(myobj.LE.BB[[24]],plot=TRUE,both.strands=FALSE,ylim=c(0,200000),labels = FALSE,breaks=10)

dev.off()
#meth=unite(myobj, destrand=FALSE)
#meth.ANS.LE=unite(myobj.ANS.LE, destrand=FALSE)
meth.LE.BB=unite(myobj.LE.BB, destrand=FALSE)
meth.LE.BB.min=unite(myobj.LE.BB,min.per.group=5L)

#meth.LE.BB=unite(myobj.LE.BB, destrand=FALSE)
#meth.ANS.BB.min=unite(myobj.LE.BB, min.per.group=2L)

pdf("LE.BB.COR.pdf")
######getCorrelation(meth,plot=TRUE)
getCorrelation(meth.LE.BB,plot=TRUE)
#getCorrelation(meth.ANS.BB,plot=FALSE)
dev.off()



#clusterSamples(meth.ANS.LE, dist="correlation", method="ward.D", plot=TRUE)
pdf("LE.BB.Cluster.pdf")
op <- par(mfrow=c(1,1), mar=c(2,2,3,1))
clusterSamples(meth.LE.BB.min, dist="correlation", method="ward", plot=TRUE)
dev.off()

#hc = clusterSamples(meth.ANS.LE, dist="correlation", plot=FALSE)
pdf("Hclust.pdf")
hc = clusterSamples(meth.LE.BB, method="ward", dist="correlation", plot=TRUE )
dev.off
#hc = clusterSamples(meth.ANS.LE.min, method="ward", dist="correlation", plot=TRUE )
#hc = clusterSamples(meth.ANS.BB, method="ward", dist="correlation", plot=TRUE )
#hc = clusterSamples(meth.ANS.BB.min, method="ward", dist="correlation", plot=TRUE )


#PCASamples(meth.ANS.LE, screeplot=TRUE)
#PCASamples(meth.ANS.LE)
pdf("PCA.LE.BB.pdf")
op <- par(mfrow=c(1,2), mar=c(2,2,3,1))
PCASamples(meth.LE.BB)
PCASamples(meth.LE.BB.min)
dev.off()
#PCASamples(meth.ANS.BB)
#PCASamples(meth.ANS.BB.min)


### Tilling Window
tiles.LE.BB = tileMethylCounts(meth.LE.BB,win.size=1000,step.size = 1000)

####
myDiff.LE=calculateDiffMeth(meth.LE.BB.min)
#myDiff.LE.min=calculateDiffMeth(meth.ANS.LE.min)
#myDiff.BB=calculateDiffMeth(meth.ANS.BB)
#myDiff.BB.min=calculateDiffMeth(meth.ANS.BB.min)

# get hyper methylated bases
sink("Hyper.LE.BB.txt")
myDiff25p.hyper.LE=getMethylDiff(myDiff.LE,difference=0.25,qvalue=1,type="hyper")
sink()
#myDiff25p.hyper.BB=getMethylDiff(myDiff.BB,difference=25,qvalue=0.01,type="hyper")
# get hypo methylated bases
sink("Hypo.LE.BB.txt")
myDiff25p.hypo.LE=getMethylDiff(myDiff.LE,difference=0.25,pvalue=0.01,type="hypo")
sink()
#myDiff25p.hypo.BB=getMethylDiff(myDiff.BB,difference=25,qvalue=0.01,type="hypo")











######################
require(bsseqData)
library(bsseq)
library(DSS)
library(devtools)
library(parallel)
require(graphics)

A07_Ans <- read.table("A-07_Ancestral", header=TRUE)
A08_Ans <- read.table("A-08_Ancestral", header=TRUE)
A10_Ans <- read.table("A-10_Ancestral", header=TRUE)
A15_Ans <- read.table("A-15_Ancestral", header=TRUE)
A21_Ans <- read.table("A-21_Ancestral", header=TRUE)
B30_Ans <- read.table("B-30_Ancestral", header=TRUE)
B34_Ans <- read.table("B-34_Ancestral", header=TRUE)

A01_LE <- read.table("A-01_LE", header=TRUE)
A03_LE <- read.table("A-03_LE", header=TRUE)
A11_LE <- read.table("A-11_LE", header=TRUE)
A12_LE <- read.table("A-12_LE", header=TRUE)
A14_LE <- read.table("A-14_LE", header=TRUE)
A16_LE <- read.table("A-16_LE", header=TRUE)
A20_LE <- read.table("A-20_LE", header=TRUE)
A22_LE <- read.table("A-22_LE", header=TRUE)
A23_LE <- read.table("A-23_LE", header=TRUE)
B29_LE <- read.table("B-29_LE", header=TRUE)
B33_LE <- read.table("B-33_LE", header=TRUE)
B35_LE <- read.table("B-35_LE", header=TRUE)
B36_LE <- read.table("B-36_LE", header=TRUE)

A02_BB <- read.table("A-02_BB", header=TRUE)
A04_BB <- read.table("A-04_BB", header=TRUE)
A05_BB <- read.table("A-05_BB", header=TRUE)
A06_BB <- read.table("A-06_BB", header=TRUE)
A09_BB <- read.table("A-09_BB", header=TRUE)
A13_BB <- read.table("A-13_BB", header=TRUE)
A17_BB <- read.table("A-17_BB", header=TRUE)
A18_BB <- read.table("A-18_BB", header=TRUE)
A19_BB <- read.table("A-19_BB", header=TRUE)
B25_BB <- read.table("B-25_BB", header=TRUE)
B27_BB <- read.table("B-27_BB", header=TRUE)
B28_BB <- read.table("B-28_BB", header=TRUE)
B31_BB <- read.table("B-31_BB", header=TRUE)
B32_BB <- read.table("B-32_BB", header=TRUE)














