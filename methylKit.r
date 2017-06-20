
library(methylKit)

file.list=list( system.file("extdata", 
                            "ANS1", package = "methylKit"),
                system.file("extdata",
                            "ANS2", package = "methylKit"),
                system.file("extdata", 
                            "ANS3", package = "methylKit"),
                system.file("extdata", 
                            "ANS4", package = "methylKit"),
                system.file("extdata", 
                            "ANS5", package = "methylKit"),
                system.file("extdata", 
                            "ANS6", package = "methylKit"),
                system.file("extdata", 
                            "ANS7", package = "methylKit"),
                system.file("extdata", 
                            "LE1", package = "methylKit"),
                system.file("extdata", 
                            "LE2", package = "methylKit"),
                system.file("extdata", 
                            "LE3", package = "methylKit"),
                system.file("extdata", 
                            "LE4", package = "methylKit"),
                system.file("extdata", 
                            "LE5", package = "methylKit"),
                system.file("extdata", 
                            "LE6", package = "methylKit"),
                system.file("extdata", 
                            "LE7", package = "methylKit"),
                system.file("extdata", 
                            "LE8", package = "methylKit"),
                system.file("extdata", 
                            "LE9", package = "methylKit"),
                system.file("extdata", 
                            "LE10", package = "methylKit"),
                system.file("extdata", 
                            "LE11", package = "methylKit"))

myobj.ANS.LE=methRead(file.list,
                  sample.id=list("ANS1","ANS2","ANS3","ANS4","ANS5","ANS6","ANS7","LE1","LE2","LE3","LE4","LE6","LE7","LE8","LE9","LE10","LE11"),
                  assembly="tair10",
                  treatment=c(1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
                  context="CpG"
)

file.list=list( system.file("extdata", 
                            "ANS1", package = "methylKit"),
                system.file("extdata",
                            "ANS2", package = "methylKit"),
                system.file("extdata", 
                            "ANS3", package = "methylKit"),
                system.file("extdata", 
                            "ANS4", package = "methylKit"),
                system.file("extdata", 
                            "ANS5", package = "methylKit"),
                system.file("extdata", 
                            "ANS6", package = "methylKit"),
                system.file("extdata", 
                            "ANS7", package = "methylKit"),
                system.file("extdata", 
                            "BB1", package = "methylKit"),
                system.file("extdata", 
                            "BB2", package = "methylKit"),
                system.file("extdata", 
                            "BB3", package = "methylKit"),
                system.file("extdata", 
                            "BB4", package = "methylKit"),
                system.file("extdata", 
                            "BB5", package = "methylKit"),
                system.file("extdata", 
                            "BB6", package = "methylKit"),
                system.file("extdata", 
                            "BB7", package = "methylKit"),
                system.file("extdata", 
                            "BB9", package = "methylKit"),
                system.file("extdata", 
                            "BB10", package = "methylKit"),
                system.file("extdata", 
                            "BB11", package = "methylKit"),
                system.file("extdata", 
                            "BB13", package = "methylKit"),
                system.file("extdata", 
                            "BB14", package = "methylKit"),
                system.file("extdata", 
                            "BB15", package = "methylKit"),
                system.file("extdata", 
                            "BB16", package = "methylKit"))



myobj.ANS.BB=methRead(file.list,
                      sample.id=list("ANS1","ANS2","ANS3","ANS4","ANS5","ANS6","ANS7","BB1","BB2","BB3","BB4","BB5","BB6","BB7","BB9","BB10","BB11","BB13","BB14","BB15","BB16"),
                      assembly="tair10",
                      treatment=c(1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                      context="CpG"
)






op <- par(mfrow=c(3,3), mar=c(2,2,3,1))
#getMethylationStats(myobj[[1]],plot=FALSE,both.strands=FALSE)
#getMethylationStats(myobj.ANS.LE[[1]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.ANS.LE.BB[[1]],plot=FALSE,both.strands=FALSE)


### Plot 1
getMethylationStats(myobj.ANS.LE.BB[[1]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,chunk.size = 1e+06)
getMethylationStats(myobj.ANS.LE.BB[[2]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,chunk.size = 1e+06,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[3]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,chunk.size = 1e+06)
getMethylationStats(myobj.ANS.LE.BB[[4]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[5]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[6]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[7]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[8]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[9]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[10]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[11]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[12]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[13]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[14]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[15]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)
getMethylationStats(myobj.ANS.LE.BB[[16]],plot=TRUE,ylim=c(0,200000),both.strands=FALSE,labels = FALSE,breaks=10)

### Plot 2
getCoverageStats(myobj.ANS.LE.BB[[1]],plot=TRUE,both.strands=FALSE,ylim=c(0,300000),labels = FALSE,breaks=10)
getCoverageStats(myobj.ANS.LE.BB[[2]],plot=TRUE,both.strands=FALSE,ylim=c(0,300000),labels = FALSE,breaks=10)
getCoverageStats(myobj.ANS.LE.BB[[3]],plot=TRUE,both.strands=FALSE,ylim=c(0,300000),labels = FALSE,breaks=10)
getCoverageStats(myobj.ANS.LE.BB[[4]],plot=TRUE,both.strands=FALSE,ylim=c(0,300000),labels = FALSE,breaks=10)
getCoverageStats(myobj.ANS.LE.BB[[5]],plot=TRUE,both.strands=FALSE,ylim=c(0,300000),labels = FALSE,breaks=10)
getCoverageStats(myobj.ANS.LE.BB[[6]],plot=TRUE,both.strands=FALSE,ylim=c(0,300000),labels = FALSE,breaks=10)
getCoverageStats(myobj.ANS.LE.BB[[7]],plot=TRUE,both.strands=FALSE,ylim=c(0,300000),labels = FALSE,breaks=10)
getCoverageStats(myobj.ANS.LE.BB[[8]],plot=TRUE,both.strands=FALSE,ylim=c(0,300000),labels = FALSE,breaks=10)

#meth=unite(myobj, destrand=FALSE)
#meth.ANS.LE=unite(myobj.ANS.LE, destrand=FALSE)
meth.ANS.LE=unite(myobj.ANS.LE, destrand=FALSE)
meth.ANS.LE.min=unite(myobj.ANS.LE,min.per.group=2L)

meth.ANS.BB=unite(myobj.ANS.BB, destrand=FALSE)
meth.ANS.BB.min=unite(myobj.ANS.BB, min.per.group=2L)




######getCorrelation(meth,plot=TRUE)
#getCorrelation(meth.ANS.LE,plot=FALSE)
#getCorrelation(meth.ANS.BB,plot=FALSE)



op <- par(mfrow=c(1,1), mar=c(2,2,3,1))
#clusterSamples(meth.ANS.LE, dist="correlation", method="ward.D", plot=TRUE)
clusterSamples(meth.min.ANS.LE.BB, dist="correlation", method="ward", plot=TRUE)
#hc = clusterSamples(meth.ANS.LE, dist="correlation", plot=FALSE)
hc = clusterSamples(meth.ANS.LE, method="ward", dist="correlation", plot=TRUE )
hc = clusterSamples(meth.ANS.LE.min, method="ward", dist="correlation", plot=TRUE )

hc = clusterSamples(meth.ANS.BB, method="ward", dist="correlation", plot=TRUE )
hc = clusterSamples(meth.ANS.BB.min, method="ward", dist="correlation", plot=TRUE )


#PCASamples(meth.ANS.LE, screeplot=TRUE)
#PCASamples(meth.ANS.LE)
PCASamples(meth.ANS.LE)
PCASamples(meth.ANS.LE.min)

PCASamples(meth.ANS.BB)
PCASamples(meth.ANS.BB.min)


### Tilling Window
tiles.LE = tileMethylCounts(meth.ANS.LE,win.size=1000,step.size = 1000)


####
myDiff.LE=calculateDiffMeth(meth.ANS.LE)
myDiff.LE.min=calculateDiffMeth(meth.ANS.LE.min)

myDiff.BB=calculateDiffMeth(meth.ANS.BB)
myDiff.BB.min=calculateDiffMeth(meth.ANS.BB.min)

# get hyper methylated bases
myDiff25p.hyper.LE=getMethylDiff(myDiff.LE,difference=0.5,qvalue=1,type="hyper")
myDiff25p.hyper.BB=getMethylDiff(myDiff.BB,difference=25,qvalue=0.01,type="hyper")

#
# get hypo methylated bases
myDiff25p.hypo.LE=getMethylDiff(myDiff.LE,difference=0.25,pvalue=0.01,type="hypo")
myDiff25p.hypo.BB=getMethylDiff(myDiff.BB,difference=25,qvalue=0.01,type="hypo")

#
#
# get all differentially methylated bases
myDiff25p.LE=getMethylDiff(meth.ANS.LE,difference=25,qvalue=0.01)
myDiff25p.BB=getMethylDiff(meth.ANS.BB,difference=25,qvalue=0.01)




# read-in transcript locations to be used in annotation
# IMPORTANT: annotation files that come with the package (version >=0.5) are a subset of full annotation
# files. Download appropriate annotation files from UCSC (or other sources) in BED format
biocLite("genomation")
library(genomation) # install from BioC



