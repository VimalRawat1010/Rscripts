source("http://bioconductor.org/biocLite.R")
biocLite("bsseq")
biocLite("DSS")
biocLite("bsseqData")
require(bsseqData)
library(bsseq)
library(DSS)
library(devtools)
library(parallel)
require(graphics)

setwd("~/Desktop/")

dat1.1 <- read.table("A-07_Ancestral_0_824", header=TRUE)
dat1.2 <- read.table("A-08_Ancestral_0_817", header=TRUE)
dat1.3 <- read.table("A-10_Ancestral_0_829", header=TRUE)
dat2.1 <- read.table("A-01_LE_13_713", header=TRUE)
dat2.2 <- read.table("A-03_LE_13_716", header=TRUE)
dat2.3 <- read.table("A-11_LE_4_743", header=TRUE)

## make BSseq objects
BSobj <- makeBSseqData( list(dat1.1, dat1.2, dat1.3, dat2.1, dat2.2, dat2.3), c("C1","C2", "C3" ,"N1", "N2", "N3") )
##  DML test
dmlTest.smoothed <- DMLtest(BSobj, group1=c("C1", "C2", "C3"), group2=c("N1","N2","N3"),smoothing=T)
## call DML
dmls <- callDML(dmlTest)
dmls.smoothed <- callDML(dmlTest.smoothed)
write.table(dmls.smoothed,'Ans_LE_DML.dss.txt',sep='\t',row.names=F)

## call DMR
dmrs=callDMR(dmlTest,p.threshold=0.01)
dmrs.smoothed=callDMR(dmlTest.smoothed,p.threshold=0.01)
write.table(dmrs.smoothed,'Ans_LE_DMR.dss.txt',sep='\t',row.names=F)







## For whole-genome BS-seq data, perform DML test with smoothing
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2", "C3"), group2=c("N1","N2","N3"), smoothing=TRUE, smoothing.span=500)
dmls <- callDML(dmlTest)
head(dmls)
