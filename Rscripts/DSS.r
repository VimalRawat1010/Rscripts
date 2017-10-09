
#source("http://bioconductor.org/biocLite.R")
#biocLite("bsseq")
#biocLite("DSS")
#biocLite("bsseqData")
require(bsseqData)
library(bsseq)
library(DSS)
library(devtools)
library(parallel)
require(graphics)

setwd("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DSS_Analysis/")



################################################################################
################### Loading Data
################################################################################



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
#A16_LE <- read.table("A-16_LE_13_710", header=TRUE)
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
#A18_BB <- read.table("A-18_BB_16_211", header=TRUE)
A19_BB <- read.table("A-19_BB", header=TRUE)
B25_BB <- read.table("B-25_BB", header=TRUE)
B27_BB <- read.table("B-27_BB", header=TRUE)
B28_BB <- read.table("B-28_BB", header=TRUE)
B31_BB <- read.table("B-31_BB", header=TRUE)
B32_BB <- read.table("B-32_BB", header=TRUE)


################################################################################
################### ANS vs LE
################################################################################

## make BSseq objects
BSobj_Ans_LE <- makeBSseqData( list(A07_Ans, A08_Ans, A10_Ans, A15_Ans, A21_Ans, B30_Ans, B34_Ans, A01_LE, A03_LE, 
                             A11_LE, A12_LE, A14_LE, A20_LE, A22_LE, A23_LE, B29_LE, B33_LE, B35_LE, B36_LE), 
                        c("C1","C2", "C3", "C4","C5", "C6", "C7","N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                          "N10" ,"N11", "N12") )


dmlTest_Ans_LE <- DMLtest(BSobj_Ans_LE, group1=c("C1","C2", "C3", "C4","C5", "C6", "C7"), group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                                                                     "N10" ,"N11", "N12"))

dmlTest_Ans_LE.smoothed <- DMLtest(BSobj_Ans_LE, group1=c("C1","C2", "C3", "C4","C5", "C6", "C7"), group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                                                                              "N10" ,"N11", "N12"),smoothing=T)

## call DML
dmls_Ans_LE <- callDML(dmlTest_Ans_LE)
dmls_Ans_LE.smoothed <- callDML(dmlTest_Ans_LE.smoothed)
write.table(dmls_Ans_LE,'Ans_LE_DML.dss.txt',sep='\t',row.names=F)
write.table(dmls_Ans_LE.smoothed,'Ans_LE_DML.smooth.dss.txt',sep='\t',row.names=F)

## call DMR
dmrs_Ans_LE =callDMR(dmlTest_Ans_LE,p.threshold=0.01)
dmrs_Ans_LE.smoothed=callDMR(dmlTest_Ans_LE.smoothed,p.threshold=0.01)
write.table(dmrs_Ans_LE,'Ans_LE_DMR.dss.txt',sep='\t',row.names=F)
write.table(dmrs_Ans_LE.smoothed,'Ans_LE_DMR.smooth.dss.txt',sep='\t',row.names=F)


################################################################################
################### ANS vs BB
################################################################################

BSobj_Ans_BB <- makeBSseqData( list(A07_Ans, A08_Ans, A10_Ans, A15_Ans, A21_Ans, B30_Ans, B34_Ans, A02_BB, A04_BB, 
                                    A05_BB, A06_BB, A09_BB, A13_BB, A17_BB, A19_BB, B25_BB, B27_BB, B28_BB, B31_BB, B32_BB), 
                               c("C1","C2", "C3", "C4","C5", "C6", "C7","N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                 "N10" ,"N11", "N12" ,"N13") )


dmlTest_Ans_BB <- DMLtest(BSobj_Ans_BB, group1=c("C1","C2", "C3", "C4","C5", "C6", "C7"), group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                                                                            "N10" ,"N11", "N12" ,"N13"))

dmlTest_Ans_BB.smoothed <- DMLtest(BSobj_Ans_BB, group1=c("C1","C2", "C3", "C4","C5", "C6", "C7"), group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                                                                                     "N10" ,"N11", "N12" ,"N13"),smoothing=T)

## call DML
dmls_Ans_BB <- callDML(dmlTest_Ans_BB)
dmls_Ans_BB.smoothed <- callDML(dmlTest_Ans_BB.smoothed)
write.table(dmls_Ans_BB,'Ans_BB_DML.dss.txt',sep='\t',row.names=F)
write.table(dmls_Ans_BB.smoothed,'Ans_BB_DML.smooth.dss.txt',sep='\t',row.names=F)

## call DMR
dmrs_Ans_BB =callDMR(dmlTest_Ans_BB,p.threshold=0.01)
dmrs_Ans_BB.smoothed=callDMR(dmlTest_Ans_BB.smoothed,p.threshold=0.01)
write.table(dmrs_Ans_BB,'Ans_BB_DMR.dss.txt',sep='\t',row.names=F)
write.table(dmrs_Ans_BB.smoothed,'Ans_BB_DMR.smooth.dss.txt',sep='\t',row.names=F)




################################################################################
################### LE vs BB
################################################################################

## make BSseq objects
BSobj_LE_BB <- makeBSseqData(
            list(A01_LE,A03_LE,A11_LE,A12_LE,A14_LE,A20_LE,A22_LE,A23_LE,B29_LE,B33_LE,B35_LE,B36_LE,
                 A02_BB,A04_BB,A05_BB,A06_BB,A09_BB,A13_BB,A17_BB,A19_BB,B25_BB,B27_BB,B28_BB,B31_BB,B32_BB),
            c("C1","C2", "C3", "C4","C5", "C6", "C7", "C8", "C9","C10", "C11", "C12",
              "N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9","N10" ,"N11", "N12" ,"N13"))


dmlTest_LE_BB <- DMLtest(BSobj_LE_BB, group1=c("C1","C2", "C3", "C4","C5", "C6", "C7", "C8", "C9","C10", "C11", "C12"),
                         group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9", "N10" ,"N11", "N12" ,"N13"))


dmlTest_LE_BB.smoothed <- DMLtest(
                          BSobj_LE_BB, 
                          group1=c("C1","C2", "C3", "C4","C5", "C6", "C7", "C8", "C9","C10", "C11", "C12"), 
                          group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9","N10" ,"N11", "N12","N13"),
                          smoothing=T)

## call DML
dmls_LE_BB <- callDML(dmlTest_LE_BB)
saveRDS(dmlTest_LE_BB, "/home/ubuntu/data/5TB/Aphid/BedGraph/INPUT/DSS_Analysis/All_dml_LE_BB.RDS")
#dmls_LE_BB.smoothed <- callDML(dmlTest_LE_BB.smoothed)
write.table(dmls_LE_BB,'LE_BB_DML.dss.txt',sep='\t',row.names=F)
#write.table(dmls_LE_BB.smoothed,'LE_BB_DML.smooth.dss.txt',sep='\t',row.names=F)

## call DMR
dmrs_LE_BB =callDMR(dmlTest_LE_BB,p.threshold=0.01)
#dmrs_LE_BB.smoothed=callDMR(dmlTest_LE_BB.smoothed,p.threshold=0.01)
write.table(dmrs_LE_BB,'LE_BB_DMR.dss.txt',sep='\t',row.names=F)
#write.table(dmrs_LE_BB.smoothed,'LE_BB_DMR.smooth.dss.txt',sep='\t',row.names=F)


