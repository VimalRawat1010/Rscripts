
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

setwd("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DML/")



################################################################################
################### Loading Data
################################################################################

dat1.1 <- read.table("A-07_Ancestral_0_824", header=TRUE)
dat1.2 <- read.table("A-08_Ancestral_0_817", header=TRUE)
dat1.3 <- read.table("A-10_Ancestral_0_829", header=TRUE)
dat2.1 <- read.table("A-01_LE_13_713", header=TRUE)
dat2.2 <- read.table("A-03_LE_13_716", header=TRUE)
dat2.3 <- read.table("A-11_LE_4_743", header=TRUE)

A07_Ans <- read.table("A-07_Ancestral_0_824", header=TRUE)
A08_Ans <- read.table("A-08_Ancestral_0_817", header=TRUE)
A10_Ans <- read.table("A-10_Ancestral_0_829", header=TRUE)
A15_Ans <- read.table("A-15_Ancestral_0_812", header=TRUE)
A21_Ans <- read.table("A-21_Ancestral_0_830", header=TRUE)
B30_Ans <- read.table("B-30_Ancestral_0_811", header=TRUE)
B34_Ans <- read.table("B-34_Ancestral_0_823", header=TRUE)

A01_LE <- read.table("A-01_LE_13_713", header=TRUE)
A03_LE <- read.table("A-03_LE_13_716", header=TRUE)
A11_LE <- read.table("A-11_LE_4_743", header=TRUE)
A12_LE <- read.table("A-12_LE_28_787", header=TRUE)
A14_LE <- read.table("A-14_LE_4_748", header=TRUE)
A16_LE <- read.table("A-16_LE_13_710", header=TRUE)
A20_LE <- read.table("A-20_LE_4_744", header=TRUE)
A22_LE <- read.table("A-22_LE_28_785", header=TRUE)
A23_LE <- read.table("A-23_LE_13_712", header=TRUE)
B29_LE <- read.table("B-29_LE_28_795", header=TRUE)
B33_LE <- read.table("B-33_LE_4_747", header=TRUE)
B35_LE <- read.table("B-35_LE_13_715", header=TRUE)
B36_LE <- read.table("B-36_LE_28_799", header=TRUE)


A02_BB <- read.table("A-02_BB_27_698", header=TRUE)
A04_BB <- read.table("A-04_BB_19_652", header=TRUE)
A05_BB <- read.table("A-05_BB_19_647", header=TRUE)
A06_BB <- read.table("A-06_BB_27_696", header=TRUE)
A09_BB <- read.table("A-09_BB_19_650", header=TRUE)
A13_BB <- read.table("A-13_BB_19_644", header=TRUE)
A17_BB <- read.table("A-17_BB_16_575", header=TRUE)
A18_BB <- read.table("A-18_BB_16_211", header=TRUE)
A19_BB <- read.table("A-19_BB_27_688", header=TRUE)
B25_BB <- read.table("B-25_BB_14_563", header=TRUE)
B27_BB <- read.table("B-27_BB_16_702", header=TRUE)
B28_BB <- read.table("B-28_BB_27_689", header=TRUE)
B31_BB <- read.table("B-31_BB_19_645", header=TRUE)
B32_BB <- read.table("B-32_BB_27_690", header=TRUE)


################################################################################
################### ANS vs LE
################################################################################

## make BSseq objects
BSobj_Ans_LE <- makeBSseqData( list(A07_Ans, A08_Ans, A10_Ans, A15_Ans, A21_Ans, B30_Ans, B34_Ans, A01_LE, A03_LE, 
                             A11_LE, A12_LE, A14_LE, A16_LE, A20_LE, A22_LE, A23_LE, B29_LE, B33_LE, B35_LE, B36_LE), 
                        c("C1","C2", "C3", "C4","C5", "C6", "C7","N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                          "N10" ,"N11", "N12" ,"N13") )


dmlTest_Ans_LE <- DMLtest(BSobj_Ans_LE, group1=c("C1","C2", "C3", "C4","C5", "C6", "C7"), group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                                                                     "N10" ,"N11", "N12" ,"N13"))

dmlTest_Ans_LE.smoothed <- DMLtest(BSobj_Ans_LE, group1=c("C1","C2", "C3", "C4","C5", "C6", "C7"), group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                                                                              "N10" ,"N11", "N12" ,"N13"),smoothing=T)

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
                                    A05_BB, A06_BB, A09_BB, A23_BB, A17_BB, A18_BB, A19_BB, B25_BB, B27_BB, B28_BB, B31_BB, B32_BB), 
                               c("C1","C2", "C3", "C4","C5", "C6", "C7","N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                 "N10" ,"N11", "N12" ,"N13" ,"N14") )


dmlTest_Ans_BB <- DMLtest(BSobj_Ans_BB, group1=c("C1","C2", "C3", "C4","C5", "C6", "C7"), group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                                                                            "N10" ,"N11", "N12" ,"N13"  ,"N14"))

dmlTest_Ans_BB.smoothed <- DMLtest(BSobj_Ans_BB, group1=c("C1","C2", "C3", "C4","C5", "C6", "C7"), group2=c("N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                                                                                     "N10" ,"N11", "N12" ,"N13"  ,"N14"),smoothing=T)

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
################### LE vs LE  && BB vs BB
################################################################################

## make BSseq objects
BSobj <- makeBSseqData( list(dat1.1, dat1.2, dat1.3, dat2.1, dat2.2, dat2.3), c("C1","C2", "C3" ,"N1", "N2", "N3") )
##  DML test
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2", "C3"), group2=c("N1","N2","N3"))
dmlTest.smoothed <- DMLtest(BSobj, group1=c("C1", "C2", "C3"), group2=c("N1","N2","N3"),smoothing=T)
## call DML
dmls <- callDML(dmlTest)
dmls.smoothed <- callDML(dmlTest.smoothed)
write.table(dmls,'Ans_LE_DML.dss.txt',sep='\t',row.names=F)
write.table(dmls.smoothed,'Ans_LE_DML.smooth.dss.txt',sep='\t',row.names=F)

## call DMR
dmrs=callDMR(dmlTest,p.threshold=0.01)
dmrs.smoothed=callDMR(dmlTest.smoothed,p.threshold=0.01)
write.table(dmrs,'Ans_LE_DMR.dss.txt',sep='\t',row.names=F)
write.table(dmrs.smoothed,'Ans_LE_DMR.smooth.dss.txt',sep='\t',row.names=F)

