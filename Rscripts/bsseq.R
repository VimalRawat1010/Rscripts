.libPaths("/home/vimal/R_Library/")
library(bsseqData)
library(bsseq)
library(DSS)


setwd("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DSS_Analysis/methylSig")


fileList = c(system.file("extdata", "AML_1.txt", package = "methylSig"),
             system.file("extdata", "AML_2.txt", package = "methylSig"),
             system.file("extdata", "AML_3.txt", package = "methylSig"),
             system.file("extdata", "AML_4.txt", package = "methylSig"),
             system.file("extdata", "NBM_1.txt", package = "methylSig"),
             system.file("extdata", "NBM_2.txt", package = "methylSig"),
             system.file("extdata", "NBM_3.txt", package = "methylSig"),
             system.file("extdata", "NBM_4.txt", package = "methylSig"))


A07_Ans <- read.table("A-07_Ancestral_0_824.MethylKit", header=TRUE)
A08_Ans <- read.table("A-08_Ancestral_0_817.MethylKit", header=TRUE)


A01_LE <- read.table("A-01_LE_13_713.MethylKit", header=TRUE)
A03_LE <- read.table("A-03_LE_13_716.MethylKit", header=TRUE)


BSobj_Ans_LE <- makeBSseqData( list(A07_Ans, A08_Ans, A10_Ans, A15_Ans, A21_Ans, B30_Ans, B34_Ans, A01_LE, A03_LE, 
                                    A11_LE, A12_LE, A14_LE, A20_LE, A22_LE, A23_LE, B29_LE, B33_LE, B35_LE, B36_LE), 
                               c("C1","C2", "C3", "C4","C5", "C6", "C7","N1", "N2" ,"N3", "N4", "N5" ,"N6", "N7" ,"N8", "N9",
                                 "N10" ,"N11", "N12") )
