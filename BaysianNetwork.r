## try http:// if https:// URLs are not supported
###  https://www.bioconductor.org/packages/release/data/experiment/vignettes/DREAM4/inst/doc/DREAM4.pdf

source("https://bioconductor.org/biocLite.R")
biocLite("networkBMA")


TF_Exp_FPKM <-read.table("/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/Zscore_TF.With.Annotation.txt",row.names=1,header=T)
TF_Exp_FPKM_M <-(t(data.matrix(TF_Exp_FPKM)))
TF_Inferred <- networkBMA(TF_Exp_FPKM_M, nTimePoints=nrow(TF_Exp_FPKM_M), self=FALSE)
TF_Inferred <- TF_Inferred[order(TF_Inferred$PostProb, decreasing=TRUE),]

print(head(TF_Inferred, n=10))