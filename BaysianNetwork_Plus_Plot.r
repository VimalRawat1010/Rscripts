#install.packages("igraph")
#install.packages("network") 
#install.packages("sna")
#install.packages("ndtv")


## TIME SERIES DATA


install.packages("igraph", dependencies=TRUE, repos='http://cran.rstudio.com/')

library("networkBMA")
library("igraph")
library("network") 
library("sna")
library("ndtv")





TF_Exp_FPKM <-read.table("/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/FPKM_TF_GENE_CELL_TYPES_short.txt",row.names=1,header=T)
TF_Exp_FPKM_M <-(t(data.matrix(TF_Exp_FPKM)))
TF_Inferred <- networkBMA(TF_Exp_FPKM_M, nTimePoints=nrow(TF_Exp_FPKM_M), self=FALSE)
TF_Inferred <- TF_Inferred[order(TF_Inferred$PostProb, decreasing=TRUE),]
