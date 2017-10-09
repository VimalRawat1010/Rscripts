.libPaths("/home/vimal/R_Library/")
library(minet)
library(Rgraphviz)
library(reshape)
library(reshape2)
require(data.table) 

# build.mim 
#takes the dataset as input and computes the mutual information beetween all pair of
#variables according to the mutual inforamtion estimator estimator

TF_Exp_FPKM <-read.table("/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/FPKM_AllGenes_SomeSamples.txt",row.names=1,header=T)
TF_Exp_FPKM_M <-t(data.matrix(TF_Exp_FPKM))

mim <- build.mim(TF_Exp_FPKM_M, estimator = "spearman", disc = "none", nbins = sqrt(NROW(TF_Exp_FPKM_M)))
mim_melt <- setDT(melt(mim))

mim_melt <- mim_melt[mim_melt$value >6]

#newdata <- mim_melt[$V1myvars]
# This function takes the mutual information matrix as input in order to return the infered network
# according to the Aracne algorithm.   This algorithm applies the data processing inequality to all
# triplets of nodes in order to remove the least significant edge in each triplet

network <- aracne( mim, eps=0 )


#mrnetb
#takes the mutual information matrix as input in order to infer the network using the maximum relevance/minimum redundancy
#criterion combined with a backward elimination and a sequential replacement - see references. This method is a variant of mrnet

network <- mrnetb(mim)

plot( as( network ,"graphNEL") )