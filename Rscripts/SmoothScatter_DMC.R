library(RColorBrewer)
k <- 100
my.cols <- rev(brewer.pal(k, "RdYlBu"))

setwd("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DSS_Analysis/")
Ans_LE <-read.table("Ans_LE_DML.dss.txt", header = T)
Ans_BB <-read.table("Ans_BB_DML.dss.txt", header = T)
LE_BB <-read.table("LE_BB_DML.dss.txt", header = T)

Ans_LE$mu1 <-as.numeric(as.character(Ans_LE$mu1)) *100
Ans_LE$mu2 <-as.numeric(as.character(Ans_LE$mu2)) *100
Ans_BB$mu1 <-as.numeric(as.character(Ans_BB$mu1)) *100
Ans_BB$mu2 <-as.numeric(as.character(Ans_BB$mu2)) *100
LE_BB$mu1 <-as.numeric(as.character(LE_BB$mu1)) *100
LE_BB$mu2 <-as.numeric(as.character(LE_BB$mu2)) *100

#na.rm = TRUE

ANS_LE_ALL <- readRDS("All_dml_Ans_LE.RDS", refhook = NULL)
ANS_BB_ALL <- readRDS("All_dml_Ans_BB.RDS", refhook = NULL)
LE_BB_ALL <- readRDS("All_dml_LE_BB.RDS", refhook = NULL)

ANS_LE_ALL$mu1 <-as.numeric(as.character(ANS_LE_ALL$mu1)) *100
ANS_LE_ALL$mu2 <-as.numeric(as.character(ANS_LE_ALL$mu2)) *100
ANS_BB_ALL$mu1 <-as.numeric(as.character(ANS_BB_ALL$mu1)) *100
ANS_BB_ALL$mu2 <-as.numeric(as.character(ANS_BB_ALL$mu2)) *100
LE_BB_ALL$mu1 <-as.numeric(as.character(LE_BB_ALL$mu1)) *100
LE_BB_ALL$mu2 <-as.numeric(as.character(LE_BB_ALL$mu2)) *100



par(mfrow=c(2,3))
smoothScatter(ANS_LE_ALL$mu1,ANS_LE_ALL$mu2,nrpoints = Inf,cex=1,col=adjustcolor( "black", alpha.f = 0.1))
smoothScatter(ANS_BB_ALL$mu1,ANS_BB_ALL$mu2,nrpoints = Inf,cex=1,col=adjustcolor( "black", alpha.f = 0.1))
smoothScatter(LE_BB_ALL$mu1,LE_BB_ALL$mu2,nrpoints = Inf,cex=1,col=adjustcolor( "black", alpha.f = 0.1))

smoothScatter(Ans_LE$mu1,Ans_LE$mu2,nrpoints = Inf,cex=5,col=adjustcolor( "black", alpha.f = 0.1))
smoothScatter(Ans_BB$mu1,Ans_BB$mu2,nrpoints = Inf,cex=5,col=adjustcolor( "black", alpha.f = 0.1))
smoothScatter(LE_BB$mu1,LE_BB$mu2,nrpoints = Inf,cex=5,col=adjustcolor( "black", alpha.f = 0.1))
