library(GGMselect)
library(network)


###### For Simulation Only
p=30
n=30
eta =0.11
Gr <-simulateGraph(p,eta)
X <-rmvnorm(n,mean=rep(0,p), sigma = Gr$C)
GRest <-selectFast(X, family = "EW")
g <- network(GRest$G)
network(GRest$G)
GRest <-selectFast(X, family = "EW")
g <- network(GRest$EM$G)
a <-plot(gV, usearrow=FALSE)
plot(g,coord=a, usearrow= FALSE)
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(1,2))
a <-plot(gV, usearrow=FALSE)
plot(g,coord=a, usearrow= FALSE)


### Real data
TF_Exp_FPKM <-read.table("/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/Zscore_TF.With.Annotation.txt",row.names=1,header=T)
