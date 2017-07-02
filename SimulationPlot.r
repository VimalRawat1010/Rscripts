library(LSD)
library(vioplot)
library(LSD)
setwd("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Simulations/DNemulator/BSMAP_Analysis/51/")
M51 <-read.table("MethylationLevels_F", header=T, row.names=1)
M75 <-read.table("../75/MethylationLevels_F", header=T, row.names=1)
M101 <-read.table("../101/MethylationLevels_F", header=T, row.names=1)
M <-read.table("../126/MethylationLevels_F", header=T, row.names=1)

png("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Simulations/DNemulator/BSMAP_Analysis/BSMAP_M51_M75_101_126_Error_Methylation_vs_Fold.png",width = 1000, height = 600)
par(mfrow=c(1,4))
vioplot(M51[,2],M51[,3],M51[,4],M51[,5],M51[,6],M51[,7],M51[,8],M51[,9],M51[,10],M51[,11],col="grey")
vioplot(M75[,2],M75[,3],M75[,4],M75[,5],M75[,6],M75[,7],M75[,8],M75[,9],M75[,10],M75[,11],col="red")
vioplot(M101[,2],M101[,3],M101[,4],M101[,5],M101[,6],M101[,7],M101[,8],M101[,9],M101[,10],M101[,11],col="blue")
vioplot(M[,2],M[,3],M[,4],M[,5],M[,6],M[,7],M[,8],M[,9],M[,10],M[,11],col="green")
dev.off()

png("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Simulations/DNemulator/BSMAP_Analysis/BSMAP_M51_126_Error_Methylation_vs_Fold_BOXPLOT.png",width = 1000, height = 600)
par(mfrow=c(1,4))
boxplot(M51[,2],M51[,3],M51[,4],M51[,5],M51[,6],M51[,7],M51[,8],M51[,9],M51[,10],M51[,11],outline=T,col=rainbow(10))
boxplot(M75[,2],M75[,3],M75[,4],M75[,5],M75[,6],M75[,7],M75[,8],M75[,9],M75[,10],M75[,11],outline=T,col=rainbow(10))
boxplot(M101[,2],M101[,3],M101[,4],M101[,5],M101[,6],M101[,7],M101[,8],M101[,9],M101[,10],M101[,11],outline=T,col=rainbow(10))
boxplot(M[,2],M[,3],M[,4],M[,5],M[,6],M[,7],M[,8],M[,9],M[,10],M[,11],outline=T,col=rainbow(10))
dev.off() 


png("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Simulations/DNemulator/BSMAP_Analysis/BSMAP_M51_126_Error_Methylation_vs_Fold_BOXPLOT_2.png",width = 1000, height = 600)
par(mfrow=c(1,4))
boxplot(M51[,2],M51[,3],M51[,4],M51[,5],M51[,6],M51[,7],M51[,8],M51[,9],M51[,10],M51[,11],outline=F,col=rainbow(10))
boxplot(M75[,2],M75[,3],M75[,4],M75[,5],M75[,6],M75[,7],M75[,8],M75[,9],M75[,10],M75[,11],outline=F,col=rainbow(10))
boxplot(M101[,2],M101[,3],M101[,4],M101[,5],M101[,6],M101[,7],M101[,8],M101[,9],M101[,10],M101[,11],outline=F,col=rainbow(10))
boxplot(M[,2],M[,3],M[,4],M[,5],M[,6],M[,7],M[,8],M[,9],M[,10],M[,11],outline=F,col=rainbow(10))
dev.off() 



png("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Simulations/DNemulator/BSMAP_Analysis/BSMAP_M51_M75_Error_Methylation.png",width = 1000, height = 600)
par(mfrow=c(2,2))
heatscatter(M51[,1],M51[,2],cexplot = 1, ylim=c(0,100))
heatscatter(M51[,1],M51[,11],cexplot = 1, ylim=c(0,100))

heatscatter(M75[,1],M75[,2],cexplot = 1, ylim=c(0,100))
heatscatter(M75[,1],M75[,11],cexplot = 1, ylim=c(0,100))
dev.off()

png("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Simulations/DNemulator/BSMAP_Analysis/BSMAP_M101_M126_Error_Methylation.png",width = 1000, height = 600)
par(mfrow=c(2,2))
heatscatter(M101[,1],M101[,2],cexplot = 1, ylim=c(0,100))
heatscatter(M101[,1],M101[,11],cexplot = 1, ylim=c(0,100))

heatscatter(M[,1],M[,2],cexplot = 1, ylim=c(0,100))
heatscatter(M[,1],M[,11],cexplot = 1, ylim=c(0,100))
dev.off()
