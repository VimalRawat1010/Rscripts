setwd("/projects/dep_coupland/grp_nordstrom/projects/VimalData/PhylogenomicShadowing/PromoterUpstream/DHS_Validation/DHS_correction/")
par(mfrow=c(4,4))
for (i in (15:18))
{
  
  str <- paste("Plot for Motif", i, sep=" ")
  file <- paste("INPUT_FOR_PLOT/", i, sep="")
  data <-read.table(file,header=F,row.names=1)
  a = seq(-100,99,by=1)
  
  smoothingSpline1 = smooth.spline(a, data[,1], spar=0.3)
  smoothingSpline2 = smooth.spline(a, data[,2], spar=0.3)
  plot(smoothingSpline1,col="red",main=str,type="l",ylab="Average DNase Activity",lwd=2)
  lines(smoothingSpline2,col="black",lwd=2)
}