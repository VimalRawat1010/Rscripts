#### Enrichment Plot for AT & AA

setwd("/projects/dep_coupland/grp_nordstrom/projects/VimalData/PhylogenomicShadowing/PromoterUpstream/SHIFT_ENRICHMENT/")


##### Function for Changing alpha values
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
} 


myColours = c("red", "blue")
myColoursAlpha <- add.alpha(myColours, alpha=0.4) 



par(mfrow=c(4,4))
Data <-read.table("test",row.names=1)
Data <-data.matrix(Data)
tp <-seq(2,24,by=4)
for (i in seq(1,(nrow(Data)-1),2))
{
  j <- i +1
  label <- rownames(Data)
  str <- paste("Plot",label[i],"and" ,label[j], sep=" ")
 
x <-c(0,23)
y <-c(0,0)
plot(tp,Data[i,],col=myColoursAlpha[1],ylim=c(-2,2), main =str, pch=19,ylab="Z score")
points(tp,Data[j,],col=myColoursAlpha[2],add=T, pch=19)
smoothingSpline1 = smooth.spline(tp, Data[i,], spar=0)
smoothingSpline2 = smooth.spline(tp, Data[j,], spar=0)
lines(smoothingSpline1,col=myColoursAlpha[1])
lines(smoothingSpline2,col=myColoursAlpha[2])
lines(x,y)

AT <- Data[i,]
AA <- Data[j,]
#ccfvalues = ccf(AT,AA,main=str)


}
