### Bar plot upward Downward
par(mfrow=c(1,3))
### CC
data<-matrix(c(46,-54,26,-74,23,-77),ncol=3)
### EC
data<-matrix(c(58,-42,56,-44,62,-38),ncol=3)
### SC
data<-matrix(c(30,-70,27,-73,26,-74),ncol=3)

colnames(data)=c("Seedling","Root","Arial")
data1 <- data2 <- data
data1[data1<0] <- 0
data2[data2>0] <- 0
barplot(data1,ylim=c(-150,150),col="green",width=1, main ="Synergid Cell (44,26,34)")
barplot(data2,add=TRUE,ylim=c(-100,100),col="red",width=1)
legend(2,150,  fill= c("green","red"), c("H3K4Me3", "H3K27Me3"), bty = "n"  );
#

### BG line plotting

### Seedling
segments(0.2,65,1.2,65,lwd=1,lty=2)
segments(0.2,-35,1.2,-35,lwd=1,lty=2)
segments(0.2,65,1.2,65,lwd=1,lty=2)
segments(0.2,-35,1.2,-35,lwd=1,lty=2)

### Root
segments(1.4,68,2.4,68,lwd=1,lty=2)
segments(1.4,-32,2.4,-32,lwd=1,lty=2)
segments(1.4,68,2.4,68,lwd=1,lty=2)
segments(1.4,-32,2.4,-32,lwd=1,lty=2)

### Aerial
segments(2.6,69,3.6,69,lwd=1,lty=2)
segments(2.6,-31,3.6,-31,lwd=1,lty=2)
segments(2.6,69,3.6,69,lwd=1,lty=2)
segments(2.6,-31,3.6,-31,lwd=1,lty=2)
