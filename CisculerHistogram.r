library(circular)
set.seed(1)
Omega=2*pi*A[,1]/24
OmegaAA=2*pi*B[,1]/24

A <-read.table("PRR1_AT")
B <-read.table("PRR1_AA")

Omegat=2*pi*trunc(A[,1])/24
OmegatAA=2*pi*trunc(B[,1])/24

#### A.thaliana
H=circular(Omega,type="angle",units="radians",rotation="clock")
Ht=circular(Omegat,type="angle",units="radians",rotation="clock")
#### A.alpina
HAA=circular(OmegaAA,type="angle",units="radians",rotation="clock")
HtAA=circular(OmegatAA,type="angle",units="radians",rotation="clock")





#### A.thaliana
H=circular(Omega,type="angle",units="radians",rotation="clock")
Ht=circular(Omegat,type="angle",units="radians",rotation="clock")
#### A.alpina
HAA=circular(OmegaAA,type="angle",units="radians",rotation="clock")
HtAA=circular(OmegatAA,type="angle",units="radians",rotation="clock")



plot(Ht, stack=FALSE, shrink=1.3, cex=1.03,axes=TRUE,tol=0.8,zero=c(rad(90)),bins=24,ylim=c(0,1))
points(Ht, rotation = "clock", zero =c(rad(90)),col = "1", cex=0.03, stack=TRUE )

rose.diag(Ht-pi/2,bins=24,shrink=0.5,xlim=c(-2,2),ylim=c(-2,2),axes=FALSE,prop=1.5)


#### Athaliana plot
circ.dens = density(Ht+3*pi/2,bw=20)
#### Aalpina plot
circ.densAA = density(HtAA+3*pi/2,bw=20)




#### Plotting Clock
plot(Ht, stack=TRUE, shrink=.35, cex=0, sep=0.0,axes=FALSE,tol=.8,zero=c(0),bins=24,xlim=c(-2,2),ylim=c(-2,2.3), ticks=TRUE, tcl=.075,main="PRR1 Shift")
text(0,0.8,"24", cex=2); text(0,-0.8,"12",cex=2); 
text(0.8,0,"6",cex=2); text(-0.8,0,"18",cex=2)


###### Plotting At Diurnal gene density
lines(circ.dens, col="red", lwd=3)
lines(circ.densAA, col="blue", lwd=3)

