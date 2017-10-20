#source("https://bioconductor.org/biocLite.R")
#biocLite("OmicCircos")
library(OmicCircos);

chr.n  <- 7;
chrN     <- paste0("Chr", 1:5);
chrEN     <- c("ChrC", "ChrM");
chr = c(chrN, chrEN)
val     <- c(10, 11,20,15,20,0.05,0.1) ;
val2     <- c(15,15,15,15,14) ;
seg.dat <- data.frame(chr=chr, start=rep(1, chr.n), end=c(34.964571,22.037565,25.499034,20.862711,31.270811,0.154478,0.366924), name=chr, value=val);
#seg.dat2 <- data.frame(chr=chrN, start=rep(10, 5), end=c(11,11,11,11,11), name=chrN, value=val2);
seg.dat2 <- read.table("/home/vimal/Desktop/GeneCoord.txt",header=T)
seg.c   <- segAnglePo(seg.dat, seg=chr);

cols    <- "black";

par(mar=c(2, 2, 2, 2))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="");

circos(R=300,type="chr", cir=seg.c,  W=80 , print.chr.lab=TRUE, cex = 10)
circos(R=250, cir=seg.c, W=40, mapping=seg.dat2, type="b", col.v=5, col="black", B=T, cex=abs(seg.dat[,5])*1.5);
circos(R=220, cir=seg.c, W=40, mapping=seg.dat2, type="label2", col.v=5, col=cols, B=T, cex=0.8);

#circos(R=190, cir=seg.c, W=80, mapping=seg.dat, type="b", col.v=5, col=cols, B=F, lwd=abs(seg.dat[,5])*1.5);
#circos(R=110, cir=seg.c, W=80, mapping=seg.dat, type="b2", col.v=5, col=cols[c(1,7)], cutoff=0, B=T, lwd=2);

