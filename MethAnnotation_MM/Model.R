library(data.table)
library(entropy)
library(grDevices)
#library(stats)
library(parallel)
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("cummeRbund")
library(cummeRbund)


setwd("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DSS_Analysis/")
filesLE <- list.files(path = "/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DSS_Analysis/",pattern = "_LE")
filesBB <- list.files(path = "/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DSS_Analysis/",pattern = "_BB")
filesAN <- list.files(path = "/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DSS_Analysis/",pattern = "_Ancestral")
files <-c(filesAN,filesLE,filesBB)

numberOfCores <- 4
temp <- mclapply(files, fread, sep=" ")

ALL <- temp[[1]]
counter = 0
for (i in 2:length(temp)) {
  counter <<- counter + 1
  ALL <- merge(ALL, temp[[i]], by=c("chr", "pos"), all = T)
  setnames(ALL, c(head(names(ALL), -2), paste0('N', counter), paste0('X', counter)))
}
setnames(ALL,"N.x", "N0")
setnames(ALL,"X.x", "X0")

ALL$pos <- paste(ALL$chr, ".",ALL$pos, sep = "")
ALL <- subset(ALL, select = -c(chr) )
ALL$na_count <- apply(ALL, 1, function(x) sum(is.na(x)))
ALL_work = subset(ALL, ALL$na_count < 1)
ALL_work_subset = ALL_work
myfun <- function (y, v)
{
  return(r = y*100/v)
}

#ALL_work = subset(ALL_work , ALL_work$na_count < 10)
#LE_work_subset <-subset(LE_work, N.0 >= 30 | X.0 > 30)
ALL_work_subset[,M.0:=myfun(X0,N0)]
ALL_work_subset[,M.1:=myfun(X1,N1)]
ALL_work_subset[,M.2:=myfun(X2,N2)]
ALL_work_subset[,M.3:=myfun(X3,N3)]
ALL_work_subset[,M.4:=myfun(X4,N4)]
ALL_work_subset[,M.5:=myfun(X5,N5)]
ALL_work_subset[,M.6:=myfun(X6,N6)]
ALL_work_subset[,M.7:=myfun(X7,N7)]
ALL_work_subset[,M.8:=myfun(X8,N8)]
ALL_work_subset[,M.9:=myfun(X9,N9)]
ALL_work_subset[,M.10:=myfun(X10,N10)]
ALL_work_subset[,M.11:=myfun(X11,N11)]
ALL_work_subset[,M.12:=myfun(X12,N12)]
ALL_work_subset[,M.13:=myfun(X13,N13)]
ALL_work_subset[,M.14:=myfun(X14,N14)]
ALL_work_subset[,M.15:=myfun(X15,N15)]
ALL_work_subset[,M.16:=myfun(X16,N16)]
ALL_work_subset[,M.17:=myfun(X17,N17)]
ALL_work_subset[,M.18:=myfun(X18,N18)]
ALL_work_subset[,M.19:=myfun(X19,N19)]
ALL_work_subset[,M.20:=myfun(X20,N20)]
ALL_work_subset[,M.21:=myfun(X21,N21)]
ALL_work_subset[,M.22:=myfun(X22,N22)]
ALL_work_subset[,M.23:=myfun(X23,N23)]
ALL_work_subset[,M.24:=myfun(X24,N24)]
ALL_work_subset[,M.25:=myfun(X25,N25)]
ALL_work_subset[,M.26:=myfun(X26,N26)]
ALL_work_subset[,M.27:=myfun(X27,N27)]
ALL_work_subset[,M.28:=myfun(X28,N28)]
ALL_work_subset[,M.29:=myfun(X29,N29)]
ALL_work_subset[,M.30:=myfun(X30,N30)]
ALL_work_subset[,M.31:=myfun(X31,N31)]
#ALL_work_subset[,M.32:=myfun(X32,N32)]
#ALL_work_subset[,M.33:=myfun(X33,N33)]

#ALL_work_subset <- ALL_work_subset[, c("pos","M.0","M.1","M.2","M.3","M.4","M.5","M.6","M.7","M.8","M.9","M.10","M.12","M.13","M.14"
#                                       ,"M.15","M.16","M.17","M.18","M.19","M.20","M.21","M.22","M.23","M.24","M.25","M.26","M.27"
#                                       ,"M.28","M.29","M.30","M.31","M.32","M.33"), with=FALSE]

ALL_work_subset <- ALL_work_subset[, c("pos","M.0","M.1","M.2","M.3","M.4","M.5","M.6","M.7","M.8","M.9","M.10","M.12","M.13","M.14"
                                       ,"M.15","M.16","M.17","M.18","M.19","M.20","M.21","M.22","M.23","M.24","M.25","M.26","M.27","M.28","M.29","M.30","M.31"), with=FALSE]

ANS_work_subset <- ALL_work_subset[, c("pos","M.0","M.1","M.2","M.3","M.4","M.5","M.6"), with=FALSE]
ANS_work_subset_complete <- ANS_work_subset[, -c("pos"), with =FALSE]

LE_work_subset <- ALL_work_subset[, c("pos","M.7","M.8","M.9","M.10","M.12","M.13","M.14","M.15","M.16","M.17","M.18"), with=FALSE]
LE_work_subset_complete <- LE_work_subset[,-c("pos"), with =FALSE ]

BB_work_subset <- ALL_work_subset[, c("pos","M.19","M.20","M.21","M.22","M.23","M.24","M.25","M.26","M.27","M.28","M.29","M.30","M.31","M.32","M.33"), with=FALSE]
BB_work_subset_complete <- BB_work_subset[, -c("pos"), with =FALSE]



## print a message function
f.print.message <- function(x) { cat("=== JBNJK", format(Sys.time(), "%Y %b %d %X"), paste0("===JBNJK ", x,"\n")) }



T.density_Ans <- c("H", "L", "H", "L", "L", "L", "H")
T.density_LE <- c("M", "H", "H", "H", "H", "M", "H", "H", "H", "H", "L")
T.density <- c(T.density_Ans,T.density_LE)

test1.df <- as.data.frame(ALL_work_subset)
class(test.df)
head(test1.df)
test1.df <-melt(test.df, na.rm = TRUE, value.name = "pos")
test1.df$pos <- factor(test1.df$pos, levels=test.df$pos)
test1.df <- test1.df[order(test1.df$pos),]
test1.df$selection <- ifelse(test1.df$variable %in% c("M.0","M.1","M.2","M.3","M.4","M.5","M.6"), "ANS","LEBB")
test1.df$Tri <- rep(T.density ,3233472)
colnames(test1.df) <- c("pos","randomEffect","perMeth","selection", "Tri")
## specify the model used for the analysis of an individual cytosine
basicModel <- function(x){
  mod <- as.matrix(anova(lm(terms(perMeth~selection+Tri+selection*Tri, keep.order=TRUE), data = x)))
  out <- c(x[1,2],nrow(x), mean(x[,3], na.rm = TRUE), sd(x[,3], na.rm = TRUE), mod[1:3, "F value"], mod[1:3, "Pr(>F)"], mod["Residuals","Mean Sq"])
  return(out)
}

modelOutputNames <- c("pos", "sample", "numSamples", "PerMeth",  "sdPcent", # fixed
                      paste0(rep(c("F_", "P_", "NA"), each = 1), c("selection", "Tri", "selectionXTri")), "MS_res") # except for MS_res, this depends on the model

headers <- c("pos","randomEffect","perMeth","sample","selection","Tri")
toFactor <- c("selection", "Tri")
colsForPasteMean <- c("perMeth")

usedModel <- basicModel
outfileEnding <- ".mod"
myPath <- "/home/vimal/Desktop/"
numberOfCores <- 4
outfileName <- paste0(myPath, "Output")

f.print.message(paste0("processing ", outfileName, "\n"))

# run the models
f.print.message(paste0("processing ", myPath, "\n"))
meth <- test1.df

for (tf in toFactor) {meth[[tf]] <- as.factor(meth[[tf]])}
f.print.message(paste0("splitting ", myPath, "\n"))
indices <- split(1:nrow(meth), meth$pos) #### Pretty cool function
f.print.message(paste0("testing ", length(indices), " models.\n"))


all.tests <- mclapply(indices, function(x) usedModel(meth[x,]), mc.cores = 5)
out <- do.call("rbind", all.tests)
#colnames(out) <- c("Pos", "Num", "Samples", )
out$padjust <- p.adjust(out$P_selection, method = "fdr", n = ro(out))
#colnames(out) <- c(modelOutputNames,"padjust")
f.print.message(paste0("writing ", outfileName, "\n"))
write.table(out, outfileName, sep = '\t', quote = FALSE)

# add the group means
f.print.message(paste0("adding group means ", outfileName, "\n"))
out <- read.table(outfileName, quote = '', sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
if (length(colsForPasteMean) > 1) {
  meth$forMean <- paste0("meanWithin_", apply(meth[,colsForPasteMean], 1, function(x) paste(x, collapse = '_')))
} else {
  meth$forMean <- meth[[colsForPasteMean]]
}
allGroups <- unique(meth$forMean)
temp <- aggregate(meth$perMeth, by = list(chromPos = meth$pos, grp = meth$forMean)
                  , function(x) mean(x, na.rm = TRUE), simplify = TRUE)
for (curGroup in allGroups) {
  subTemp <- subset(temp, grp == curGroup)
  curGroup <- paste0("meanWithin_", curGroup)
  out[[curGroup]] <- NA
  out[subTemp$chromPos, curGroup] <- subTemp$x
}
f.print.message(paste0("writing ", outfileName, "\n"))
write.table(out, outfileName, sep = '\t', quote = FALSE)
f.print.message("finished.\n")















#ANS_work_subset_complete <- methAvg(ANS_work_subset_complete,7)
#LE_work_subset_complete <- methAvg(LE_work_subset_complete,7)
#BB_work_subset_complete <- methAvg(BB_work_subset_complete,7)
#smoothScatter(BB_work_subset_complete$Avg,LE_work_subset_complete$Avg,nrpoints = Inf,cex=3)
#lines(x=c(40,100), y=c(20,60))



######### Plotting Entropy


ANS_work_subset_complete <- entropySJ(ANS_work_subset_complete,7)
LE_work_subset_complete <- entropySJ(LE_work_subset_complete,7)
BB_work_subset_complete <- entropySJ(BB_work_subset_complete,7)

par(mfrow=c(1,2))
boxplot(ANS_work_subset_complete$Entropy,LE_work_subset_complete$Entropy,BB_work_subset_complete$Entropy)
plot(density(BB_work_subset_complete$Entropy), col="blue")
lines(density(LE_work_subset_complete$Entropy), col="red")
lines(density(ANS_work_subset_complete$Entropy), col="black")

smoothScatter(1:68705,sort(ANS_work_subset_complete$Entropy), colramp=colorRampPalette(my.cols))
smoothScatter(1:68705,sort(BB_work_subset_complete$Entropy), colramp=colorRampPalette(my.cols))
smoothScatter(1:68705,sort(LE_work_subset_complete$Entropy), colramp=colorRampPalette(my.cols))



############# Smooth Scatter
library(RColorBrewer)
k <- 100
my.cols <- rev(brewer.pal(k, "RdYlBu"))

Ans_LE <-read.table("Ans_LE_DML.smooth.all.plotInput.txt", header = T)
Ans_BB <-read.table("Ans_BB_DML.smooth.all.plotInput.txt", header = T)

Ans_LE$mu1 <-as.numeric(as.character(Ans_LE$mu1)) *100
Ans_LE$mu2 <-as.numeric(as.character(Ans_LE$mu2)) *100
Ans_BB$mu1 <-as.numeric(as.character(Ans_BB$mu1)) *100
Ans_BB$mu2 <-as.numeric(as.character(Ans_BB$mu2)) *100

summary(Ans_LE)
par(mfrow=c(2,3))
smoothScatter(ANS_work_subset_complete$Avg,BB_work_subset_complete$Avg,nrpoints = Inf,cex=8,col=adjustcolor( "black", alpha.f = 0.1))
smoothScatter(ANS_work_subset_complete$Avg,LE_work_subset_complete$Avg,nrpoints = Inf,cex=8,col=adjustcolor( "black", alpha.f = 0.1))
smoothScatter(BB_work_subset_complete$Avg,LE_work_subset_complete$Avg,nrpoints = Inf,cex=8,col=adjustcolor( "black", alpha.f = 0.1))

smoothScatter(Ans_LE$mu1,Ans_LE$mu2,nrpoints = Inf,cex=8,col=adjustcolor( "black", alpha.f = 0.1))
smoothScatter(Ans_BB$mu1,Ans_BB$mu2,nrpoints = Inf,cex=8,col=adjustcolor( "black", alpha.f = 0.1))
smoothScatter(Ans_LE$mu1 ,Ans_LE$mu2,nrpoints = Inf,cex=8,col=adjustcolor( "black", alpha.f = 0.1))










######### Function to get smoothScatterPlot
library(functional) 
smoothScatterPlot <- function(x,y) {
  plotMeth <- data.frame(x =x, y =y)
  plotMeth$diff = log2((plotMeth$x + 10)/(plotMeth$y+10))
  plotMeth$col <- ifelse(plotMeth$diff >= 2, "red", ifelse(plotMeth$diff < -2, "red","black"))
  plotMeth[apply(plotMeth, 1, Compose(is.finite, all)),]
  head(plotMeth)
  na.omit(plotMeth)
  return(plotMeth)
}


######### Function to get SJ entropy
entropySJ <- function(df.sj,colNum) {
  df.sj[df.sj < 30] <- 0
  df.sj[df.sj >= 30 && df.sj < 67 ] <- 1
  df.sj[df.sj >= 67] <- 2
  colNum <- integer(colNum)
  
  for (i in 1:dim(df.sj)[1]){
    A <- as.numeric(df.sj[i,1:7])
    b <- table(A)
    #print(i)
    #print(colNum)
    #print(class(colNum))
    df.sj[i, Entropy:= entropy.empirical(b, unit="log2")]
  }
  return(df.sj)
}


######### Average methylation of a populetion
methAvg <- function(df.sj,colNum) {
  #df.sj[df.sj < 30] <- 0
  #df.sj[df.sj >= 30 && df.sj < 67 ] <- 1
  #df.sj[df.sj >= 67] <- 2
  colNum <- integer(colNum)
  
  for (i in 1:dim(df.sj)[1]){
    A <- as.numeric(df.sj[i,1:colNum])
    df.sj[i, Avg:= mean(A,na.rm = T) + (((50 - mean(A,na.rm = T))/50) *  (rnorm(1)**2) )]
    #print(A)
    #print(mean(A,na.rm = T))
  }
  
  return(df.sj)
}





