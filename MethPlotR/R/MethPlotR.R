#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


rm(list=ls())

#file.list=list( system.file("extdata", "test1.myCpG.txt", package = "MethPlotR"),
#                system.file("extdata", "test2.myCpG.txt", package = "MethPlotR"),
#                system.file("extdata", "control1.myCpG.txt", package = "MethPlotR"),
#                system.file("extdata", "control2.myCpG.txt", package = "MethPlotR") )

.libPaths("/home/vimal/R_Library/")
setwd("/home/vimal/Software/Github/Latest/MethPlotR/")
library(methylKit)
library(bsseq)
library(pryr)
file.list=list( test1="extdata/A-01_LE", test2="extdata/A-03_LE",
                ctrl1="extdata/A-02_BB",ctrl2="extdata/A-04_BB")



#myobj=methRead(file.list,sample.id=list("LE1","LE2","BB1","BB2"),assembly="tair10",treatment=c(1,1,0,0))


dat1.1 <- read.table(as.character(file.list[1]), header=TRUE)
dat1.2 <- read.table(as.character(file.list[2]), header=TRUE)
dat2.1 <- read.table(as.character(file.list[3]), header=TRUE)
dat2.2 <- read.table(as.character(file.list[4]), header=TRUE)

listOfDF = list(dat1.1, dat1.2, dat2.1, dat2.2)
sampleNames =list("A-01_LE","A-03_LE","A-02_BB","A-04_BB")


vimalObj =methBedRead(file.list,sampleNames)
typeof(vimalObj)
otype(vimalObj)
#bsseqObj =makeBSseqData(listOfDF,sampleNames)
typeof(bsseqObj)
otype(bsseqObj)


getMethylationStats(myobj[[2]],plot=F,both.strands=F)



##################### Not so nice Way


##### More precise way: S3 method
methBedRead <- function(file.list,sampleNames)
{
  readObj =list()
  for (i in 1:length(file.list))
  {
    file = as.character(file.list[i])
    readObj[[i]] <- data.frame(read.table(file))
  }

    ## Set the name for the class
    class(readObj) <- append(class(readObj),"MethPlotR")
    return(readObj)
}


makeBSseqData <- function(dat, sampleNames) {
  n0 <- length(dat)
  
  if(missing(sampleNames))
    sampleNames <- paste("sample", 1:n0, sep="")
  
  alldat <- dat[[1]]
  if(any(alldat[,"N"] < alldat[,"X"], na.rm=TRUE))
    stop("Some methylation counts are greater than coverage.\n")
  ix.X <- which(colnames(alldat) == "X")
  ix.N <- which(colnames(alldat) == "N")
  colnames(alldat)[ix.X] <- "X1"
  colnames(alldat)[ix.N] <- "N1"
  
  
  if(n0 > 1) { ## multiple replicates, merge data
    for(i in 2:n0) {
      thisdat <- dat[[i]]
      if(any(thisdat[,"N"] < thisdat[,"X"], na.rm=TRUE))
        stop("Some methylation counts are greater than coverage.\n")
      
      ix.X <- which(colnames(thisdat) == "X")
      ix.N <- which(colnames(thisdat) == "N")
      colnames(thisdat)[c(ix.X,ix.N)] <- paste(c("X", "N"),i, sep="")
      alldat <- merge(alldat, thisdat, all=TRUE)
    }
  }
  
  ## make BSseq object
  ix.X <- grep("X", colnames(alldat))
  ix.N <- grep("N", colnames(alldat))
  alldat[is.na(alldat)] <- 0
  M <- as.matrix(alldat[,ix.X, drop=FALSE])
  Cov <- as.matrix(alldat[,ix.N, drop=FALSE])
  colnames(M) <- colnames(Cov) <- sampleNames
  
  result <- BSseq(chr=alldat$chr, pos=alldat$pos, M=M, Cov=Cov)
  
  result
}




##### Best way: S4 method

methBedRead <- setClass(
  #Setting Name
  "methBedRead",

)












MPD <- function() {


  }
