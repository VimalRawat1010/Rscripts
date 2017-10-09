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


#file.list=list( system.file("extdata", "test1.myCpG.txt", package = "MethPlotR"),
#                system.file("extdata", "test2.myCpG.txt", package = "MethPlotR"),
#                system.file("extdata", "control1.myCpG.txt", package = "MethPlotR"),
#                system.file("extdata", "control2.myCpG.txt", package = "MethPlotR") )

.libPaths("/home/vimal/R_Library/")
setwd("/home/vimal/Software/Github/Latest/MethPlotR/")
library(methylKit)
file.list=list( test1="extdata/test1.myCpG.txt", test2="extdata/test2.myCpG.txt",
                ctrl1="extdata/control1.myCpG.txt",ctrl2="extdata/control2.myCpG.txt")



myobj=methRead( file.list,sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="tair10",treatment=c(1,1,0,0))

#typeof(myobj)
getMethylationStats(myobj[[2]],plot=F,both.strands=F)

sampleNames =list("test1","test2","ctrl1","ctrl2")


##################### Not so nice Way
methBedRead <- function(file.list,sampleNames) {
                obj =list()
                for (i in 1:length(file.list))
                {
                  file = as.character(file.list[i])
                  obj[[i]] <- data.frame(read.table(file))
                }
                return (obj)
}

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


##### Best way: S4 method

methBedRead <- setClass(
  #Setting Name
  "methBedRead",



)



  vimalObj =methBedRead(file.list,sampleNames)











MPD <- function() {


  }
