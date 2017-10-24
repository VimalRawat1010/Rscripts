#source("https://bioconductor.org/biocLite.R")
#biocLite("TCC")

## The  current  methods implicitly  assume  a  balanced  DE,  wherein  
## the  numbers  of  highly  and  lowly  expressed  DE in  samples  are  (nearly)
## equal.   As  a  result,  methods  assuming  unbiased  DE  will 
## not  work  well  on  data  with  biased  DE.  
## Although  a  major  purpose  of  data  normalization is to detect such DE entities,  their
## existence themselves consequently interferes with their opportunity to be top-ranked.  
.libPaths("/home/vimal/R_Library/")
library(TCC)

### Case 1
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",iteration = 3,
                       FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)
head(result)


### Case 2
data(hypoData)
group <- c(1,2)
tcc <- new("TCC", hypoData[,c(1,4)], group)
tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",iteration = 3,
                       FDR = 0.1, floorPDEG = 0.05)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)
head(result)


### Filtering low-count genes (optional)
library(TCC)
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)
filter <- as.logical(rowSums(hypoData) > 0)
tcc <- new("TCC", hypoData[filter,], group)
#tcc <- filterLowCountGenes(tcc, low.count = 0)
dim(tcc$count)


#######################      
####   DEGES/TbT  #####
#######################
set.seed(1000)
data(hypoData)
samplesize <- 100
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)

### This  method  estimates  an  empirical  distribution  of  the  parameters  of  the
### NB distribution by bootstrapping from the input data.  
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "bayseq",iteration = 1,
                       samplesize =samplesize,FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors
tcc$DEGES$execution.time

## Above function is a wrapper for code below
set.seed(1000)
data(hypoData)
samplesize <- 100
group <- c(1, 1, 1, 2, 2, 2)
floorPDEG <- 0.05
FDR <- 0.1
### STEP 1 ###
d <- DGEList(count = hypoData, group = group)
d <- calcNormFactors(d)
norm.factors <- d$samples$norm.factors
norm.factors <- norm.factors / mean(norm.factors)
### STEP 2 ###
  cD <- new("countData", data = hypoData, replicates = group,
             groups = list(NDE = rep(1, length = length(group)), DE = group),
             libsizes = colSums(hypoData) * norm.factors)
cD <- getPriors.NB(cD, samplesize = samplesize, estimation = "QL", cl = NULL)
cD <- getLikelihoods(cD, pET = "BIC", cl = NULL)
result <- topCounts(cD, group = "DE", number = nrow(hypoData))
result <- result[order(result$rowID), ]
pval <- result$FWER.DE
qval <- result$FDR.DE
if (sum(qval < FDR) > (floorPDEG * nrow(hypoData))) {
  is.DEG <- as.logical(qval < FDR)
  } else {
    is.DEG <- as.logical(rank(pval, ties.method = "min") <= nrow(hypoData) * floorPDEG)
    }
### STEP 3 ###
d <- DGEList(count = hypoData[!is.DEG, ], group = group)
d <- calcNormFactors(d)
norm.factors <- d$samples$norm.factors * colSums(hypoData[!is.DEG, ]) /colSums(hypoData)
norm.factors <- norm.factors / mean(norm.factors)
norm.factors





#######################      
####  DEGES/edgeR ##### 200-400 times faster than DEGES/TbT
#######################
library(TCC)
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                        iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors




#######################      
####  DEGES/edgeR ##### MultiGroup 200-400 times faster than DEGES/TbT
#######################
library(TCC)
data(hypoData_mg)
group <- c(1, 2, 3, 4, 5, 6,7,8,9)
tcc <- new("TCC", hypoData_mg, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)
tcc$norm.factors
###  DE analysis  ###
set.seed(1000)
samplesize <- 100
tcc <- estimateDE(tcc, test.method = "bayseq",
                   FDR = 0.1, samplesize = samplesize)
result <- getResult(tcc, sort = TRUE)
head(result)
table(tcc$estimatedDEG)

##################################

ExpData <- read.table("/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/RNAseq.txt",
                     header= T, row.names=1)

ExpData.mat <- as.matrix(as.data.frame(lapply(ExpData, as.integer)))
ExpData.mat <- as.matrix(as.data.frame(lapply(ExpData, as.integer)))

ExpData.Labels <- read.table(
  "/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/RNAseqSampleInfo.txt",
  header= F, row.names=1)
levels(ExpData.Labels)
length(levels(ExpData.Labels))
group <- as.numeric(ExpData.Labels)

#numCores <- detectCores()
#cl <- makeCluster(4)
library(doParallel)
library(foreach)
library(plyr)
cl <- makeCluster(4)  # Use 3 cores
registerDoParallel(cl) # register these 3 cores with the "foreach" package




tcc <- new("TCC", ExpData.mat, group)
tcc <- filterLowCountGenes(tcc, low.count = 0)


tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 1, FDR = 0.1, floorPDEG = 0.05)

library(baySeq)
tcc <- estimateDE(tcc, test.method = "bayseq",FDR = 0.1, samplesize = samplesize)
result <- getResult(tcc, sort = TRUE)
head(result)
table(tcc$estimatedDEG)










