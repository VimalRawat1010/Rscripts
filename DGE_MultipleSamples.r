#### Diff Gene Exp analysis (RNA-seq data)  for multiple samples
#Loading stored data
data("simData")
if(require("parallel")) cl <- makeCluster(6) else cl <- NULL
# Defining data: Column headers basically
replicates <-c("simA", "simA", "simA", "simA", "simA", "simB", "simB", "simB", "simB", "simB")

# Defining model for NDE & DE genes
# The key to constructing vectors corresponding to a model is to see for which groups of libraries 
# we expect equivalent expression of tags
# The ultimate aim of the baySeq package is to evaluate posterior likelihoods of each model for each row of the data.
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1),DE = c(1,1,1,1,1,2,2,2,2,2))

# Combining the count data and user-dened groups into a countData object
CD <- new("countData", data = simData, replicates = replicates, groups = groups)

#Library sizes can be inferred from the data
libsizes(CD) <- getLibsizes(CD)

# MA-plot : data  are uniformly zero (and hence the log-ratio is innite)  
plotMA.CD(CD, samplesA = "simA", samplesB = "simB",col = c(rep("red", 100), rep("black", 900)))


# Adding annotation
CD@annotation <- data.frame(name = paste("count", 1:1000, sep = "_"))

# Negitive Bionomial parameter Estimation
CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)

# posterior likelihoods, estimating the proportions of dierentially expressed counts
CD <- getLikelihoods(CD, cl = cl, bootStraps = 3, verbose = FALSE)
CD@estProps

CD@posteriors[1:10,]
CD@posteriors[101:110,]


# Top candidates for differential expression using the topCounts function
topCounts(CD, group = "DE")


plotPosteriors(CD, group = "DE", col = c(rep("red", 100), rep("black", 900)))