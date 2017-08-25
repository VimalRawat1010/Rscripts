### Diff gene expression Analysis using DESeq : for Wild type & Mutant samples

#===============================================================================
#
#         FILE: DESeq_pipeline.r
#
#        USAGE: Rscript DESeq_pipeline.r
#
#  DESCRIPTION:  Script needs input file to be stored in working directory, format of file is 
# 		INF.gene.count.final is a tab seprated file : File format
# 		
#		Gene_Id Read_count_WT_1 Read_count_WT_2 Read_count_Mutant_1 Read_count_Mutant_2 
# 		AT1G22770 10 9 209 310
# 		AT2G33770 5 11 80 110 			      
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Vimal Rawat (vimal.biochem@gmail.com or vimal.rawat@botinst.uzh.ch),
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 22/02/2016
#     REVISION: ---
#===============================================================================

### Clear the workspace
rm(list =ls())

### Set path to your R default library (or your personal R library)
.libPaths("~/R_Library/")

### I defined a function to check and install if some R packages are not present which are needed for this code to run.
install_pack <-  function(package_list) {
  # package_list: is a list of packages needs to be installed
  # Returns nothing :)
    for (package in package_list) {
     if(!require(package,character.only=TRUE)) {
       install.packages(package , dependencies = TRUE)
       suppressMessages(suppressWarnings(require(package)))
      }
    }
}

### Calling  install_pack function
install_pack(c("DESeq"))



###  set working directory : Directory where input files are located and plots and result files will be saved
setwd("~/DESeqResults/") # Replace to desired path

### Load DESeq package : install the package if you havn't installed it already with "install" command
library(DESeq)

### Load input files 
counts = read.delim("INF.gene.count.final", header=F, row.names=1)

### Defining experimental design 
my.design <- data.frame(
  row.names = colnames( counts ),
  condition = c( "wt", "wt", "mut", "mut" ),
  libType = c( "single-end", "single-end", "single-end", "single-end" )
)
conds <- factor(my.design$condition)
CDS <- newCountDataSet( counts, conds )
CDS <- estimateSizeFactors( CDS )

# Calculating size factor for library : Libraray size bias
sizeFactors( CDS )
CDS <- estimateDispersions( CDS )

# Saving PDF file for Dispersion plot
pdf("INF_WT_fdfdp_DESeq-dispersion_estimates.pdf")

# Generating Dispersion plot
plot(
  rowMeans( counts( CDS, normalized=TRUE ) ),
  fitInfo(CDS)$perGeneDispEsts,
  pch = '.', log="xy"
)
xg <- 10^seq( -.5, 5, length.out=300 )
lines( xg, fitInfo(CDS)$dispFun( xg ), col="red" )

# Saving Dispersion plot
dev.off()

# Negative binomial test for differential gene expression : Core analysis
result <- nbinomTest( CDS, "wt", "mut" )
result = result[order(result$pval), ]

# Writting Results to a csv file
write.csv(result, "DESeq-WT-vs-fdfdp.csv")

# Plotting Results to a pdf file
pdf("INF_WT_fdfdp_DESeq-MA-plot.pdf")
plot(
  result$baseMean,
  result$log2FoldChange,
  log="x", pch=20, cex=.3,
  col = ifelse( result$padj < .1, "red", "black" ) )
dev.off()


##### Now relax !! You are done with Analysis
##### Publish & become famous !!!


