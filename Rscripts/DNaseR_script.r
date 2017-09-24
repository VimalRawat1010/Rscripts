
### DNaseR script running

library(DNaseR)
setwd("/biodata/dep_coupland/common/FT_to_Vimal")
bamfile <- "DNase-Col-R1.bam"
f <- system.file("extdata", bamfile, package="DNaseR",mustWork = TRUE)
dgf <- footprints(bam=f, chrN="chr1", chrL=34.9e6, p=1e-5, width=c(6,40), N=2e6)

head(dgf$footprint.events)
nrow(dgf$footprint.events)
setwd(owd)