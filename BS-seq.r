library("SRAdb")

# Specify several file paths and directories in which
# you would like to work in.
readsDir <- "/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/BS-seq/"
metaDB <- "/media/vimal/DATA_only/Read_DATA/SRAmetadb.sqlite"
mySamplesFile <- "/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/BS-seq/BS-seq.txt"
autoAnnotation <- "/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/BS-seq/annotation/automaticallyGeneratedAnnotation.csv"

# Check for a local copy of the database and
# download it if it does not exist.
if (!file.exists(metaDB)) { getSRAdbFile(dirname(metaDB)) }

# Open a connection to the database
sra_con <- dbConnect(SQLite(), metaDB)


rs <- getSRA( search_terms ='"bisulfite-seq" AND "Col" AND "Arabidopsis" AND "thaliana"',out_types = c('sample', 'experiment'), sra_con )



# Here we assume that you have a precompiled list of 
# samples/experiments/studies which are of interest to you.
# See Note 3 for full text searches in the SRA metadata.
# Read in the list of experiments you are interested in.
samples <- scan(mySamplesFile, what = "character")

# You may have a collection of experiments and individual runs. 
# However, only individual runs can be downloaded. Note that 
# the conversion of a certain SRA ID type to another requires
# the different query types not to be mixed.
sampleTypes <- unique(substr(samples, 1, 3))
splitSamples <- lapply(sampleTypes, function(x)
  grep(paste0("^", x), samples, value = TRUE))
conversions <- lapply(splitSamples, function(x)
  sraConvert(x, sra_con = sra_con ))
allSRAids <- do.call("rbind", conversions)

# Retrieve the remaining metadata to determine
# the sequencing platform, the read length, etc.
SRAdata <- lapply(allSRAids$run, function(x)
  getSRA(search_terms = x, out_types = "sra", sra_con))
SRAdata <- do.call("rbind", SRAdata)
rownames(SRAdata) <- SRAdata$run # for data access
SRAdata$approxRL <- SRAdata$bases/SRAdata$spots
summary(SRAdata$spots)    # the total number of reads
summary(SRAdata$approxRL) # the average read lengths
table(SRAdata$platform)   # the sequencing platform
SRAdata$sample_name[1:15] # the names of the first 15 samples
# The sample names are not always too informative
# and I recommend renaming them.

# Search for strand-specific samples with regular
# expression. Note that this gives no guarantee
# for strand specificity - in any case, this
# should be verified manually.
posExp <- "strand[[:space:]{0,1}|-]specific"
negExp <- paste0("no[[:alpha:]{0,1}][[:space:]{0,1}|-]", posExp)
strSpec <- apply(SRAdata, 1, function(x)
  length(grep(posExp, x)) > 0)
notStSp <- apply(SRAdata, 1, function(x)
  length(grep(negExp, x)) > 0)
strSpec[notStSp] <- FALSE
sum(strSpec)  # the number of strand specific libraries
SRAdata$strSpec <- as.numeric(strSpec)

# Try to download FASTQ files (this may fail).
# If necessary, check the availability with 
# the function getFASTQinfo(). For details
# see ?getFASTQinfo or the SRAdb manual


#getFASTQfile(allSRAids$run, sra_con, readsDir, 'ftp')

# Alternatively download the SRA files and convert them.
# If necessary, check the availability with 
# the function getSRAinfo(). For details
# see ?getSRAinfo or the SRAdb manual
getSRAfile(allSRAids$run, sra_con, readsDir, 'sra' )
commands <- sapply(allSRAids$run, function(x)
  paste0("fastq-dump.2.8.0 --gzip --split-files ",file.path(readsDir, x), ".sra"))
sapply(commands, system)

# Get a vector of all FASTQ file names and
# remove the reverse reads of the PE samples.
allFastq <- list.files(readsDir, ".fastq.gz")
revFastq <- list.files(readsDir, "_2.fastq.gz")
onlyForwardFastq <- setdiff(allFastq, revFastq)
sampleIDs <- gsub(".fastq", "", onlyForwardFastq)

# Create a table with information for the alignment
# and the preprocessing. Store the sampleID, read length,
# platform type, FASTQ file names, sample/library names, 
# strand specificity and additional attributes in a csv file. 
# Even though the sample/library names and the attributes
# are not complete for each sample, they may help renaming
# all samples. 
sortedSRAdata <- SRAdata[sampleIDs,]
tabForProcessing <- data.frame(
  sampleID   = sampleIDs,
  readLen    = round(sortedSRAdata$approxRL, 0),
  platform   = sortedSRAdata$platform,
  fastqFile  = onlyForwardFastq,
  #strandSpec = sortedSRAdata$strSpec,
  SRAsamName = sortedSRAdata$sample_name,
  SRAlibName = sortedSRAdata$library_name,
  SRAattrib  = sortedSRAdata$sample_attribute
)
write.csv(tabForProcessing, file.path(autoAnnotation))
