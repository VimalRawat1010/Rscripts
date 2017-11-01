library(doParallel)result
library(foreach)
library(plyr)
library(baySeq)
library(TCC)

library(BiocParallel)
register(MulticoreParam(1))

ExpData <- read.table("/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/RNAseq.txt",
                      header= T, row.names=1)


filter <- as.logical(rowSums(ExpData) > 0)
ExpData.new <- ExpData[filter,]


ExpData.mat <- as.matrix(as.data.frame(lapply(ExpData.new, as.integer)))

ExpData.Labels <- read.table("/media/vimal/DATA_only/FG_Transcriptome_Project/VR/Results_Gene_Expression/RNAseqSampleInfo.txt",
                             header= F, row.names=1)
ExpData.Labels <- factor(ExpData.Labels$V2)
levels(ExpData.Labels)
length(levels(ExpData.Labels))
group <- as.numeric(ExpData.Labels)

names(group) <- colnames(ExpData)


tcc <- new("TCC", ExpData.mat, group)
tcc <- filterLowCountGenes(tcc, low.count = 0)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 3, FDR = 0.01, floorPDEG = 0.05)

tcc <- estimateDE(tcc, test.method = "edger",FDR = 0.1, samplesize = 100)
#tcc <- estimateDE(tcc, test.method = "deseq2",FDR = 0.1, samplesize = 32827)

result <- getResult(tcc, sort = FALSE)

filter <- as.logical(rowSums(ExpData) > 10)
ExpData.new <- ExpData[filter,]

result$gene_id <- rownames(ExpData.new)
result.sorted <- result[order(result$rank),]

head(result.sorted)
table(tcc$estimatedDEG)


install.packages("viridis")
library(viridis)


top100 <- result.sorted[result.sorted$rank < 5000,]
top100_heatmap <- ExpData[rownames(ExpData) %in% top100$gene_id,]
png("A.png")
heatmap.2(log10(data.matrix(top100_heatmap)+1), Colv=NA, scale="row", col=inferno(1000), ColSideColors=rainbow(29)[group] , trace= "none")
dev.off()
bottom100 <- result.sorted[result.sorted$rank > 28000,]
bottom100_heatmap <- ExpData[rownames(ExpData) %in% bottom100$gene_id,]
png("B.png")
heatmap.2(log10(data.matrix(bottom100_heatmap) +1), Colv=NA, scale="row", col=inferno(1000), ColSideColors=rainbow(29)[group] , trace= "none")
dev.off()

library(ggfortify)
pdf("PCA_TOP.pdf")
autoplot(prcomp(sqrt(data.matrix(top100_heatmap) `+1)))
dev.off()

pdf("PCA_BOTTOM.pdf")
autoplot(prcomp(sqrt(data.matrix(bottom100_heatmap) +1)))
dev.off()



