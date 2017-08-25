



# perl GeneSynteny.pl ../Genome_Data/Athal_TAIR10/TAIR10_GFF3_genesOnly.gff
# perl GeneSynteny.pl ../Genome_Data/CR GFF file location
# perl GeneSyntenyCompiler.pl ../Gene_Duplication_Project/AL.synteny  ../Gene_Duplication_Project/AT.synteny ../FG_Transcriptome_Project/NETWORK/New_GeneList/CC_overlapping  ../Gene_Duplication_Project/AT_AL.out  > AL_AT.synteny_overlap.input

setwd("/media/vimal/DATA_only/Gene_Duplication_Project")
AL_AT <-read.table("AL_AT.synteny.input")
CR_AT <-read.table("CR_AT.synteny.input")
plot(AL_AT[,2],AL_AT[,4], col=AL_AT[,5], pch=19,cex=(AL_AT[,6]),main="Synteny Map", xlab="AL gene order", ylab="AT gene order")
text(AL_AT[,2],AL_AT[,4],labels=AL_AT[,7],cex=AL_AT[,6], col = c(AL_AT[,5]))
plot(CR_AT[,2],CR_AT[,4], col=CR_AT[,5], pch=19,cex=(CR_AT[,6]),main="Synteny Map", xlab="CR gene order", ylab="CR gene order")
text(CR_AT[,2],CR_AT[,4],labels=CR_AT[,7],cex=CR_AT[,6], col = c(CR_AT[,5]))