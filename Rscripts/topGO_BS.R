# try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("topGO")


library(topGO)
library(SparseM)
##collect gene names from biomart
mart <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
##Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
#examine result
head (GTOGO)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
#convert from table format to list format
geneID2GO <- by(GTOGO$go_id, GTOGO$ensembl_gene_id, function(x) as.character(x))
#examine result
head (geneID2GO)


all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
#int.genes <- sample(x = all.genes, size = 2000) # some random genes 
AnsLE.genes <- read.table("/home/vimal/Desktop/AnsLe.genes")
AnsBB.genes <- read.table("/home/vimal/Desktop/AnsBB.genes")
LEBB.genes <- read.table("/home/vimal/Desktop/LEBB.genes")

AnsLE.genes = AnsLE.genes[,1]
AnsBB.genes = AnsBB.genes[,1]
LEBB.genes = LEBB.genes[,1]

AnsLE.genes <- factor(as.integer(all.genes %in% AnsLE.genes))
AnsBB.genes<- factor(as.integer(all.genes %in% AnsBB.genes))
LEBB.genes<- factor(as.integer(all.genes %in% LEBB.genes))

names(AnsLE.genes) = all.genes
names(AnsBB.genes) = all.genes
names(LEBB.genes) = all.genes

go.ANS.LE.obj <- new("topGOdata", ontology='BP'
              , allGenes = AnsLE.genes
              , annot = annFUN.gene2GO
              , gene2GO = geneID2GO
)
go.ANS.BB.obj <- new("topGOdata", ontology='BP'
                    , allGenes = AnsBB.genes
                    , annot = annFUN.gene2GO
                    , gene2GO = geneID2GO
)
go.LE.BB.obj <- new("topGOdata", ontology='BP'
                    , allGenes = LEBB.genes
                    , annot = annFUN.gene2GO
                    , gene2GO = geneID2GO
)



results.Ans.LE <- runTest(go.ANS.LE.obj, algorithm = "elim", statistic = "fisher")
results.Ans.BB <- runTest(go.ANS.BB.obj, algorithm = "elim", statistic = "fisher")
results.LE.BB <- runTest(go.LE.BB.obj, algorithm = "elim", statistic = "fisher")

results.tab1 <- GenTable(object = go.ANS.LE.obj, elimFisher = results)
results.tab2 <- GenTable(object = go.ANS.BB.obj, elimFisher = results)
results.tab3 <- GenTable(object = go.LE.BB.obj, elimFisher = results)




