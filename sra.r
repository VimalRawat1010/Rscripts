library(SRAdb)
library(sra)


readsDir <- "/media/vimal/DATA_only/Read_DATA/"
metaDB <- "/media/vimal/DATA_only/Read_DATA/SRAmetadb.sqlite"
mySamplesFile <- "/media/vimal/DATA_only/Read_DATA/mySamples.txt"
autoAnnotation <- "/media/vimal/DATA_only/Read_DATA/annotation/automaticallyGeneratedAnnotation.csv"


if (!file.exists(metaDB)) { getSRAdbFile(dirname(metaDB)) }

sra_con <- dbConnect(SQLite(),metaDB)
sra_tables <- dbListTables(sra_con)

##### IMPORTANT for quesry !!!
sra_tables
dbListFields(sra_con,"sample")
#########################

dbGetQuery(sra_con,'PRAGMA TABLE_INFO(study)')
rs <- dbGetQuery(sra_con,"select * from study limit 1")
rs[, 1:3]


rs <- dbGetQuery(sra_con, paste( "select * from sample, experiment where",
                                 "experiment.library_strategy like 'RNA-Seq%' AND sample.scientific_name in ('Arabidopsis thaliana')",sep=" "))
rs[1:3,]

getTableCounts <- function(tableName,conn) 
                  { 
                    sql <- sprintf("select count(*) from %s",tableName) 
                    return(dbGetQuery(conn,sql)[1,1]) 
                    }

do.call(rbind,sapply(sra_tables[c(2,4,5,11,12)], getTableCounts, sra_con, simplify=FALSE))



rs <- getSRA( search_terms ='"RNA-seq" AND "Col" AND "Arabidopsis" AND "thaliana"',out_types = c('sample', 'experiment'), sra_con )
rs <- getSRA( search_terms ='"central cell"',out_types = c('sample', 'experiment'), sra_con )


############################3
con <- dbConnect(SQLite(),'SRAmetadb.sqlite')
query <- dbGetQuery(con,
                    paste(
                      "SELECT
                            sample.sample_accession,
                            sample.scientific_name,
                            experiment.experiment_accession,
                            experiment.library_strategy
                          FROM sample
                          JOIN experiment ON
                            sample.sample_accession = experiment.sample_accession
                          WHERE experiment.library_strategy in ('RNA-Seq') AND sample.scientific_name in ('Arabidopsis thaliana')
                          GROUP by sample.sample_accession
                          HAVING COUNT(DISTINCT experiment.library_strategy) = 1"
                    )
)

######################33

query <- dbGetQuery(con,
                    paste(
                      "SELECT *

                          FROM sample
                          JOIN experiment ON
                            sample.sample_accession = experiment.sample_accession
                          WHERE experiment.library_strategy in ('RNA-Seq') AND sample.scientific_name in ('Arabidopsis thaliana')
                          GROUP by sample.sample_accession
                          HAVING COUNT(DISTINCT experiment.library_strategy) = 1"
                    )
)


query[1]

