

A <- scan(file="RNA_header.txt",what=character())
B <- vector(mode="character", length=116)

for (i in 1:116)
{
q= paste("select * from sra_ft where run_accession='", b ,"' OR experiment_accession='", A[i] ,"'" ,sep="")
res <- dbGetQuery(sra_con, q)
B[i] = res[22]
}