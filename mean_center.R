pheno <- snakemake@input[["pheno"]]
outfile <- snakemake@output[["pheno"]]


inf<-read.table (pheno) [, -c(4)]
inf[,3] <- scale (inf[,3])

write.table (inf, outfile, col.names=FALSE, row.names=FALSE,quote=FALSE,sep="\t")

