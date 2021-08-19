#!/bin/Rscript

library(seqinr)

pc2vc <- read.table("TESTPCVCLIST")

for(i in 1:length(pc2vc$V1)) { 
p <- paste(pc2vc$V1[i], pc2vc$V2[i], "aligned.faa", sep=".")
q <- paste("TESTPUTDIR", pc2vc$V2[i], p, sep="/")

myseqs <- read.alignment(q, format = "fasta")
mat <- dist.alignment(myseqs, matrix = "identity")
changetotab <- as.matrix(mat)

o <- paste(pc2vc$V1[i], pc2vc$V2[i], "dist",sep=".")
outf <- paste("TESTPUTDIR", pc2vc$V2[i], o, sep="/")

write.table(changetotab, file=outf, sep="\t")

}
