## Short R script to process table of SRM counts in shuiffled genomes
## Author Luke B Harrison

shuffled.table<-read.csv("refseq_GCF_eference_genomes.20230829_shuffled.genomes.SRMs.tsv",sep="\t")
table.n<-shuffled.table[,-c(1,2)]
write.table(cbind(shuffled.table[,c(1,2)],rowMeans(table.n),(rowMeans(table.n)/shuffled.table[,2])*1000000),file="refseq_GCF_eference_genomes.20230829_shuffled.genomes.SRMs.means.tsv",sep="\t",quote=F,row.names=F,col.names=F)
