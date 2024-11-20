## R Script to plot NTS_like in relation to homology
# Author: Luke B Harrison (luke.harrison@mail.mcgill.ca)
#

nts.homology<-read.csv("SsnA_NTS_homologs.PSI.BLAST.refseq_reference_genomes_flanks.csv")
### Make a graph with bins for the genomes in 100 NTS_like windwos
bin_1_100<-subset(nts.homology,NTS_LIKE <= 100)
bin_101_200<-subset(nts.homology,NTS_LIKE <= 200 & NTS_LIKE > 100)
bin_201_plus<-subset(nts.homology,NTS_LIKE > 200)
nts.homology.matrix<-cbind(c((sum(bin_1_100$HOMOLOG60 == "YES")/length(bin_1_100$HOMOLOG60)),(sum(bin_1_100$HOMOLOG == "YES")/length(bin_1_100$HOMOLOG)),(sum(bin_1_100$HOMOLOG_PSI == "YES")/length(bin_1_100$HOMOLOG_PSI)),(sum(bin_1_100$HOMOLOG_CDD == "YES")/length(bin_1_100$HOMOLOG_CDD))),c((sum(bin_101_200$HOMOLOG60 == "YES")/length(bin_101_200$HOMOLOG60)),(sum(bin_101_200$HOMOLOG == "YES")/length(bin_101_200$HOMOLOG)),(sum(bin_101_200$HOMOLOG_PSI == "YES")/length(bin_101_200$HOMOLOG_PSI)),(sum(bin_101_200$HOMOLOG_CDD == "YES")/length(bin_101_200$HOMOLOG_CDD))),c((sum(bin_201_plus$HOMOLOG60 == "YES")/length(bin_201_plus$HOMOLOG60)),(sum(bin_201_plus$HOMOLOG == "YES")/length(bin_201_plus$HOMOLOG)),(sum(bin_201_plus$HOMOLOG_PSI == "YES")/length(bin_201_plus$HOMOLOG_PSI)),(sum(bin_201_plus$HOMOLOG_CDD == "YES")/length(bin_201_plus$HOMOLOG_CDD))))*100
colnames(nts.homology.matrix)<-c("0-99","100-199","> 200")
rownames(nts.homology.matrix)<-c("Tblastn_60","Tblastn","PSIblast","RPRblast")


svg("NTSlike.vs.frac.homologs.svg",width=10,height=10)
baseplot<-barplot(height=nts.homology.matrix,beside=T,legend=T,xlab="Number of NTS-like elements",ylab="Fraction of Genomes with an ssnA homolog",ylim=c(0,100))
dev.off()

