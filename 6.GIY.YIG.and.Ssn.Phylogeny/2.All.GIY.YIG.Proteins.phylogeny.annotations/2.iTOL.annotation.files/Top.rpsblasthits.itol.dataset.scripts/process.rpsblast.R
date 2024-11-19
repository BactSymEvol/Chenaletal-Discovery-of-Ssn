#### R script to process rpsblast output and return the top hit for each protein/genome pair
#### Author: Luke B Harrison
#### Note: if a specific hit is superseeded by the superfamily hit, report support the sepcific hit
#### Input argu
cdd.key<-read.csv("giy.yig.updated.smp.bitscores_newcolors_final.tsv",sep="\t",header=F)
prots<-read.csv("../../../1.Extract.All.GIY.YIG.Proteins/3.CDHIT.Cluster.Proteins/refseq_GCF_eference_genomes.20230829.all.giy.yig.updated.lessthan750.c06.GCFprot.list.txt",header=F,sep="\t")

output.table<-NULL
rpsblast.output<-read.csv("../../../1.Extract.All.GIY.YIG.Proteins/2.Extract.Proteins/all.giy.yig.updated.rpsblasthits.tsv",header=F,sep=" ")
unique.prots<-prots$V1

for(i in 1:length(unique.prots)) {
 #print(unique.prots[i])
 t.prot<-strsplit(unique.prots[i],"|",fixed=T)[[1]][2]
 #print(t.prot)
 t.table<-rpsblast.output[which(rpsblast.output[,1] %in% t.prot),]
 t.table<-t.table[-which(t.table[,2] %in% "CDD:198380"),]
 if(dim(t.table)[1]==0) {
  t.table<-rpsblast.output[which(rpsblast.output[,1] %in% t.prot),]
  #next
  
 }
 cat(unique.prots[i],"\t",sep="")
 top.hit<-t.table[which.min(t.table[,10]),]
 pssm<-gsub("CDD:","",top.hit[,2])
 cdd<-cdd.key[which(cdd.key[,1] %in% pssm),2]
 cdd.name<-cdd.key[which(cdd.key[,1] %in% pssm),4]
 colour<-cdd.key[which(cdd.key[,1] %in% pssm),5]
 cat(pssm,"\t",cdd,"\t",cdd.name,"\t",colour,"\n",sep="")
}

