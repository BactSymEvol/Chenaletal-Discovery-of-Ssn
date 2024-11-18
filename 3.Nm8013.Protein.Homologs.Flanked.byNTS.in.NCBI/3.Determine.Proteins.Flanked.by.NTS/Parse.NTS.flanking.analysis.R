# R script to analyse flanking analysis + make graphs.
# Author: Luke B Harrison
# Expected headings in input file: GENE,GENOME,NTS_5,NTS_3,NTS_LIKE_5,NTS_LIKE_3

flanking.anal<-read.csv("Flanking.analysis.NTS.counts.csv")

# Exclude Neisseriaceae
neisseria.to.exclude<-as.matrix(read.csv("Neisseriaceae.to.exclude.from.refseq.reference.csv",header=F))
EXCLUDE_NEISSERIA<-TRUE

genes<-levels(as.factor(flanking.anal$GENE))
whole.protein.list<-as.vector(as.matrix(read.csv("../../2.Nm8013.Proteins.flanked.by.NTS.Annotations/1.NM8013.Proteome/Nm8013proteome.filelist.txt",header=F)))

output.mat<-c()
## Now iterate over loci
for(i in 1:length(genes)) {
 output.v<-rep(NA,times=11)
 dim(output.v)<-c(1,11)
 colnames(output.v)<-c("NUM_HOMOLOGS","NTS_MEAN","NTS_PRESENCE","NTS_PRES_%","NTS_ANDPRESENCE","NTS_ANDPRES_%","NTSlike_MEAN","NTSlike_PRESENCE","NTSlike_PRES_%","NTSlike_ANDPRESENCE","NTSlike_ANDPRES_%")
 t.subset<-subset(flanking.anal,GENE==genes[i])
 # Now exlcude nesseria if desired
 if(EXCLUDE_NEISSERIA) {
  clean.genome<-unlist(lapply(strsplit(t.subset$GENOME,"_"),function(x){paste(x[1],x[2],sep="_")}))
  t.subset<-t.subset[which(!clean.genome %in% neisseria.to.exclude),]
 }
 rownames(output.v)<-genes[i]
 # Case where no homologs were found
 if(length(t.subset$NTS_5) == 0) {
  output.v[1,1]<-0
  output.mat<-rbind(output.mat,output.v)
  next
 }
 ## NTS
 output.v[1,1]<-length(t.subset$NTS_5)
 output.v[1,2]<-mean(t.subset$NTS_5+t.subset$NTS_3)
 output.v[1,3]<-sum((t.subset$NTS_5+t.subset$NTS_3) > 0)
 output.v[1,4]<-(sum((t.subset$NTS_5+t.subset$NTS_3) > 0)/length(t.subset$NTS_5))*100
 output.v[1,5]<-sum(t.subset$NTS_5 > 0 & t.subset$NTS_3 > 0)
 output.v[1,6]<-(sum(t.subset$NTS_5 > 0 & t.subset$NTS_3 > 0)/length(t.subset$NTS_5))*100
 
 ## NTS_LIKE
 output.v[1,7]<-mean(t.subset$NTS_LIKE_5+t.subset$NTS_LIKE_3)
 output.v[1,8]<-sum((t.subset$NTS_LIKE_5+t.subset$NTS_LIKE_3) > 0)
 output.v[1,9]<-(sum((t.subset$NTS_LIKE_5+t.subset$NTS_LIKE_3) > 0)/length(t.subset$NTS_LIKE_5))*100
 output.v[1,10]<-sum(t.subset$NTS_LIKE_5 > 0 & t.subset$NTS_LIKE_3 > 0)
 output.v[1,11]<-(sum(t.subset$NTS_LIKE_5 > 0 & t.subset$NTS_LIKE_3 > 0)/length(t.subset$NTS_LIKE_5))*100
 
 output.mat<-rbind(output.mat,output.v)

}


## Now add in the proteins for which there was no blast hit (and thus not included in the .csv)
prot.present<-unique(as.matrix(flanking.anal$GENE))
prot.not.present<-whole.protein.list[which(!whole.protein.list%in% prot.present)]
extra.prot<-cbind(rep(0,times=length(prot.not.present)),NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
row.names(extra.prot)<-prot.not.present
output.mat<-rbind(output.mat,extra.prot)

write.csv(output.mat,file="NTS.flanking.analysis.summary.csv")

## Generate a figure 
#par(mfrow=c(2,1))
pdf("NTS.flanking.analysis.summary.figure1.pdf",width=10,height=5.5)

## NTS
output.mat.clean<-output.mat[which(output.mat[,1] > 0),]
plot(output.mat.clean[,1],output.mat.clean[,6],ylim=c(0,20),xlim=c(1,10000),xlab="Number of homologous proteins (log)",ylab="Percentage of homologous proteins flanked by NTS elements",log="x",xaxt="n",pch=19,col="blue",cex=0.4,bty="l")
at.x <- outer(1:9, 10^(0:4))
lab.x <- ifelse(log10(at.x) %% 1 == 0, at.x, NA)
axis(1, at=at.x, labels=lab.x, las=1)
dev.off()

pdf("NTS.flanking.analysis.summary.figure2.pdf",width=10,height=5.5)
## NTS-like
output.mat.clean<-output.mat[which(output.mat[,1] > 0),]
plot(output.mat.clean[,1],output.mat.clean[,11],ylim=c(0,50),xlim=c(1,10000),xlab="Number of homologous proteins (log)",ylab="Percentage of homologous proteins flanked by NTS-like elements",log="x",xaxt="n",pch=19,col="blue",cex=0.4,bty="l")
at.x <- outer(1:9, 10^(0:4))
lab.x <- ifelse(log10(at.x) %% 1 == 0, at.x, NA)
axis(1, at=at.x, labels=lab.x, las=1)
dev.off()


