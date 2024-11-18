## R script 
## Author: Luke B Harrison
## 
## Input: output of bash script containg the following headings 
## Expected headings: GENE,NTS_5,NTS_3 in CSV
## With one protein per line and infomration on NTS motifs in 5' and 3' flanking regions

## Output: added column about whether it is flanked by NTS elements or now as well as list of protein so flanked

flanking.anal<-read.csv("../3.Extract.flanking.regions/Nm8013_proteome_NTS.flanking.analysis.csv",header=TRUE,row.names=1)
pseudo<-flanking.anal[,1]
flanking.anal<-flanking.anal[,-1]
flanking.anal<-cbind(pseudo,flanking.anal,apply(flanking.anal,1,function(x) { if(x[1] > 0 && x[2] > 0) return(1) else return(0) }))
colnames(flanking.anal)<-c("PSEUDO","NTS_5","NTS_5","NTS_FLANKED")
write.csv(flanking.anal,file="./Nm8013_proteome_NTS.flanking.analysis.annotated.csv")
proteins.flanked<-rownames(flanking.anal)[flanking.anal[,4] > 0]
write.table(proteins.flanked,file="Nm8013_proteome_cdnas_flankedby.NTS.csv",row.names=F,col.names=F,quote=F)
