#!/bin/bash
## Author: Luke B Harrison
## Script to extract flanking regions for a given set of proteins in a given genome using tblastn
## And then calculate the presence and number of NTS regions in a ~10kb region of genomic sequence centered around a given protein
## PARAMETERS
## Size of flanking region to extract (bp in each direction)
FLANK_SIZE=5000

NTSSEQ="CGTCATTCCCGCG[AC]A[ACG]GCGGGAATC[CT][AG]G"
NTSVARSEQ="CGTCATACTCGGGNNNNNCCCGAGTATC"
## INPUTS
# $1 = The list of protein fasta files to extract flanking regions for (one per line, fully qualify names)
# Note the directory structures below

CDS="../1.NM8013.cDNAs/Nm8013_GCA_000026965.1_CDNAs.fastaheaders.txt"
GENOME="../2.NM8013.Genome/GCA_000026965.1_ASM2696v1_genomic.fna"

## OUTPUT:
## .csv file with the NTS statistics
OUTDIR="."
OUTPUT_FILE="${OUTDIR}/Nm8013_proteome_NTS.flanking.analysis.csv"
echo "GENE,PSEUDO,NTS_5,NTS_3" > $OUTPUT_FILE

cat ${CDS} | while read i || [[ -n $i ]];
do
 # Some basic environment / path parameters - I find using full path names is helpful on the DRCA clusters
 # Extract the gene name from the protein fasta input file
 GENE_NAME=$(echo $i | sed -e 's/.*\[locus_tag=\([^]]*\).*/\1/g')
 CONTIG=$(echo $i | sed -e 's/>.*lcl|\(.*\)_cds_.*/\1/g')
 LOCATION=$(echo $i | sed -e 's/.*\[location=\(.*\)\] .*/\1/g')
 if ( echo $i | fgrep -q "[pseudo=") 
 then
  PSEUDO=$(echo $i | sed -e 's/.*\[pseudo=\([^]]*\).*/\1/g')
 else 
  PSEUDO="false"
 fi
 echo -n "Extracing flanking regions for ${GENE_NAME}, contig = ${CONTIG}, location = ${LOCATION}, PSEUDO = ${PSEUDO}"
 echo -n "${GENE_NAME}," >> $OUTPUT_FILE
 echo -n "${PSEUDO}," >> $OUTPUT_FILE

 if ( echo $i | grep -q "complement" )
 then
 	echo $LOCATION | sed -e 's/complement(//g' -e 's/)//g' -e 's/>//g' -e 's/<//g' | awk -v contig="$CONTIG" 'BEGIN { FS = "\\.\\." }; {print contig,$1,$2,"rev",1,"-";}' OFS='\t' | awk '{a=$2-1; $2 = a; print;}' OFS='\t' > ${CONTIG}_${GENE_NAME}.bed
 else 
	echo $LOCATION | sed -e 's/>//g' -e 's/<//g' | awk -v contig="$CONTIG" 'BEGIN { FS = "\\.\\." }; {print contig,$1,$2,"fwd",1,"+";}' OFS='\t' | awk '{a=$2-1; $2 = a; print;}' OFS='\t' > ${CONTIG}_${GENE_NAME}.bed
 fi
 bedtools flank -i ${CONTIG}_${GENE_NAME}.bed -g ${GENOME}.fai -b ${FLANK_SIZE} > ${GENE_NAME}.${FLANK_SIZE}.bed
 bedtools getfasta -s -fi ${GENOME} -bed ${GENE_NAME}.${FLANK_SIZE}.bed -fo ${GENE_NAME}.${FLANK_SIZE}.fasta
 HEADER=$(grep ">" ${GENE_NAME}.${FLANK_SIZE}.fasta | tr '\n' '|' | sed -e 's/|>/|/g' -e 's/|$//g')
 echo "${HEADER}" | sed -e "s/>/>${CONTIG} /g" >> ${OUTDIR}/flanks/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta
 grep -v ">" ${GENE_NAME}.${FLANK_SIZE}.fasta | tr -d '\n' | sed -e 's/$/\n/g'>> ${OUTDIR}/flanks/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta
 rm ${GENE_NAME}.${FLANK_SIZE}.bed ${CONTIG}_${GENE_NAME}.bed ${GENE_NAME}.${FLANK_SIZE}.fasta
 #NTS  
 SEQEND=$((${FLANK_SIZE} * 2))

 echo -n "`fuzznuc -sequence ${OUTDIR}/flanks/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta -pattern $NTSSEQ -sbegin1 0 -send $((${FLANK_SIZE} - 1)) -complement -pmismatch 1 -pname NTS -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`," >> $OUTPUT_FILE
 echo "`fuzznuc -sequence ${OUTDIR}/flanks/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta -pattern $NTSSEQ -sbegin1 ${FLANK_SIZE} -send $((${SEQEND} -1)) -complement -pmismatch 1 -pname NTS -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`" >> $OUTPUT_FILE
 echo " done."
done



