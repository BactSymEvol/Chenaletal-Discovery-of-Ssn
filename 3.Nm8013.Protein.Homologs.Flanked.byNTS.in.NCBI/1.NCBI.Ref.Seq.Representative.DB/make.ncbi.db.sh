#!/bin/bash

WORK_DIR="RefSeq_genomes_20230829"

for i in $(cat refseq_GCF_eference_genomes.20230829.txt)
do
GENOME_NAME=`echo $i | sed -e 's/.*\/\(GCF.*\)/\1/g'`
echo -n "$GENOME_NAME"
GCF=`echo ${GENOME_NAME} | sed -e 's/\(GCF_[0-9]*\.[0-9]\)_.*/\1/g'`

if test -f ${WORK_DIR}/${GCF}/protein.faa.fai
then
	echo -n "DB+index already created."
else 
        makeblastdb -dbtype nucl -in ${WORK_DIR}/${GCF}/${GENOME_NAME}_genomic.fna -out ${WORK_DIR}/${GCF}/${GENOME_NAME}_genomic.fna
        samtools faidx ${WORK_DIR}/${GCF}/${GENOME_NAME}_genomic.fna
	fastaindex -f ${WORK_DIR}/${GCF}/${GENOME_NAME}_genomic.fna -i ${WORK_DIR}/${GCF}/${GENOME_NAME}_genomic.fna.faidx
	makeblastdb -dbtype nucl -in ${WORK_DIR}/${GCF}/cds_from_genomic.fna -out ${WORK_DIR}/${GCF}/cds_from_genomic.fna
        samtools faidx ${WORK_DIR}/${GCF}/cds_from_genomic.fna
	fastaindex -f ${WORK_DIR}/${GCF}/cds_from_genomic.fna -i ${WORK_DIR}/${GCF}/cds_from_genomic.fna.faidx
	makeblastdb -dbtype prot -in ${WORK_DIR}/${GCF}/protein.faa -out ${WORK_DIR}/${GCF}/protein.faa
	samtools faidx ${WORK_DIR}/${GCF}/protein.faa
	fastaindex -f ${WORK_DIR}/${GCF}/protein.faa -i ${WORK_DIR}/${GCF}/protein.faa.faidx
fi
echo "done"
done
