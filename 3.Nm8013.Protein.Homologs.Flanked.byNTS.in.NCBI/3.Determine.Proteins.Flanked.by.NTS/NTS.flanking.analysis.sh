#!/bin/bash
#######
## Script to calculate the presence and number of NTS regions in a ~10kb region of genomic sequence
## That is centered around an arbitrary homolog of a given protein
## Input: - list of NM8013 proteins from which to calculate flanking NTS regions of homologs
## Output: CSV file with genome, number of NTS in the region, and a column just describing presence/absence of a NTS in that genome

FLANK_DIR="../2.Extract.Flanking.Regions"
FASTA_SUFFIX="_refseq_rGCF_eference_genomes.20230829.full.corrected.txt_flanks_5000.fasta"
FLANK_SIZE=5000
NM8013_PROTEOME_LIST="../../2.Nm8013.Proteins.flanked.by.NTS.Annotations/1.NM8013.Proteome/Nm8013proteome.filelist.txt"

NTSSEQ="CGTCATTCCCGCG[AC]A[ACG]GCGGGAATC[CT][AG]G"
NTSVARSEQ="CGTCATACTCGGGNNNNNCCCGAGTATC"
NTS_LIKE="[CTA]-G-T-C-[AT]-[TC]-N(0,2)-[CT]-[CT]-[CGA]-[CG]-[CTG]-[GA]-[CAT]-N(0,5)-[AT]-N-[GAC]-[CG]-[GCT]-[GA]-G-N(0,2)-[AG]-N-C-[CT]-N-[TCG]"

OUTPUT_FILE="Flanking.analysis.NTS.counts.csv"

## Loop over genomes
echo "GENE,GENOME,NTS_5,NTS_3,NTS_LIKE_5,NTS_LIKE_3" > $OUTPUT_FILE

for i in $(cat $NM8013_PROTEOME_LIST)
do
	GENE_NAME="$i"
	echo -n "Calculating number of repeated elements (NTS and NTS_like) elements in 5prime and 3prime flanking regions of holomos of ${GENE_NAME}..."
	grep ">" ${FLANK_DIR}/${GENE_NAME}${FASTA_SUFFIX} | while read j || [[ -n $j ]]
	do
		GENOME_NAME=$(echo $j | cut -d ' ' -f 1 | sed -e 's/>//g')
		fastafetch -f ${FLANK_DIR}/${GENE_NAME}${FASTA_SUFFIX} -i ${FLANK_DIR}/${GENE_NAME}${FASTA_SUFFIX}.faidx -q ${GENOME_NAME} > /tmp/${GENE_NAME}_${GENOME_NAME}.fasta
		echo -n "${GENE_NAME}," >> $OUTPUT_FILE
		echo -n "${GENOME_NAME}," >> $OUTPUT_FILE
		SEQ_LEN=$(fastalength /tmp/${GENE_NAME}_${GENOME_NAME}.fasta | awk '{print $1}')
		FLANK_START=$(($SEQ_LEN - 1 - $FLANK_SIZE))
		echo -n "`fuzznuc -sequence /tmp/${GENE_NAME}_${GENOME_NAME}.fasta -pattern $NTSSEQ -sbegin1 0 -send ${FLANK_SIZE} -complement -pmismatch 1 -pname NTS -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`," >> $OUTPUT_FILE
		echo -n "`fuzznuc -sequence /tmp/${GENE_NAME}_${GENOME_NAME}.fasta -pattern $NTSSEQ -sbegin1 $FLANK_START -send $((SEQ_LEN - 1)) -complement -pmismatch 1 -pname NTS -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`," >> $OUTPUT_FILE
		echo -n "`fuzznuc -sequence /tmp/${GENE_NAME}_${GENOME_NAME}.fasta -pattern $NTS_LIKE -sbegin1 0 -send ${FLANK_SIZE} -complement -pmismatch 0 -pname NTSlike -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`," >> $OUTPUT_FILE 
		echo "`fuzznuc -sequence /tmp/${GENE_NAME}_${GENOME_NAME}.fasta -pattern $NTS_LIKE -sbegin1 $FLANK_START -send $((SEQ_LEN - 1)) -complement -pmismatch 0 -pname NTSDGEN -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`" >> $OUTPUT_FILE
		rm /tmp/${GENE_NAME}_${GENOME_NAME}.fasta
	done
	echo "done."
done

