### Script to search for NTS-like elements in the flanking regions of representative/reference genome database
### Author: Luke B. Harrison
### Input: $1 a list of NCBI genome assembly IDs corresponding to the DB in question

## PARAMETERS
GENOMES_DBDIR="../3.Nm8013.Protein.Homologs.Flanked.byNTS.in.NCBI/1.NCBI.Ref.Seq.Representative.DB/RefSeq_genomes_20230829/GCFs"
NTSSEQ="CGTCATTCCCGCG[AC]A[ACG]GCGGGAATC[CT][AG]G"
NTS_LIKE="[CTA]-G-T-C-[AT]-[TC]-N(0,2)-[CT]-[CT]-[CGA]-[CG]-[CTG]-[GA]-[CAT]-N(0,5)-[AT]-N-[GAC]-[CG]-[GCT]-[GA]-G-N(0,2)-[AG]-N-C-[CT]-N-[TCG]"

FLANK_SIZE=500

OUTPUT_DIR="."
OUTPUT_FILE="${OUTPUT_DIR}/giy.yig.all.flanking.${FLANK_SIZE}.NTS.csv"
ITOL_OUTPUT="${OUTPUT_DIR}/giy.yig.all.flanking.${FLANK_SIZE}.NTS.foritol.csv"
echo "GENOME,PROT,NTS_5,NTS_3,NTS_FLANKED,NTS_LIKE_5,NTS_LIKE_3,NTS_LIKE_FLANKED" > $OUTPUT_FILE

COUNTER=1
NUM_GENOMES=$(wc -l ${1} | awk '{print $1}')

## Loop over genomes
cat ${1} | while read i || [[ -n $i ]];
do
	echo -n "Processing ${i}..... , #${COUNTER}/${NUM_GENOMES}..."
	GCF=$(echo ${i} | cut -f 1 -d '|')
	PROT=$(echo ${i} | cut -f 2 -d '|')
	echo -n "${GCF}," >> $OUTPUT_FILE
	echo -n "${PROT}," >> $OUTPUT_FILE
	## Find the CDNA that corresponds to this protein
	CDLINE=$(grep ${PROT} ${GENOMES_DBDIR}/${GCF}/cds_from_genomic.fna)
	GENOME="${GENOMES_DBDIR}/${GCF}/${GCF}_*_genomic.fna"
	## now extract the genomic region to calculate NTS numbers
	GENE_NAME=$(echo $CDLINE | sed -e 's/.*\[locus_tag=\([^]]*\).*/\1/g')
	CONTIG=$(echo $CDLINE | sed -e 's/>.*lcl|\(.*\)_cds_.*/\1/g')
	LOCATION=$(echo $CDLINE | sed -e 's/.*\[location=\(.*\)\] .*/\1/g')
	if ( echo $CDLINE | fgrep -q "[pseudo=")
	then
		PSEUDO=$(echo $i | sed -e 's/.*\[pseudo=\([^]]*\).*/\1/g')
	else
		PSEUDO="false"
	fi
	echo -n "Extracing flanking regions for ${GENE_NAME}, contig = ${CONTIG}, location = ${LOCATION}, PSEUDO = ${PSEUDO}"
	if ( echo $CDLINE | grep -q "complement" )
	then
        	echo $LOCATION | sed -e 's/complement(//g' -e 's/)//g' -e 's/>//g' -e 's/<//g' | awk -v contig="$CONTIG" 'BEGIN { FS = "\\.\\." }; {print contig,$1,$2,"rev",1,"-";}' OFS='\t' | awk '{a=$2-1; $2 = a; print;}' OFS='\t' > ${SLURM_TMPDIR}/${CONTIG}_${GENE_NAME}.bed
	else
        	echo $LOCATION | sed -e 's/>//g' -e 's/<//g' | awk -v contig="$CONTIG" 'BEGIN { FS = "\\.\\." }; {print contig,$1,$2,"fwd",1,"+";}' OFS='\t' | awk '{a=$2-1; $2 = a; print;}' OFS='\t' > ${SLURM_TMPDIR}/${CONTIG}_${GENE_NAME}.bed
	fi
	bedtools flank -i ${SLURM_TMPDIR}/${CONTIG}_${GENE_NAME}.bed -g ${GENOME}.fai -b ${FLANK_SIZE} > ${SLURM_TMPDIR}/${GENE_NAME}.${FLANK_SIZE}.bed
	bedtools getfasta -s -fi ${GENOME} -bed ${SLURM_TMPDIR}/${GENE_NAME}.${FLANK_SIZE}.bed -fo ${SLURM_TMPDIR}/${GENE_NAME}.${FLANK_SIZE}.fasta
	HEADER=$(grep ">" ${SLURM_TMPDIR}/${GENE_NAME}.${FLANK_SIZE}.fasta | tr '\n' '|' | sed -e 's/|>/|/g' -e 's/|$//g')
	echo "${HEADER}" | sed -e "s/>/>${CONTIG} /g" >> ${SLURM_TMPDIR}/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta
	grep -v ">" ${SLURM_TMPDIR}/${GENE_NAME}.${FLANK_SIZE}.fasta | tr -d '\n' | sed -e 's/$/\n/g'>> ${SLURM_TMPDIR}/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta
	rm ${SLURM_TMPDIR}/${GENE_NAME}.${FLANK_SIZE}.bed ${SLURM_TMPDIR}/${CONTIG}_${GENE_NAME}.bed ${SLURM_TMPDIR}/${GENE_NAME}.${FLANK_SIZE}.fasta
	#NTS  
	SEQEND=$((${FLANK_SIZE} * 2))
	NTS_5=`fuzznuc -sequence ${SLURM_TMPDIR}/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta -pattern $NTSSEQ -sbegin1 0 -send $((${FLANK_SIZE} - 1)) -complement -pmismatch 1 -pname NTS -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`
	NTS_3=`fuzznuc -sequence ${SLURM_TMPDIR}/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta -pattern $NTSSEQ -sbegin1 ${FLANK_SIZE} -send $((${SEQEND} -1)) -complement -pmismatch 1 -pname NTS -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`
	echo -n "${NTS_5},${NTS_3}," >> $OUTPUT_FILE
	if [[ $NTS_5 -gt 0 && $NTS_3 -gt 0 ]]
	then
		echo -n "YES," >> $OUTPUT_FILE
	else
		echo -n "NO," >> $OUTPUT_FILE
	fi
	# NTS_like
	NTSLIKE_5=`fuzznuc -sequence ${SLURM_TMPDIR}/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta -pattern $NTS_LIKE -sbegin1 0 -send $((${FLANK_SIZE} - 1)) -complement -pmismatch 1 -pname NTS -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`
        NTSLIKE_3=`fuzznuc -sequence ${SLURM_TMPDIR}/${GENE_NAME}_flanks_${FLANK_SIZE}.fasta -pattern $NTS_LIKE -sbegin1 ${FLANK_SIZE} -send $((${SEQEND} -1)) -complement -pmismatch 1 -pname NTS -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`
        echo -n "${NTSLIKE_5},${NTSLIKE_3}," >> $OUTPUT_FILE
        if [[ $NTSLIKE_5 -gt 0 && $NTSLIKE_3 -gt 0 ]]
        then
                echo "YES" >> $OUTPUT_FILE
		echo "${GCF}|${PROT} #FF0000" >> $ITOL_OUTPUT
        else
                echo "NO" >> $OUTPUT_FILE
        fi
        ((COUNTER=COUNTER+1))
	echo "done."
done


