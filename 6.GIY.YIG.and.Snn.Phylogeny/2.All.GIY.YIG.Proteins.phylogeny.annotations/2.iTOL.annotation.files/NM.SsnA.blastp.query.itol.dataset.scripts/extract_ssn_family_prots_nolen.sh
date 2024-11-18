## Script to run blastp against all reference/representative genomes
## Author: Luke B Harrison
## ${1} = list of genome accession (GCF_aseembly) format

## PARAMETERS
EVALUE_CUTOFF_BLAST="1e-10"
HOMOLOGY_CUTOFF=60
MIN_ALIGNMENT_LENGTH=75
MAX_LENGTH=2000

GENOMES_DBDIR="../3.Nm8013.Protein.Homologs.Flanked.byNTS.in.NCBI/1.NCBI.Ref.Seq.Representative.DB/RefSeq_genomes_20230829/GCFs"
SSNA="/project/def-fveyrier/luke/SSna.Paper.Revisisions/SsnA.Nm.8013.faa"

SCRATCH_DIR="/tmp"

OUTDIR="."
OUTPUT_FILE="${OUTDIR}/NM_ssnA_homologs.stats.csv"

echo "GENOME,HOMOLOG_BLASTP" > $OUTPUT_FILE

COUNTER=1
NUM_GENOMES=$(wc -l ${1} | awk '{print $1}')
rm -f ${OUTDIR}/SSn_protein_homologs.faa

## Loop over genomes
cat ${1} | while read i || [[ -n $i ]];
do
	echo -n "Processing ${i}..... , #${COUNTER}/${NUM_GENOMES}..."
	GCF=$(echo ${i} | sed -e 's/\(GCF_[0-9]*\.[0-9]\)_.*/\1/g')
	echo -n "${i}," >> $OUTPUT_FILE
	rm -f ${SCRATCH_DIR}/${i}.temp.protlist.txt	
	## First blastp homologs
	blastp -query $SSNA -db ${GENOMES_DBDIR}/${GCF}/protein.faa -seg no -max_target_seqs 100 -max_hsps 1 -qcov_hsp_perc $MIN_ALIGNMENT_LENGTH -evalue ${EVALUE_CUTOFF_BLAST} -outfmt 6 2> /dev/null > ${SCRATCH_DIR}/${i}.blast.results.temp.tsv
 	rm -f ${SCRATCH_DIR}/${i}.temp.protlist.txt
        cat ${SCRATCH_DIR}/${i}.blast.results.temp.tsv | while read k || [[ -n $k ]];
	do
	 # Iterate over blast results
	 if [ -n "${k}" ]
	 then
	  PROT=$(echo $k | cut -d ' ' -f 2)
	  LEN=$(grep "${PROT}" ${GENOMES_DBDIR}/${GCF}/protein.faa.fai | cut -f 2)
          if awk "BEGIN {exit !(${LEN} >= ${MAX_LENGTH})}"
	  then
	   echo -n "found a protein but its length (${LEN}) is too long.."
	  else
	   #echo -n "found a homolog of sufficient length"
	   echo ${k} | cut -d ' ' -f 2 >> ${SCRATCH_DIR}/${i}.temp.protlist.txt
	  fi
	 fi
	done
	if [ -f "${SCRATCH_DIR}/${i}.temp.protlist.txt" ]
	then
	 echo -n "YES," >> $OUTPUT_FILE
	 seqkit grep --quiet --pattern-file ${SCRATCH_DIR}/${i}.temp.protlist.txt ${GENOMES_DBDIR}/${GCF}/protein.faa | sed -e "s/>/>${GCF}|/" >> ${OUTDIR}/SSn_protein_homologs.faa

	else 
	 echo -n "NO," >> $OUTPUT_FILE
	fi

	rm -f ${SCRATCH_DIR}/${i}.temp.protlist.txt ${SCRATCH_DIR}/${i}.blast.results.temp.tsv
	((COUNTER=COUNTER+1))
	echo " done."
done


