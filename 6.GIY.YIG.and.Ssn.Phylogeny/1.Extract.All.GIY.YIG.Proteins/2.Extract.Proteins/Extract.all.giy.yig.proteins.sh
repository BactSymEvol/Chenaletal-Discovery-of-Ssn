## Script to run rpsblast against all custom CDD database to extract ALL giy.yig proteins in the reference/respresentative genome database
## Author: Luke B Harrison
## ${1} = list of genome accession (GCF_aseembly) format
## Output: 1. RPSblast results (all.giy.yig.updated.rpsblasthits.tsv) 2: extracted GIY.YIG protein (refseq_GCF_eference_genomes.20230829.all.giy.yig.updated.faa)

## PARAMETERS
EVALUE_CUTOFF_RPSBLAST="0.1"

GENOMES_DBDIR="../3.Nm8013.Protein.Homologs.Flanked.byNTS.in.NCBI/1.NCBI.Ref.Seq.Representative.DB/RefSeq_genomes_20230829/GCFs"
CDD_DB="../1.Custom.CDD.Database/giy.yig.updated"
SSNA="../../../5.SRM.search.in.NCBI.refseq/1_Enumerate_SRM_in_NCBI/SsnA.Nm.8013.faa"
GIY_CDDs="../1.Custom.CDD.Database/giy.yig.updated.pssm.list.txt"
GIY_CDDs_BITSCORE="../1.Custom.CDD.Database/giy.yig.updated.smp.bitscores.lowered.tsv"

SCRATCH_DIR="/tmp"

OUTDIR="."
LOCAL_GIY_CDDs="${CDD_DB}/giy.yig.updated.pssm.list.txt"
LOCAL_GIY_CDDs_BITSCORE="${CDD_DB}/giy.yig.updated.smp.bitscores.lowered.tsv"

COUNTER=1
NUM_GENOMES=$(wc -l ${1} | awk '{print $1}')
rm -f ${OUTDIR}/RefSeq_genomes_20230829_all.giyyig.faa

## Loop over genomes
cat ${1} | while read i || [[ -n $i ]];
do
	echo -n "Processing ${i}..... , #${COUNTER}/${NUM_GENOMES}..."
	GCF=$(echo ${i} | sed -e 's/\(GCF_[0-9]*\.[0-9]\)_.*/\1/g')
	rpsblast -query ${GENOMES_DBDIR}/${GCF}/protein.faa -db $CDD_DB -seg no -evalue ${EVALUE_CUTOFF_RPSBLAST} -outfmt 6 > ${SCRATCH_DIR}/${i}.rpsblast.output.tsv
	rm -f ${SCRATCH_DIR}/${i}.temp.giyprotlist.txt


	cat $LOCAL_GIY_CDDs_BITSCORE | while read k || [[ -n $k ]];
	do
	 # iterate over CDDs
	 T_CDD=$(echo $k | cut -d ' ' -f2)
	 BIT_SCORE=$(echo $k | cut -d ' ' -f3)
	 PSSM=$(echo $k | cut -d ' ' -f1)
	 fgrep $PSSM ${SCRATCH_DIR}/${i}.rpsblast.output.tsv | while read j || [[ -n $j ]];
	 do
	  # Now iteratve over results for this CDD
	  SCORE=$(echo $j | cut -d ' ' -f 12)
	  if awk "BEGIN {exit !($SCORE >= $BIT_SCORE)}"
	  then
	   # yes we have foudn one
	   echo $j >> ${OUTDIR}/all.giy.yig.updated.rpsblasthits.tsv
	   echo $j | cut -d ' ' -f1 >> ${SCRATCH_DIR}/${i}.temp.giyprotlist.txt
	  fi
	 done
	done
        NUM_HITS=$(wc -l ${SCRATCH_DIR}/${i}.temp.giyprotlist.txt | cut -d ' ' -f 1)
	echo -n "found $NUM_HITS giy-yig proteins in genome ${i}..., extracting proteins..."
	if [ "$NUM_HITS" -gt 0 ]
	then
 	 seqkit grep --quiet --pattern-file ${SCRATCH_DIR}/${i}.temp.giyprotlist.txt ${GENOMES_DBDIR}/${GCF}/protein.faa | sed -e "s/>/>${GCF}|/" >> ${OUTDIR}/refseq_GCF_eference_genomes.20230829.all.giy.yig.updated.faa
	 echo -n "extracted.."
	else
	 echo -n "no hits..."
	fi
	((COUNTER=COUNTER+1))
	echo "done."
	rm -f ${SCRATCH_DIR}/${i}.rpsblast.output.tsv ${SCRATCH_DIR}/${i}.temp.giyprotlist.txt
done


