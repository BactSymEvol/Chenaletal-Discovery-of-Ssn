
## Script to run rpsblast against all custom CDD database to extract ALL giy.yig proteins in the reference/respresentative genome database
## Author: Luke B Harrison
## ${1} = list of genome accession (GCF_aseembly) format

## PARAMETERS
EVALUE_CUTOFF_RPSBLAST="0.1"

GENOMES_DBDIR="../3.Nm8013.Protein.Homologs.Flanked.byNTS.in.NCBI/1.NCBI.Ref.Seq.Representative.DB/RefSeq_genomes_20230829/GCFs"
CDD_DB="/home/luke/Dropbox/Postdoc.INRS/Projects/SsnA.Paper.2023/Data.supplement/6.GIY.YIG.and.Snn.Phylogeny/1.Extract.All.GIY.YIG.Proteins/1.Custom.CDD.Database"
SSNA="/project/def-fveyrier/luke/SSna.Paper.Revisisions/SsnA.Nm.8013.faa"
CDD_PSSM="/project/def-fveyrier/luke/SSna.Paper.Revisisions/cd10448.rescaled.smp"
GIY_CDDs="/project/def-fveyrier/luke/SSna.Paper.Revisisions/giy.yig.cdd.database.updated/giy.yig.updated.pssm.list.txt"
GIY_CDDs_BITSCORE="/project/def-fveyrier/luke/SSna.Paper.Revisisions/giy.yig.cdd.database.updated/giy.yig.updated.smp.bitscores.lowered.tsv"

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


