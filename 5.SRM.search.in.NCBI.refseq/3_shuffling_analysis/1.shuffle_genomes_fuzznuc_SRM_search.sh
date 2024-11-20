### Script to shuggle genomes of the bacterial reference / representative database and conduct NTS/SRM searchs on the shuffled genomes
### Author: Luke B Harrison
### Input: control file with genome identifiers (GCF_Assembly format)
### Output: 
## PARAMETERS
## kmer size = 2 (di-nucleotide shuffle)

GENOMES_DBDIR="../3.Nm8013.Protein.Homologs.Flanked.byNTS.in.NCBI/1.NCBI.Ref.Seq.Representative.DB"

NTS_LIKE="[CTA]-G-T-C-[AT]-[TC]-N(0,2)-[CT]-[CT]-[CGA]-[CG]-[CTG]-[GA]-[CAT]-N(0,5)-[AT]-N-[GAC]-[CG]-[GCT]-[GA]-G-N(0,2)-[AG]-N-C-[CT]-N-[TCG]"
KMER_LEN=2
N_SHUFFLES=100
SCRATCH_DIR="/tmp"

OUTDIR="."

COUNTER=1
NUM_GENOMES=$(wc -l ${1} | awk '{print $1}')
OUTPUT_FILE=${OUTDIR}/NCBI.Ref.Seq.Representative.DB.shuffled.genomes.SRM.counts.tsv

echo -n "GENOME	LENGTH" >> $OUTPUT_FILE
for j in $(seq 1 ${N_SHUFFLES})
do
	echo -n "	SRM_SHUF_${j}" >> $OUTPUT_FILE
done
echo -en "\n" >> $OUTPUT_FILE

## Loop over genomes
cat ${1} | while read i || [[ -n $i ]];
do
	echo -n "Processing ${i}..... , #${COUNTER}/${NUM_GENOMES}..."
	echo -n ${i} >> $OUTPUT_FILE
	GCF=$(echo $i | cut -d '_' -f 1,2 )
	cp ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna ${SCRATCH_DIR}
	echo -n "	`seqkit stats -T ${SCRATCH_DIR}/${i}_genomic.fna | cut -f 5 | tail -n +2`" >> $OUTPUT_FILE
	for j in $(seq 1 ${N_SHUFFLES})
	do
		fasta-shuffle-letters -kmer ${KMER_LEN} -dna -seed $RANDOM ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna > ${SCRATCH_DIR}/${GCF}_shuffled_kmer${KMER_LEN}_copy${j}.fna
		echo -n "	`fuzznuc -sequence ${SCRATCH_DIR}/${GCF}_shuffled_kmer${KMER_LEN}_copy${j}.fna -pattern $NTS_LIKE -complement -pmismatch 0 -pname NTSLIKE -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`" >> $OUTPUT_FILE
		rm ${SCRATCH_DIR}/${GCF}_shuffled_kmer${KMER_LEN}_copy${j}.fna
	done
	echo -en "\n" >> $OUTPUT_FILE
	((COUNTER=COUNTER+1))
	echo "done."
done


