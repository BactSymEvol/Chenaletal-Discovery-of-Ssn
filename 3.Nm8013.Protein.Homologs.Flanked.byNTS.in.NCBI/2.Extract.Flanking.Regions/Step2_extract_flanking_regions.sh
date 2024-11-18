#!/bin/bash

## Script to extract flanking regions from the blast results from step1 of the script.
## This script will extract flanking genomic regions from hits that satisfies homology criteria
## Author: Luke B Harrison (luke.harrison@mail.mcgill.ca)
## This script is intended to be parallelized on a HPC cluster 

# PARAMETERS
#### Minimum percent identity to include in output
HOMOLOGY_CUTOFF=60

### At what identity pair-wise identity cutoff do we consider sequences identical (and thus dropped)
IDENTICAL_CUTOFF=0.9

## Minimum alignment length to the protein query
MIN_ALIGNMENT_LENGTH=75

## Size of flanking region to extract (bp in each direction)
FLANK_SIZE=5000

# Some basic environment / path parameters - I find using full path names is helpful on the DRCA clusters
GENOMES_DBDIR="../1.NCBI.Ref.Seq.Representative.DB/RefSeq_genomes_20230829/RefSeq_genomes_20230829"
GENOME_LIST="../1.NCBI.Ref.Seq.Representative.DB/refseq_GCF_eference_genomes.20230829.txt"
OUTDIR="."
BLAST_FILE="refseq_GCF_reference_genomes.20230828_blast.results.txt"
BLASTDB="."
NM8013_PROTEOME_LIST="../../2.Nm8013.Proteins.flanked.by.NTS.Annotations/1.NM8013.Proteome/Nm8013proteome.filelist.txt"

# Extract the name of this protein
GENOME_FILE=$(basename "${GENOME_LIST}")
echo "Copy blast results to local scratch"
cp ${BLASTDB}/${BLAST_FILE} ${SLURM_TMPDIR}
echo "done."

for z in $(cat ${NM8013_PROTEOME_LIST})
do
 GENE_NAME=$z

 echo -n "Searching for flanking regions of the following protein: "
 echo -n $GENE_NAME

 NUM_GENOMES=$(wc -l ${GENOME_LIST} | awk '{print $1}')
 echo " in a total of ${NUM_GENOMES} genomes"
 if test -f ${OUTDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}${IDENTICAL_CUTOFF}.fastaheaders.txt
 then
 	echo "Homologs already procssed...computed"
 else 
	## This for loop loops over all genomes, blasting the protein of interested against each one's CDNAs using tblastn
	COUNTER=1
	fgrep $GENE_NAME ${SLURM_TMPDIR}/${BLAST_FILE} > ${SLURM_TMPDIR}/temp.blast.${GENE_NAME}.txt
	for i in $(cat ${SLURM_TMPDIR}/temp.blast.${GENE_NAME}.txt | cut -f 2 | sed -e 's/\(GCF_.*\)|lcl.*/\1/g' | sort | uniq)
	do
		GCF=$i
		BLAST_OUT=$(grep $GCF ${SLURM_TMPDIR}/temp.blast.${GENE_NAME}.txt | head -n 1)
		if [ -n "$BLAST_OUT" ]
		then
			PER_IDENT=$(echo $BLAST_OUT | awk '{print $3}')
			if awk "BEGIN {exit !($PER_IDENT >= $HOMOLOGY_CUTOFF)}"
			then
				CDS=$(echo $BLAST_OUT | awk '{print $2}' | sed -e "s/${GCF}|//g" )
				fastafetch --fasta ${GENOMES_DBDIR}/${GCF}/cds_from_genomic.fna --index ${GENOMES_DBDIR}/${GCF}/cds_from_genomic.fna.faidx -q "$CDS" | sed -e "s/>/>${i} /" >> ${SLURM_TMPDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}.fasta
			fi
		fi
		 ((COUNTER=COUNTER+1))
	done
	## Now use cd-hit-est to generate sequence sluters and remove duplications with identity above thresholds
	if test -f ${SLURM_TMPDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}.fasta
	then
		echo -n "We have blast results.... now removing duplicates"
		cd-hit-est -i ${SLURM_TMPDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}.fasta -o ${SLURM_TMPDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}${IDENTICAL_CUTOFF}.fasta -c ${IDENTICAL_CUTOFF} -n 8 -d 0 -M 1024 -T 1
		rm ${SLURM_TMPDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}.fasta 
		grep ">" ${SLURM_TMPDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}${IDENTICAL_CUTOFF}.fasta > ${OUTDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}${IDENTICAL_CUTOFF}.fastaheaders.txt
		rm ${SLURM_TMPDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}${IDENTICAL_CUTOFF}.fasta
		echo "done"
	else
		echo "no blast results, therefore nothing to do, next protein"
		continue
	fi
fi

## Now loop over genomes where we have a hit that satisfies identity criteria
COUNTER=1
NUM_GENOMES=$(wc -l ${OUTDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}${IDENTICAL_CUTOFF}.fastaheaders.txt | cut -d " " -f 1)
echo "Extracing flanking regions in ${NUM_GENOMES} genomes"

#for i in $(cat ${OUTDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}${IDENTICAL_CUTOFF}.fastaheaders.txt)
cat ${OUTDIR}/${GENE_NAME}_${GENOME_FILE}_bacterial_DB_homology${HOMOLOGY_CUTOFF}${IDENTICAL_CUTOFF}.fastaheaders.txt | while read i || [[ -n $i ]];
do
	GENOME_NAME=`echo $i | awk '{print $1}' | sed -e 's/>//g'`
	GCF=$(echo $GENOME_NAME | sed -e 's/\(GCF_[0-9]*\.[0-9]\)_.*/\1/g')
	CONTIG=$(echo $i | sed -e 's/>.*lcl|\(.*\)_cds_.*/\1/g')
	LOCATION=$(echo $i | sed -e 's/.*\[location=\(.*\)\] .*/\1/g')
	if ( echo $i | grep -q "complement" )
	then
	 echo $LOCATION | sed -e 's/complement(//g' -e 's/)//g' -e 's/>//g' -e 's/<//g' | awk -v contig="$CONTIG" 'BEGIN { FS = "\\.\\." }; {print contig,$1,$2,"rev",1,"-";}' OFS='\t' | awk '{a=$2-1; $2 = a; print;}' OFS='\t' > ${SLURM_TMPDIR}/${CONTIG}_${GENE_NAME}.bed
	
	else 
	 echo $LOCATION | sed -e 's/>//g' -e 's/<//g' | awk -v contig="$CONTIG" 'BEGIN { FS = "\\.\\." }; {print contig,$1,$2,"fwd",1,"+";}' OFS='\t' | awk '{a=$2-1; $2 = a; print;}' OFS='\t' > ${SLURM_TMPDIR}/${CONTIG}_${GENE_NAME}.bed
 	fi
 	bedtools flank -i ${SLURM_TMPDIR}/${CONTIG}_${GENE_NAME}.bed -g ${GENOMES_DBDIR}/${GCF}/${GENOME_NAME}*genomic.fna.fai -b ${FLANK_SIZE} > ${SLURM_TMPDIR}/${GENOME_NAME}_${GENE_NAME}.${FLANK_SIZE}.bed
	bedtools getfasta -s -fi ${GENOMES_DBDIR}/${GCF}/${GENOME_NAME}*genomic.fna -bed ${SLURM_TMPDIR}/${GENOME_NAME}_${GENE_NAME}.${FLANK_SIZE}.bed -fo ${SLURM_TMPDIR}/${GENOME_NAME}_${GENE_NAME}.${FLANK_SIZE}.fasta
	HEADER=$(grep ">" ${SLURM_TMPDIR}/${GENOME_NAME}_${GENE_NAME}.${FLANK_SIZE}.fasta | tr '\n' '|' | sed -e 's/|>/|/g' -e 's/|$//g')
	echo "${HEADER}" | sed -e "s/>/>${GENOME_NAME} /g" >> ${OUTDIR}/${GENE_NAME}_${GENOME_FILE}_flanks_${FLANK_SIZE}.fasta
        grep -v ">" ${SLURM_TMPDIR}/${GENOME_NAME}_${GENE_NAME}.${FLANK_SIZE}.fasta | tr -d '\n' | sed -e 's/$/\n/g'>> ${OUTDIR}/${GENE_NAME}_${GENOME_FILE}_flanks_${FLANK_SIZE}.fasta
	rm ${SLURM_TMPDIR}/${GENOME_NAME}_${GENE_NAME}.${FLANK_SIZE}.bed ${SLURM_TMPDIR}/${CONTIG}_${GENE_NAME}.bed ${SLURM_TMPDIR}/${GENOME_NAME}_${GENE_NAME}.${FLANK_SIZE}.fasta
        ((COUNTER=COUNTER+1))
done
echo "done."
done


