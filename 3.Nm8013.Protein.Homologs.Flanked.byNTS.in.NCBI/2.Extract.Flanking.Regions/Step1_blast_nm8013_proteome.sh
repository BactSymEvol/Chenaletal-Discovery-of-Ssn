#!/bin/bash

## Script to call tblastn for the Nm8013 proteome against the NCBI reference/represenative database
## Author: Luke B Harrison (luke.harrison@mail.mcgill.ca)
## This script is intended to be parallelized on a HPC cluster
## This is step one of the proceedure to extract flanking regions

# PARAMETERS
#### Minimum percent identity to include in output
HOMOLOGY_CUTOFF=60

## Minimum alignment length to the protein query
MIN_ALIGNMENT_LENGTH=75

## Size of flanking region to extract (bp in each direction)
FLANK_SIZE=5000

# Some basic environment / path parameters - I find using full path names is helpful on the DRCA clusters
GENOMES_DBDIR="../1.NCBI.Ref.Seq.Representative.DB/RefSeq_genomes_20230829/RefSeq_genomes_20230829"
GENOME_LIST="../1.NCBI.Ref.Seq.Representative.DB/refseq_GCF_eference_genomes.20230829.txt"
NM8013_PROTEOME="Nm_8013_GCF_000026965.1_proteome_nonpseudogenes.faa"
OUTDIR="."

# Extract the name of this protein
GENOME_FILE=$(basename "${GENOME_LIST}")
NUM_GENOMES=$(wc -l ${GENOME_LIST} | awk '{print $1}')
echo " in a total of ${NUM_GENOMES} genomes"
if test -f ${OUTDIR}/${GENOME_FILE}_blast.results.txt
then
	echo "blast already done already procssed...computed"
else 
	COUNTER=1
	echo "Blasting ${j} against ${NUM_GENOMES} genomes...."
	for i in $(cat ${GENOME_LIST})
	do
		GCF=$(echo $i | sed -e 's/\(GCF_[0-9]*\.[0-9]\)_.*/\1/g')
		echo -n "tblastn against genome ${i}, #${COUNTER}/${NUM_GENOMES}..."
	        tblastn -query $NM8013_PROTEOME -num_threads 4 -db ${GENOMES_DBDIR}/${GCF}/cds_from_genomic.fna -qcov_hsp_perc $MIN_ALIGNMENT_LENGTH -seg no -max_hsps 1 -max_target_seqs 5 -evalue 1e-10 -outfmt 6 | sed -e "s/lcl|/${GCF}|lcl|/g" >> ${OUTDIR}/${GENOME_FILE}_blast.results.txt
	done
fi




