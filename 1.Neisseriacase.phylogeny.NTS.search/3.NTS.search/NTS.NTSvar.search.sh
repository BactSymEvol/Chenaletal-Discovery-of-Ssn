#!/bin/bash
## Author: Luke B Harrison (luke.harrison@mail.mcgill.ca/luke.harrison@inrs.ca)
## Script to tabulate the number of NTS, NTSvar and DUS sequences in a given set of genomes
## Depends: EmBOSS tools / fuzznuc installed in an executable path
## Inputs = root name of genomes to iterate over (not including .fna / .fasta suffix) given on command line
## 	    Can be helpful to pass a control file (1 / line, with `cat genomelist.txt`
## Output = NTS.NTSvar.motif.search.csv

GENOME_DIRECTORY="../1.Neisseriacae.genomes/"
FASTA_SUFFIX="fna"

NTSSEQ="CGTCATTCCCGCG[AC]A[ACG]GCGGGAATC[CT][AG]G"
NTSVARSEQ="CGTCATACTCGGGNNNNNCCCGAGTATC"

echo "GENOME,DUS,NTS,NTSvar" > NTS.NTSvar.motif.search.csv

for i in $*
do
	echo "Processing ${i}..."
	echo -n "${i}," >> NTS.NTSvar.motif.search.csv
	## DUS
	AT_DUS=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern ATGCCGTCTGAA -complement -pmismatch 0 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	TG_wadDUS=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern TGCCTGTCTGAA -complement -pmismatch 0 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	AG_DUS=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern AGGCCGTCTGAA -complement -pmismatch 0 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	AG_musDUS=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern AGGTCGTCTGAA -complement -pmismatch 0 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	AG_simDUS=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern AGGCTGCCTGAA -complement -pmismatch 0 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	AG_kingDUS=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern AGGCAGCCTGAA -complement -pmismatch 0 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	AG_king3DUS=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern AAGCAGCCTGCA -complement -pmismatch 0 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	AG_eikDUS=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern AGGCTACCTGAA -complement -pmismatch 0 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	DUS_tot=$(($AT_DUS + $TG_wadDUS + $AG_DUS +$AG_musDUS +$AG_simDUS +$AG_kingDUS + $AG_king3DUS + $AG_eikDUS))
	echo -n "${DUS_tot}," >> NTS.NTSvar.motif.search.csv

	## NTS
	NTS=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern $NTSSEQ -complement -pmismatch 1 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	echo -n "${NTS}," >> NTS.NTSvar.motif.search.csv

	## NTSvar
	NTSvar=`fuzznuc -sequence ${GENOME_DIRECTORY}/${i}.$FASTA_SUFFIX -pattern $NTSVARSEQ -complement -pmismatch 1 -stdout -auto|grep "Reported_hitcount" | awk '{print $3}'`
	echo "${NTSvar}" >> NTS.NTSvar.motif.search.csv
done
