#!/bin/bash
### Script to search for NTS-like elements in the entire ref-seq respresentative genomes
### Input: $1 a list of NCBI genome assembly IDs corresponding to teh DB in question
### This script will check if there is a homolog of SSNa present
###

## PARAMETERS
HOMOLOGY_CUTOFF=60
MIN_ALIGNMENT_LENGTH=75

GENOMES_DBDIR="../3.Nm8013.Protein.Homologs.Flanked.byNTS.in.NCBI/1.NCBI.Ref.Seq.Representative.DB/RefSeq_genomes_20230829"
CDD_DB="cdd.database/Cdd"

SSNA="SsnA.Nm.8013.faa"
CDD_PSSM="cd10448.rescaled.smp"

NTSSEQ="CGTCATTCCCGCG[AC]A[ACG]GCGGGAATC[CT][AG]G"
NTS_LIKE="[CTA]-G-T-C-[AT]-[TC]-N(0,2)-[CT]-[CT]-[CGA]-[CG]-[CTG]-[GA]-[CAT]-N(0,5)-[AT]-N-[GAC]-[CG]-[GCT]-[GA]-G-N(0,2)-[AG]-N-C-[CT]-N-[TCG]"

OUTPUT_FILE="SsnA_NTS_homologs.PSI.BLAST.refseq_reference_genomes_flanks.csv"
OUTDIR="Regions.flanking.SSna.Homologs"

# scratch space / temporary file space
SLURM_TMPDIR="/tmp"


echo "GENOME,NTS,NTS_LIKE,HOMOLOG,HOMOLOG${HOMOLOGY_CUTOFF},HOMOLOG_PSI,HOMOLOG_CDD" > $OUTPUT_FILE

FLANK_SIZE=5000

COUNTER=1
NUM_GENOMES=`wc -l $1| awk '{print $1}'`

## Loop over genomes
for i in `cat ../3.Nm8013.Protein.Homologs.Flanked.byNTS.in.NCBI/1.NCBI.Ref.Seq.Representative.DB/refseq_GCF_ref_rep_genomes.20230829.txt`
do
	echo -n "Processing ${i}..... , #${COUNTER}/${NUM_GENOMES}..."
	GCF=`echo ${i} | sed -e 's/\(GCF_[0-9]*\.[0-9]\)_.*/\1/g'`
	echo -n "${i}," >> $OUTPUT_FILE
	## Search for NTS elements
	echo -n "`fuzznuc -sequence ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna -pattern $NTSSEQ -complement -pmismatch 1 -pname NTS -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`," >> $OUTPUT_FILE
	## Search for NTS-like elements
	echo -n "`fuzznuc -sequence ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna -pattern $NTS_LIKE -complement -pmismatch 0 -pname NTSLIKE -stdout -warning N -auto|grep "Reported_hitcount" | awk '{print $3}'`," >> $OUTPUT_FILE 
	## Now look for a homolog at our current level of cutoff
        #BLAST_OUT=`tblastn -query $SSNA -db ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna -qcov_hsp_perc $MIN_ALIGNMENT_LENGTH -seg no -max_target_seqs 1 -max_hsps 1 -evalue 1e-10 -outfmt 6 2> /dev/null`
	BLAST_OUT=`tblastn -query $SSNA -db ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna -qcov_hsp_perc $MIN_ALIGNMENT_LENGTH -seg no -max_target_seqs 5 -max_hsps 1 -evalue 1e-10 -outfmt 6 2> /dev/null| head -n1`
        if [ -n "$BLAST_OUT" ]
        then
                echo -n "YES," >> $OUTPUT_FILE
                PER_IDENT=`echo $BLAST_OUT | awk '{print $3}'`
		### Now extract flanking region from the corresponding
		echo $BLAST_OUT | cut -d " " -f 2,9,10 | awk '{if($2 > $3) {print $1,$3,$2,"rev",1,"-";} else if ($2 < $3) {print $1,$2,$3,"fwd",1,"+";}}' OFS='\t' | awk '{a=$2-1; $2 = a; print;}' OFS='\t' > ${SLURM_TMPDIR}/${i}_tblastn.bed
		bedtools flank -i ${SLURM_TMPDIR}/${i}_tblastn.bed -g ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna.fai -b ${FLANK_SIZE} > ${SLURM_TMPDIR}/${i}_tblastn_${FLANK_SIZE}.bed
		bedtools getfasta -s -fi ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna -bed ${SLURM_TMPDIR}/${i}_tblastn_${FLANK_SIZE}.bed > ${SLURM_TMPDIR}/${i}_tblastn_${FLANK_SIZE}.fasta
		HEADER=`grep ">" ${SLURM_TMPDIR}/${i}_tblastn_${FLANK_SIZE}.fasta | tr '\n' '|' | sed -e 's/|>/|/g' -e 's/|$//g'`
		echo ">${HEADER}" >> ${OUTDIR}/Flanks_tblastn_${FLANK_SIZE}.fasta
		grep -v ">" ${SLURM_TMPDIR}/${i}_tblastn_${FLANK_SIZE}.fasta | tr -d '\n' | sed -e 's/$/\n/g'>> ${OUTDIR}/Flanks_tblastn_${FLANK_SIZE}.fasta
		rm ${SLURM_TMPDIR}/${i}_tblastn.bed 
                if awk "BEGIN {exit !($PER_IDENT >= $HOMOLOGY_CUTOFF)}"
		then
			echo -n "YES," >> $OUTPUT_FILE
			## Here we have a matching homolog with good identity match
	                #bedtools getfasta -s -fi ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna -bed ${SLURM_TMPDIR}/${i}_tblastn_${FLANK_SIZE}.bed 
			echo ">${HEADER}" >> ${OUTDIR}/Flanks_tblastn_${HOMOLOGY_CUTOFF}_${FLANK_SIZE}.fasta
	 		grep -v ">" ${SLURM_TMPDIR}/${i}_tblastn_${FLANK_SIZE}.fasta | tr -d '\n' | sed -e 's/$/\n/g' >> ${OUTDIR}/Flanks_tblastn_${HOMOLOGY_CUTOFF}_${FLANK_SIZE}.fasta

		else 
			echo -n "NO," >> $OUTPUT_FILE
		fi
		rm ${SLURM_TMPDIR}/${i}_tblastn_${FLANK_SIZE}.bed ${SLURM_TMPDIR}/${i}_tblastn_${FLANK_SIZE}.fasta
	else
		echo -n "NO,NO," >> $OUTPUT_FILE
        fi
	PSI_BLAST_OUT=`psiblast -in_pssm ${CDD_PSSM} -db ${GENOMES_DBDIR}/${GCF}/protein.faa -seg no -max_target_seqs 5 -max_hsps 1 -evalue 1e-5 -outfmt 6 | head -n1`
	#PSI_BLAST_OUT=`psiblast -in_pssm ${CDD_PSSM} -db ${GENOMES_DBDIR}/${GCF}/protein.faa -qcov_hsp_perc $MIN_ALIGNMENT_LENGTH -seg no -max_target_seqs 5 -max_hsps 1 -evalue 1e-2 outfmt 6  2> /dev/null| head -n1`
	if [ -n "$PSI_BLAST_OUT" ]
        then
		### Ensure the protein is not too long
		PROT=`echo $PSI_BLAST_OUT | cut -d " " -f 2`
		LEN=`grep ${PROT} ${GENOMES_DBDIR}/${GCF}/protein.faa.fai | cut -f 2`
                echo "Lenght of matched prot = $LEN"
		if awk "BEGIN {exit !($LEN >= 165)}"
                then
			echo -n "NO," >> $OUTPUT_FILE
		else 
			echo -n "YES," >> $OUTPUT_FILE
		fi
		### Now script to reverse domain blast using NCBI-supplied cutoff to confirm membership inthe protein family
		fastafetch -f ${GENOMES_DBDIR}/${GCF}/protein.faa -i ${GENOMES_DBDIR}/${GCF}/protein.faa.faidx -q ${PROT} > ${SLURM_TMPDIR}/${PROT}.faa

		### Now pull out the extracted flanks
		RE_PSI_BLAST_OUT=`tblastn -query ${SLURM_TMPDIR}/${PROT}.faa -db ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna -qcov_hsp_perc 90 -seg no -max_target_seqs 5 -max_hsps 1 -evalue 1e-10 -outfmt 6 | head -n1`
		if [ -n "$RE_PSI_BLAST_OUT" ]
                then
			echo $RE_PSI_BLAST_OUT | cut -d " " -f 2,9,10 | awk '{if($2 > $3) {print $1,$3,$2,"rev",1,"-";} else if ($2 < $3) {print $1,$2,$3,"fwd",1,"+";}}' OFS='\t' | awk '{a=$2-1; $2 = a; print;}' OFS='\t' > ${SLURM_TMPDIR}/${i}_psiblast.bed
        	        bedtools flank -i ${SLURM_TMPDIR}/${i}_psiblast.bed -g ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna.fai -b ${FLANK_SIZE} > ${SLURM_TMPDIR}/${i}_psiblast_${FLANK_SIZE}.bed
                	bedtools getfasta -s -fi ${GENOMES_DBDIR}/${GCF}/${i}_genomic.fna -bed ${SLURM_TMPDIR}/${i}_psiblast_${FLANK_SIZE}.bed > ${SLURM_TMPDIR}/${i}_psiblast_${FLANK_SIZE}.fasta
			PSI_HEADER=`grep ">" ${SLURM_TMPDIR}/${i}_psiblast_${FLANK_SIZE}.fasta | tr '\n' '|' | sed -e 's/|>/|/g' -e 's/|$//g'`
			echo ">${PSI_HEADER}" >> ${OUTDIR}/Flanks_psiblast_${FLANK_SIZE}.fasta
	                grep -v ">" ${SLURM_TMPDIR}/${i}_psiblast_${FLANK_SIZE}.fasta | tr -d '\n' | sed -e 's/$/\n/g' >> ${OUTDIR}/Flanks_psiblast_${FLANK_SIZE}.fasta
		else
			echo "ERROR: couldn't find the psi-blast match in teh actual genome"
		fi
		RPS_BLAST_OUT=`rpsblast -query ${SLURM_TMPDIR}/${PROT}.faa -db $CDD_DB -seg no -evalue 10 -outfmt 6 | grep "CDD:198395"`
		if [ -n "$RPS_BLAST_OUT" ]
		then
			## Check bitscore
			RPSBLAST_SCORE=`echo $RPS_BLAST_OUT | cut -d " " -f 12`
			if awk "BEGIN {exit !($RPSBLAST_SCORE >= 100)}"
			then
				## Yes, we have found an appropriate member of our family
				echo "YES" >> $OUTPUT_FILE
				echo ">${PSI_HEADER}" >> ${OUTDIR}/Flanks_rpsblast_${FLANK_SIZE}
        	                grep -v ">" ${SLURM_TMPDIR}/${i}_psiblast_${FLANK_SIZE}.fasta |tr -d '\n' | sed -e 's/$/\n/g' >> ${OUTDIR}/Flanks_rpsblast_${FLANK_SIZE}
			else
				echo "We had a hit to domain CDD:198395 / cd10448 but the bitscore was lower than the threshold ${RPSBLAST_SCORE} vs 104.885"
				echo "NO" >> $OUTPUT_FILE
			fi
		else
			echo "We had a PSI-Blast hit but no match to the domain DB"
			echo "NO" >> $OUTPUT_FILE
		fi
	       	rm ${SLURM_TMPDIR}/${i}_psiblast_${FLANK_SIZE}.bed ${SLURM_TMPDIR}/${i}_psiblast_${FLANK_SIZE}.fasta
               	rm ${SLURM_TMPDIR}/${i}_psiblast.bed
		rm ${SLURM_TMPDIR}/${PROT}.faa
	else
		echo "NO,NO" >> $OUTPUT_FILE
	fi
        ((COUNTER=COUNTER+1))
	echo "done."
done


