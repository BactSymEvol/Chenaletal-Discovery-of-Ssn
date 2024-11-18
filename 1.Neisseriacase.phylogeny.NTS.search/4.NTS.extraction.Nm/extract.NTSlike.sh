#!/bin/bash
### Script to search for NTS-like elements and extract them

NTS="CGTCATTCCCGCG[AC]A[ACG]GCGGGAATC[CT][AG]G"
NTS_LIKE="[CTA]-G-T-C-[AT]-[TC]-N(0,2)-[CT]-[CT]-[CGA]-[CG]-[CTG]-[GA]-[CAT]-N(0,5)-[AT]-N-[GAC]-[CG]-[GCT]-[GA]-G-N(0,2)-[AG]-N-C-[CT]-N-[TCG]"
SEQS=`fuzznuc -sequence ${1} -pattern $NTS_LIKE -complement -pmismatch $2 -pname NTS_like -stdout -warning N -auto | grep "NTS_like:" | awk '{print $6}'`
COUNTER=1
for i in $SEQS
do
	echo ">${1}_NTS_LIKE_No${COUNTER}" 
	echo $i
	COUNTER=$((COUNTER+1))
done
