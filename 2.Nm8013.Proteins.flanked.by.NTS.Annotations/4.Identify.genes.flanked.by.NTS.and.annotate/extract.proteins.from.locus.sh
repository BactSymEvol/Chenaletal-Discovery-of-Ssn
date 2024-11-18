#!/bin/bash

### Author: Luke B Harrison
## Script to extract and translate CDNS that are flanked by NTS regions in the Nm8013 proteome

CDNA_LIST="Nm8013_proteome_cdnas_flankedby.NTS.csv"
CDS="../2.NM8013.Genome/GCA_000026965.1_ASM2696v1_cds_from_genomic.fna"
seqkit grep --use-regexp --by-name --pattern-file ${CDNA_LIST} $CDS > Nm8013_proteome_cdnas_flankedby.NTS.fasta
fastatranslate -f Nm8013_proteome_cdnas_flankedby.NTS.fasta -F 1 > Nm8013_proteome_cdnas_flankedby.NTS.translated.fasta

