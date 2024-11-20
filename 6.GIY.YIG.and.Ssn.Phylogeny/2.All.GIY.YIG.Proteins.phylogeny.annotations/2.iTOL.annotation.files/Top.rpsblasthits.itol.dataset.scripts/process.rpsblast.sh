#!/bin/bash
### Master script to call R script to process RPSBLAST output and produce a file for itol
### Author: Luke B Harrison
###
Rscript process.rpsblast.R | awk '{print $1" "$5" "$4}' > refseq_GCF_eference_genomes.20230829.all.giy.yig.updated.lessthan750.c06.all.giy.yig.updated.rpsblasthits.foritol.tsv
