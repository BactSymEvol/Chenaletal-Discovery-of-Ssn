#!/bin/bash
emapper.py --cpu 12 -i Nm8013_proteome_cdnas_flankedby.NTS.translated.fasta --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 -o Nm8013_proteome_cdnas_flankedby.NTS.translated
