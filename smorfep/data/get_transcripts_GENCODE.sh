#!/bin/bash

## Script to obtain transcript information from GENCODE
## Adapt to other sources

## Uses wget

now="$(date +'%Y-%m-%d')"
##echo ${now}

## 1- Creates a new dir for the reference genome
mkdir gencode

## 2- Download GENCODE
wget -P gencode/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gff3.gz

## 3- Uncompress reference genome and delete compressed version (to keep the compressed file too add -k option )
gzip -d gencode/gencode.v41.annotation.gff3.gz

## 4- Precomputations on genecode file 
## Pre-process gff3 file  - single header
python3 preProcess_gff.py gencode/gencode.v41.annotation.gff3 gencode/gencode.v41.annotation_columnNames.gff3

## transcript coordinates (15.seconds M1 ship)
python3 compute_transcripts_GENCODE.py gencode/gencode.v41.annotation_columnNames.gff3 gencode/gencode.v41.annotation_transcriptCoord_${now}.tsv
## intron coordinates
python3 compute_introns_GENCODE_perTransc.py gencode/gencode.v41.annotation_columnNames.gff3 gencode/gencode.v41.annotation_introns_${now}.tsv