#!/bin/bash

## Human reference genome download and pre-processing
## uses wget and gzip

## 1- Creates a new dir for the reference genome
mkdir ref_genome

## 2- Download reference from NCBI to the ref_genome/repository -- For GRC38.p13
wget -P ref_genome/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109.20211119/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz

## 3- Uncompress reference genome and delete compressed version (to keep the compressed file too add -k option )
gzip -d ref_genome/GCF_000001405.39_GRCh38.p13_genomic.fna.gz

## 4- Split reference per chromosome
python3 ref_per_chr.py GCF_000001405.39_GRCh38.p13_genomic.fna _GRCh38.p13_genomic.fna
## Excludes 'scaffold', 'patch', 'mitochondrion' sequences 
## If others to exclude add condition to line 24 
## outputname format: chr<num>_<suffix_user_given> (line above: chr<num>_GRCh38.p13_genomic.fna)