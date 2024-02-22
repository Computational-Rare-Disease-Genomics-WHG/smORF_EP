#!/bin/bash

## script to generate the report per smORF -- variants we want to look up 

date_today=$(date '+%Y-%m-%d')
##echo ${date_today}

create_report variants.vcf \ ##SMAD2_noheader.bed \
--fasta ../../../../ref_genome/hg38.fa \
--flanking 1000 \
--info-columns PRED_EFFECT SMORF CLNSIG CLNREVSTAT \
--tracks SMAD2_noheader.bed variants_pathogenic.vcf variants_uncertain.vcf variants_benign.vcf overlap.bed \
--output test_example_genome_${date_today}.html

##--fasta https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg38/hg38.fa \
