#! /bin/bash

date=$(date '+%Y-%m-%d')

## root to the dir
cd /Users/mariaf/Desktop/GitHub/smORF_EP/smorfep/test


##
## run smorfep for donor and acceptor, forward and reverse strand -- testing examples
smorfep \
-r /Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/ \
-t /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_transcriptCoord_2023-04-25.tsv \
-i /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_introns_2023-04-25.tsv \
-f SMAD2_1-3nts_acceptor.vcf \
-o SMAD2_1-3nts_acceptor_consequences_${date}.vcf

smorfep \
-r /Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/ \
-t /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_transcriptCoord_2023-04-25.tsv \
-i /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_introns_2023-04-25.tsv \
-f SMAD2_1-3nts_donor.vcf \
-o SMAD2_1-3nts_donor_consequences_${date}.vcf

## forward strand
smorfep \
-r /Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/ \
-t /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_transcriptCoord_2023-04-25.tsv \
-i /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_introns_2023-04-25.tsv \
-f CNBD1_1-3nts_acceptor.vcf \
-o CNBD1_1-3nts_acceptor_consequences_${date}.vcf

smorfep \
-r /Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/ \
-t /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_transcriptCoord_2023-04-25.tsv \
-i /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_introns_2023-04-25.tsv \
-f CNBD1_1-3nts_donor.vcf \
-o CNBD1_1-3nts_donor_consequences_${date}.vcf


## 
## Compare the results with VEP annotations
compare2vep \
-s SMAD2_1-3nts_donor_consequences_${date}.vcf \
-v VEP_analysis_SMAD2_donor.vcf \
-o compare_smorf_VEP_SMAD2_1-3nts_donor_${date}.vcf

compare2vep \
-s SMAD2_1-3nts_acceptor_consequences_${date}.vcf \
-v SMAD2_1-3nts_acceptor.vcf \
-o compare_smorf_VEP_SMAD2_1-3nts_acceptor_${date}.vcf


## forward strand
compare2vep \
-s CNBD1_1-3nts_donor_consequences_${date}.vcf \
-v VEP_analysis_CNBD1_donor.vcf \
-o compare_smorf_VEP_CNBD1_1-3nts_donor_${date}.vcf

compare2vep \
-s CNBD1_1-3nts_acceptor_consequences_${date}.vcf \
-v VEP_analysis_CNBD1_acceptor.vcf \
-o compare_smorf_VEP_CNBD1_1-3nts_acceptor_${date}.vcf