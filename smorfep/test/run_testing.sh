#! /bin/bash

date=$(date '+%Y-%m-%d')

echo $date

## forward strand
## CNBD1 donor [no missing variants testing file used]
echo "+ strand donor"
smorfep -r /Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/ -t /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_transcriptCoord_2023-04-25.tsv  -i /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_introns_2023-04-25.tsv -f CNBD1_1-13nts_donor.vcf -o CNBD1_1-13nts_donor_consequences_${date}.tsv > testing_donor${date}.log

compare2vep -s CNBD1_1-13nts_donor_consequences_${date}.tsv -v VEP_analysis_CNBD1_donor_1-13nts.vcf -o compare_CNBD1_donor_1-13nts_VEP_${date}.tsv
echo " "

## CNBD1 acceptor 
echo "+ strand acceptor"
smorfep -r /Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/ -t /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_transcriptCoord_2023-04-25.tsv  -i /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_introns_2023-04-25.tsv -f CNBD1_1-13nts_acceptor.vcf -o CNBD1_1-13nts_acceptor_consequences_${date}.tsv > testing_acceptor${date}.log

compare2vep -s CNBD1_1-13nts_acceptor_consequences_${date}.tsv -v VEP_analysis_CNBD1_acceptor_1-13nts.vcf -o compare_CNBD1_acceptor_1-13nts_VEP_${date}.tsv
echo " "



## reverse strand
## SMAD2 donor
echo "- strand donor"
smorfep -r /Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/ -t /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_transcriptCoord_2023-04-25.tsv  -i /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_introns_2023-04-25.tsv -f SMAD2_1-13nts_donor.vcf -o SMAD2_1-13nts_donor_consequences_${date}.tsv > testing_SMAD2_donor${date}.log

compare2vep -s SMAD2_1-13nts_donor_consequences_${date}.tsv -v VEP_analysis_SMAD2_donor_1-13nts.vcf -o compare_SMAD2_donor_1-13nts_VEP_${date}.tsv
echo " "


## SMAD2 acceptor
echo "- strand acceptor"
smorfep -r /Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/ -t /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_transcriptCoord_2023-04-25.tsv  -i /Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation_introns_2023-04-25.tsv -f SMAD2_1-13nts_acceptor.vcf -o SMAD2_1-13nts_acceptor_consequences_${date}.tsv > testing_SMAD2_acceptor${date}.log

compare2vep -s SMAD2_1-13nts_acceptor_consequences_${date}.tsv -v VEP_analysis_SMAD2_acceptor_1-13nts.vcf -o compare_SMAD2_acceptor_1-13nts_VEP_${date}.tsv
echo " "

echo "DONE" 
