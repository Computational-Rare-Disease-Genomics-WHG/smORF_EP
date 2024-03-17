#!/usr/bin/python3

## Author: M Fernandes

##  collect the information from the GENCODE file: 
## transcript coordinates -- start and end


import sys
import time
import pandas as pd

from smorfep.utils.functions import read_file
from smorfep.utils.functions import compute_start_end_coordinate

##print('python3 compute_transcripts_GENCODE.py <inputfilemane> <outputfilename_transcript_coodrinates>')
##print('Example: python3 compute_transcripts_GENCODE.py gencode/gencode.v41.annotation_columnNames.gff3 gencode/gencode.v41.annotation_transcriptCoord.tsv')

start_time = time.time()

gencode_file = sys.argv[1]
output_transcript_coord = sys.argv[2]

gencode = read_file(gencode_file, '\t', 0)
# print(gencode)
# print(gencode.keys())
##print(gencode.shape) ## total lines: 3 373 604
##print(gencode.shape[0], ' transcripts')

##How may different categories are in the file
#print(gencode['type'].unique())

##gencode_g = gencode[gencode.type == 'gene']
##print(gencode_g.shape)
## total genes: 61 852

## 1- remove MT annotations
gencode = gencode[gencode.chr != 'chrM']
## find substring in the column 
#check = gencode.loc[gencode['info'].str.contains("ENST00000456328.2",case=False)]

## they have different formating than others and we are not able to split in the same way
gencode = gencode[gencode.type != 'gene']

## 2 - transcripts save transcripts coordinates
gencode_transcripts = gencode[gencode.type == 'transcript']
##print(gencode_transcripts.shape)
## total transcripts: 251 236
## transcripts no MT: 251 199

## NOTE: compute the 5'utr, cds, and 3'utr is not possible here as the transcript ID is not fetched yet
## we can't use the gene ID as is shared sometimes by multiple transcripts

## split last column and keep only ID, gene_id, transcript_id and transcrypt_type
## ensmbl entris and havana entries have different format 

## As we want just: ID, gene_id, transcript_id and transcript_type
## we remove everything after, so the last column is the same.
gencode_transcripts['new_info'] = gencode_transcripts['info'].str.split(';transcript_name').str[0] ## 0 as we keep the first half
## remove info column -- repeated info
gencode_transcripts = gencode_transcripts.drop(columns='info')

gencode_transcripts[['ID', 'Parent', 'gene_id', 'transcript_id', 'gene_type' , 'gene_name', 'transcript_type']] = gencode_transcripts['new_info'].str.split(';',expand=True)

gencode_transcripts = gencode_transcripts.drop(columns=['Parent', 'transcript_id','gene_name'])
## NOTE: transcript_id is removed as we are saving per transcript, so the ID column is the same as transcript_id
## NOTE2: We keep gene type now as there is differnce in CDS and exon, if case of pesudogene there are only exon annotations

## remove info column -- repeated info
gencode_transcripts = gencode_transcripts.drop(columns='new_info')

## remove prefix on the split columns 
for i in ['ID', 'gene_id', 'gene_type', 'transcript_type']: 
    gencode_transcripts[i] = gencode_transcripts[i].str.split('=').str[1]

##print(gencode_transcripts)
print(gencode_transcripts.shape[0], 'transcripts')

gene_ids_list = gencode_transcripts['gene_id'].unique()
print("gene_ids:", len(gencode_transcripts['gene_id'].unique()))
print("much smaller than transcript ids number")



## TODO: Add 5'UTR, CDS and 3'UTR columns to transcript df - default value = 0
print("new block")


## TODO: Editing here
# print(gencode["type"].unique())
# test_5utr = gencode[gencode.type == 'five_prime_UTR']
# for index, row in test_5utr.iterrows():
#     print(row.info)

## Keep only CDS, exon, 5'UTR and 3'UTR annotations
gencode_new = gencode[gencode.type != 'transcript']
gencode_new = gencode_new[gencode_new.type != 'start_codon']
gencode_new = gencode_new[gencode_new.type != 'stop_codon']
gencode_new = gencode_new[gencode_new.type != 'stop_codon_redefined_as_selenocysteine']
##print(gencode_new['type'].unique())

## TODO: format info 
## filter to transcript 
## run compute start and end function
gencode_new['new_info'] = gencode_new['info'].str.split(';transcript_name').str[0] ## 0 as we keep the first half
## remove info column -- repeated info
gencode_new = gencode_new.drop(columns='info')

gencode_new[['ID', 'Parent', 'gene_id', 'transcript_id', 'gene_type' , 'gene_name', 'transcript_type']] = gencode_new['new_info'].str.split(';',expand=True)

gencode_new = gencode_new.drop(columns=['Parent', 'gene_type', 'gene_name','gene_id', 'gene_type']) ## here we need to keep transcript_id as the ID will be different for non-"transcipt" annotations
## remove info column -- repeated info
gencode_new = gencode_new.drop(columns='new_info')
##print(gencode_new.shape) ## Not empty

## remove prefix on the split columns 
for i in ['ID', 'transcript_id', 'transcript_type']: ## we match per transcript ID and the remaining information is already in the gencode_transcriptd DF
    gencode_new[i] = gencode_new[i].str.split('=').str[1]

## collect the unique ids for the transcripts
transcript_ids_list = list(gencode_transcripts['ID'])
print(transcript_ids_list)

##for each_id in transcript_ids_list
transcript_df = gencode_new[gencode_new['transcript_id'] == "ENST00000379410.8"] ## "ENST00000456328.2"]

cds_df, fiveprime_df, threeprime_df = compute_start_end_coordinate(transcript_df)

print(gencode_transcripts.columns)
gencode_transcripts['CDS_exon'] = 'ND'
gencode_transcripts['five_prime'] = 'ND'
gencode_transcripts['three_prime'] = 'ND'
print(gencode_transcripts.columns)

## write transcripts to file
gencode_transcripts.to_csv(output_transcript_coord, sep='\t', lineterminator='\n', index=False)

## DONE!


end_time = time.time() - start_time

print(end_time, ' seconds')

print('DONE')
