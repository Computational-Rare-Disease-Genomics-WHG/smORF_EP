#!/usr/bin/python3

## Author: M Fernandes

##  collect the information from the GENCODE file: 
## transcript coordinates -- start and end


import sys
import time
import pandas as pd

from smorfep.utils.functions import read_file

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

## 2 - transcripts save transcripts coordinates
gencode_transcripts = gencode[gencode.type == 'transcript']
##print(gencode_transcripts.shape)
## total transcripts: 251 236
## transcripts no MT: 251 199


## split last column and keep only ID, gene_id, transcript_id and transcrypt_type
## ensmbl entris and havana entries have different format 

## As we want just: ID, gene_id, transcript_id and transcript_type
## we remove everything after, so the last column is the same.
gencode_transcripts['new_info'] = gencode_transcripts['info'].str.split(';transcript_name').str[0] ## 0 as we keep the first half
## remove info column -- repeated info
gencode_transcripts = gencode_transcripts.drop(columns='info')

gencode_transcripts[['ID', 'Parent', 'gene_id', 'transcript_id', 'gene_type' , 'gene_name', 'transcript_type']] = gencode_transcripts['new_info'].str.split(';',expand=True)
gencode_transcripts = gencode_transcripts.drop(columns=['Parent', 'gene_type' , 'gene_name'])
## remove info column -- repeated info
gencode_transcripts = gencode_transcripts.drop(columns='new_info')

## remove prefix on the split columns 
for i in ['ID', 'gene_id', 'transcript_id', 'transcript_type',]:
    gencode_transcripts[i] = gencode_transcripts[i].str.split('=').str[1]

##print(gencode_transcripts)
print(gencode_transcripts.shape[0], 'transcripts')


## write transcripts to file
gencode_transcripts.to_csv(output_transcript_coord, sep='\t', lineterminator='\n', index=False)

## DONE!


end_time = time.time() - start_time

print(end_time, ' seconds')

print('DONE')
