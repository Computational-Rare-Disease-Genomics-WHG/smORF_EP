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

## NOTE: compute the 5'utr, cds, and 3'utr is not possible here as the transcript ID is not fetched yet
## we can't use the gene ID as is shared sometimes by multiple transcripts

## split last column and keep only ID, gene_id, transcript_id and transcrypt_type
## ensmbl entris and havana entries have different format 

## As we want just: ID, gene_id, transcript_id and transcript_type
## we remove everything after, so the last column is the same.
gencode['new_info'] = gencode['info'].str.split(';transcript_name').str[0] ## 0 as we keep the first half
## remove info column -- repeated info
gencode = gencode.drop(columns='info')

gencode[['ID', 'Parent', 'gene_id', 'transcript_id', 'gene_type' , 'gene_name', 'transcript_type']] = gencode['new_info'].str.split(';',expand=True)

gencode = gencode.drop(columns=['Parent','gene_name']) ## don't drop transcipt_id as cds and other elements have their typre prefix
## NOTE: transcript_id is removed as we are saving per transcript, so the ID column is the same as transcript_id
## NOTE2: We keep gene type now as there is differnce in CDS and exon, if case of pesudogene there are only exon annotations

## remove info column -- repeated info
gencode = gencode.drop(columns='new_info')

## remove prefix on the split columns 
for i in ['ID', 'gene_id', 'gene_type', 'transcript_id', 'transcript_type']: 
    gencode[i] = gencode[i].str.split('=').str[1]

## 2 - transcripts save transcripts coordinates
gencode_transcripts = gencode[gencode.type == 'transcript']
##print(gencode_transcripts.shape)
## total transcripts: 251 236
## transcripts no MT: 251 199


##print(gencode_transcripts)
print(gencode_transcripts.shape[0], 'transcripts')

## NEW BLOCK - March 2024
##print("new block")

## add new columns for the coordinates of CDS/exon, 5'UTR and 3'UTR
gencode_transcripts.loc[:, ['CDS/exon', 'five_prime', 'three_prime']] = 'ND'

## Keep only CDS, exon, 5'UTR and 3'UTR annotations
gencode_filtered = gencode[gencode['type'].isin(['CDS','five_prime_UTR', 'three_prime_UTR'])]
##print(gencode_filtered['type'].unique())
##print(gencode_filtered)

grouped_df = gencode_filtered.groupby('transcript_id')

# Define the types you want to check for
types_to_check = {'CDS', 'five_prime_UTR', 'three_prime_UTR'}

trasncripts_processed = 1
# Filter groups that have all three types
groups_with_all_types = []
for group_name, group_data in grouped_df:
    if set(group_data['type']) == types_to_check:
        groups_with_all_types.append(group_name)
        ## make the computations directly
        fiveprime_coord, cds_coord, threeprime_coord = compute_start_end_coordinate(group_data)
        ##print(cds_coord,fiveprime_coord,threeprime_coord, transcript_df['transcript_type'])

        rowIndex = int(gencode_transcripts.index[gencode_transcripts['ID'] == group_name].tolist()[0])
        gencode_transcripts.at[rowIndex, 'CDS/exon'] = cds_coord
        gencode_transcripts.at[rowIndex, 'five_prime'] = fiveprime_coord
        gencode_transcripts.at[rowIndex, 'three_prime'] = threeprime_coord
        ##print(gencode_transcripts.loc[rowIndex])
    
    if trasncripts_processed%10000 ==0:
        print(trasncripts_processed, 'transcripts processed')
    trasncripts_processed+= 1

# Print the groups that have all three types
print("Groups with all three types:", len(groups_with_all_types))



# trasncripts_processed = 1
# for each_id_t in transcript_ids_list:

#     ## check if the transcript is a lnRNA - those do not have 5'UTr, CDS, 3'UTR annotations
#     transc_type = gencode_new[gencode_new['transcript_id'] == each_id_t]['transcript_type'].iloc[0]

#     transcript_df = gencode_new[gencode_new['transcript_id'] == each_id_t] ## "ENST00000456328.2"]
#     ##print(transcript_df)
#     ##print(transcript_df['type'].unique())
#     regions_transc = transcript_df['type'].unique()
#     ## NOTE: To check, next line is case sensitive 
#     if 'three_prime_UTR' in regions_transc and 'five_prime_UTR' in regions_transc: 

#         fiveprime_coord, cds_coord, threeprime_coord = compute_start_end_coordinate(transcript_df)
#         ##print(cds_coord,fiveprime_coord,threeprime_coord, transcript_df['transcript_type'])

#         ## only adds the coordinates to the dataframe when all the 3 regions are defined
#         if cds_coord != None and fiveprime_coord != None and threeprime_coord != None: 
#             ##print('all 3 defined')
#             ## add coordinates to the specific line
#             rowIndex = int(gencode_transcripts.index[gencode_transcripts['ID'] == each_id_t].tolist()[0])
#             gencode_transcripts.at[rowIndex, 'CDS/exon'] = cds_coord
#             gencode_transcripts.at[rowIndex, 'five_prime'] = fiveprime_coord
#             gencode_transcripts.at[rowIndex, 'three_prime'] = threeprime_coord
#             ##print(gencode_transcripts.loc[rowIndex])
            



## write transcripts to file
gencode_transcripts.to_csv(output_transcript_coord, sep='\t', lineterminator='\n', index=False)

## DONE!


end_time = time.time() - start_time

print(end_time, ' seconds')

print('DONE')
