#!/usr/bin/python3

## Author: M Fernandes


##  collect the information from the GENCODE file: 
## introns locations -- pertranscript


import sys
import time
import pandas as pd

from smorfep.utils.functions import read_file

##print('python3 compute_introns_GENCODE.py <inputfilemane> <outputfilename_introns>')
##print('Example: python3 compute_introns_GENCODE.py gencode/gencode.v41.annotation_columnNames.gff3 gencode/gencode.v41.annotation_introns.tsv')

start_time = time.time()

gencode_file = sys.argv[1]
output_introns = sys.argv[2]

gencode = read_file(gencode_file, '\t', 0)

##print(gencode.shape) ## total lines: 3 373 604

## 1- remove MT annotations
gencode = gencode[gencode.chr != 'chrM']

## 2 - for these computation we only need the exons - EXONS ONLY 
## differently of MANE, there is 'exon' annotation on the 5'UTR and 3'UTR
## ex. if the 5'UTR and CDS do not have a gap -- there is a 'exon' entry from 5'UTR start until the CDS end
gencode = gencode[gencode.type == 'exon']


## split last column and keep only ID, gene_id, transcript_id and transcrypt_type
## ensmbl entris and havana entries have different format 

## As we want just: ID, gene_id, transcript_id and transcript_type
## we remove everything after, so the last column is the same.
gencode['new_info'] = gencode['info'].str.split(';exon_id').str[0] ## 0 as we keep the first half
## remove info column -- repeated info
gencode = gencode.drop(columns='info')

gencode[['ID', 'Parent', 'gene_id', 'transcript_id', 'gene_type' , 'gene_name', 'transcript_type', 'transcript_name', 'exon_number']] = gencode['new_info'].str.split(';',expand=True)
gencode = gencode.drop(columns=['Parent', 'gene_type' , 'gene_name', 'transcript_name'])
## remove info column -- repeated info
gencode = gencode.drop(columns='new_info')

## remove prefix on the split columns 
for i in ['ID', 'gene_id', 'transcript_id', 'transcript_type', 'exon_number',]:
    gencode[i] = gencode[i].str.split('=').str[1]

##print(gencode)
#print(gencode.shape)

## all unique transcript IDS
all_transcripts = gencode['transcript_id'].unique()
print('unique transcript IDs ', len(all_transcripts))

##print(gencode['strand'].unique())

all_chromosomes = gencode['chr'].unique()
##print(all_chromosomes)

print('Computing introns...')

total_transcript_introns = 0
progress = 1

header = True ## to write the header only the first time in the output

try: 
    transcripts_with_introns = [] ## unique 
    for c in all_chromosomes:
        transcripts_chr = gencode[gencode.chr == c] ## transcripts per chromosome
        transcripts_ids_chr = transcripts_chr['transcript_id'].unique()
        for each in transcripts_ids_chr:
            ##print(each)

            ## new dataframe per chromosome as we write the output per transcript
            introns_df = pd.DataFrame(data=None, columns=gencode.columns)   ## empty dataframe with the same columns as MANE
            ## rename last column into intron_number
            introns_df = introns_df.rename({'exon_number': 'intron_number'}, axis='columns')


            num_exons = transcripts_chr[transcripts_chr.transcript_id == each].shape[0]
            if num_exons > 1: ## number of rows/exons ## if there are introns
                transc = transcripts_chr[transcripts_chr.transcript_id == each]
                ##print(transc)

                strand = transc['strand'].iloc[0]
                ##print(strand)

                introns = num_exons - 1
                #print(introns)
                i = 1

                if strand == '+': 
                    first_exon = int(transc.iloc[0]['exon_number']) ## takes the first line as in gencode the exons are ordered
                    while i <= introns:

                        index_1 = transc.index[transc['exon_number'] == str(first_exon)].tolist() ## exons_number is char
                        ##print('index ', index_1)
                        intron_start = transc.loc[index_1,'end'].item() + 1 # position after exon ends

                        index_2 = transc.index[transc['exon_number'] == str(first_exon+1)].tolist()
                        intron_end = transc.loc[index_2, 'start'].item() - 1  ## position before next exon starts

                        ##print(int(intron_start), int(intron_end))
                        #print(transc.loc[index_1,'chr'].item())
                        
                        intron_computed = pd.DataFrame(
                            {'chr': [transc.loc[index_1,'chr'].item()],
                            'source': [transc.loc[index_1,'source'].item()],
                            'type': ['intron'],
                            'start': [int(intron_start)],
                            'end': [int(intron_end)],
                            'score': [transc.loc[index_1,'score'].item()],
                            'strand': [transc.loc[index_1,'strand'].item()],
                            'gen_phase': [transc.loc[index_1,'gen_phase'].item()],
                            'ID': [transc.loc[index_1,'ID'].item()], 
                            'gene_id': [transc.loc[index_1,'gene_id'].item()],
                            'transcript_id': [transc.loc[index_1,'transcript_id'].item()], 
                            'transcript_type': [transc.loc[index_1,'transcript_type'].item()], 
                            'intron_number': [int(i)]
                            })

                        introns_df = pd.concat([introns_df, intron_computed])

                        i += 1 
                        first_exon +=1

                        # print(introns_df)
                        # sys.exit()
                    
                    ##print(introns_df)
                    ##sys.exit()


                elif strand == '-':
                    ##print('strand -')
                    first_exon = int(transc.iloc[0]['exon_number'])
                    while i <= introns:

                        ## index1 for exon 1
                        index_1 = transc.index[transc['exon_number'] == str(first_exon)].tolist() ## exons_number is char
                        index_2 = transc.index[transc['exon_number'] == str(first_exon+1)].tolist() ## exon2
                        ##print('index_1 ', index_1)
                        ##print('index_2 ', index_2)
                        intron_start = transc.loc[index_2,'end'].item() + 1 # position after exon ends
                        intron_end = transc.loc[index_1, 'start'].item() - 1  ## position before next exon starts
                        ##print(intron_start, intron_end)

                        intron_computed = pd.DataFrame(
                            {'chr': [transc.loc[index_1,'chr'].item()],
                            'source': [transc.loc[index_1,'source'].item()],
                            'type': ['intron'],
                            'start': [int(intron_start)],
                            'end': [int(intron_end)],
                            'score': [transc.loc[index_1,'score'].item()],
                            'strand': [transc.loc[index_1,'strand'].item()],
                            'gen_phase': [transc.loc[index_1,'gen_phase'].item()],
                            'ID': [transc.loc[index_1,'ID'].item()],  
                            'gene_id': [transc.loc[index_1,'gene_id'].item()],
                            'transcript_id': [transc.loc[index_1,'transcript_id'].item()], 
                            'transcript_type': [transc.loc[index_1,'transcript_type'].item()], 
                            'intron_number': [int(i)]
                            })

                        introns_df = pd.concat([introns_df, intron_computed])
                        ##print(introns_df)

                        i += 1 
                        first_exon += 1


                    ##print(introns_df[['chr','transcript_id','start', 'end', 'intron_number']])
                    ##sys.exit()

                if each not in transcripts_with_introns:
                    transcripts_with_introns.append(each)

            progress += 1

            if progress % 10000 == 0: 
                print(progress, ' trasncripts')


            ## write the output -- per transcript
            introns_df.to_csv(output_introns, sep='\t', lineterminator='\n', index=False,
            header=header, mode='a')

            header = False

except ValueError:  ## made to catch the error -- see comments bellow
    print(each)
    print(num_exons)
    print(transc)

    ## the problem was that there was 2 entries in different chromosomes with the same transcript ID
    ##see: compute_GENCODE_introns_2022-09-21_errorCatch.log

    ## SOLUTION: run per chromosome, so this 


    ##     Traceback (most recent call last):
    ##  File "/Users/mariaf/Desktop/smorfs/5.check_var_effect/python/1-compute_introns_GENCODE_perTransc.py", line 99, in <module>
    ##    intron_start = transc.loc[index_1,'end'].item() + 1 # position after exon ends
    ##  File "/opt/homebrew/lib/python3.10/site-packages/pandas/core/base.py", line 349, in item  
    ##    raise ValueError("can only convert an array of size 1 to a Python scalar")
    ## ValueError: can only convert an array of size 1 to a Python scalar """


print(introns_df.shape)
#all_t_introns = introns_df['transcript_id'].unique()
#print(len(all_t_introns))

print('')
print('transcripts with introns: ', len(transcripts_with_introns))


end_time = (time.time() - start_time)/3600.0

print(end_time, ' hours')

print('DONE')