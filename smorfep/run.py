#!/usr/bin/python3

## Author: M Fernandes

## Script to check variant consequence with relation to smORF regions.

## Possible effects:
## start lost (genomic)
## synonymous variant (protein)
## missense variant (protein)
## inframe insertion (genomic)
## inframe deletion (genomic)
## stop gain (genomic)
## frameshift variant (genomic)
## splice donor variant
## splice acceptor variant
## splice region variant
## intron variant
## stop lost (genomic)


import sys
import os
import time

import pandas as pd
import argparse 

from os.path import exists
from datetime import date
from numpy import require

from smorfep.utils.functions import *
from smorfep.utils.tool_script import *




def run_smorfep(ref_path, transcripts_filename, introns_filename, splice_site, filename, outputname):

    ## 1- reads the input file
    variants_df = read_variants_file(filename, '\t', 0)

    ## creates new empty dataframe for the consequences
    vars_cons_df = pd.DataFrame(data=None, columns=variants_df.columns)

    ## checks which cromosomes are in it
    all_chromosomes = variants_df.chrm.unique()

    ## read reference
    ## stores reference per chromosome into a dictionary
    reference_genome = {} ## stores the sequence per chromosome

    ## check files prefix and suffix 
    files_prefix, files_suffix = check_prefix_sufix_ref_files(ref_path)
    
    for chrom_ref in all_chromosomes: 
        r = read_single_fasta(str(chrom_ref), ref_path, files_prefix, files_suffix)

        reference_genome[chrom_ref] = r

    print('reference ready')

    ## 1- Import transcripts
    transcripts_df = read_file(transcripts_filename, '\t', 0)
    print('transcripts ready')

    ## 2- Import introns 
    introns_df = read_file(introns_filename, '\t', 0) ## computed at the begining of the script
    print('introns ready')
    print('')


    ## 4- Check variant effect per transcript
    for each_chrom in all_chromosomes: ## runs per chromosome
        ## TODO: OPTIMIZE --> allow cache freeing after each chromosome -- remove chromosome from the ref_genome dictionary

        ## variants/smorfs per chromosome
        small_df = variants_df.loc[variants_df['chrm'] == each_chrom]
        ## TODO: OPTIMIZe --> After each smORF, we can also remove the smORF from the analysis table -- Write outputfile before

        ## transcripts and introns in the chromosome
        transcripts_chr = transcripts_df.loc[transcripts_df['chr'] == 'chr'+str(each_chrom)]
        introns_small = introns_df.loc[introns_df['chr'] == 'chr'+str(each_chrom)]
        ## TODO: OPTIMIZE --> allow cache freeing after each chromosome -- remove chromosome from the ref_genome dictionary

        ## per smORF
        ## 1- Collect hte list of smORFs
        ## 2 - iterate per smORFs
        ## 3- collect the transcript info per smORF -- check_smorf_transcript
        ## 4- Check each variant in the smORF -- 4.1 check trancript info first and exclude the unmatching transcripts first; 4.2 run the tool for the matching transcripts

        list_smorfs = small_df[]

        ## per variant
        for index, row in small_df.iterrows(): ## iterates per line 
            
            ##4.1 - find transcripts the region falls in:
            variant_position = small_df.loc[index]['var_pos']
            seq_start = small_df.loc[index]['start']
            seq_end = small_df.loc[index]['end']
            seq_strand = small_df.loc[index]['strand']
            variant_id = small_df.loc[index]['var_id']


            transcripts_small = transcripts_chr.loc[(transcripts_chr.start <= seq_start) & (transcripts_chr.end >= seq_end) & (transcripts_chr.strand == seq_strand)]
            ## transcript needs to cover the full sequence region
            ## transcript in the same strand

            ## check if the variant is within the region of interest, only run the tool if it is within
            if variant_position >= seq_start and variant_position <= seq_end: 

                if not transcripts_small.empty:
                    
                    ##4.2 - check the consequence per transcript
                    for index_t, row_t in transcripts_small.iterrows():

                        ## introns per transcript
                        introns_transcript = introns_small.loc[introns_small['transcript_id'] == row_t.transcript_id]

                        ## row_t is the info about the transctipt
                        consequence, change, prot_cons, prot_change = tool(
                            reference_genome[each_chrom], 
                            row_t, 
                            introns_transcript, 
                            row.start,
                            row.end,
                            row.strand,
                            row.ref, 
                            row.alt, 
                            row.var_pos, 
                            splice_site)

                        r_index = variants_df.index[variants_df['var_id'] == row.var_id].tolist()
                        ## variant IDs are unique - so we only get one index out of this
                        
                        ## adds to the dataframe the protein consequences
                        consequence_computed = pd.DataFrame(
                            {
                            'chrm': variants_df.iloc[r_index]['chrm'],
                            'var_pos' : variants_df.iloc[r_index]['var_pos'],
                            'ref' : variants_df.iloc[r_index]['ref'],
                            'alt' : variants_df.iloc[r_index]['alt'],
                            'start' : variants_df.iloc[r_index]['start'],
                            'end' : variants_df.iloc[r_index]['end'],
                            'strand' : variants_df.iloc[r_index]['strand'],
                            'var_id' : variants_df.iloc[r_index]['var_id'],
                            'smorf_id': variants_df.iloc[r_index]['smorf_id'], ## TODO: generalise if used for other regions
                            'transcript_id' : row_t.transcript_id, 
                            'transcript_type' : row_t.transcript_type,
                            'DNA_consequence' : consequence,
                            'DNA_seq' : change,
                            'prot_consequence' : prot_cons,
                            'prot_seq' : prot_change
                            }

                        )

                        vars_cons_df = pd.concat([vars_cons_df, consequence_computed])

                else:  ## write in the output when NO transcript overlaps the region in study
                    r_index = variants_df.index[variants_df['var_id'] == row.var_id].tolist()
                        
                    ## adds to the dataframe the protein consequences
                    consequence_computed = pd.DataFrame(
                        {
                        'chrm': variants_df.iloc[r_index]['chrm'],
                        'var_pos' : variants_df.iloc[r_index]['var_pos'],
                        'ref' : variants_df.iloc[r_index]['ref'],
                        'alt' : variants_df.iloc[r_index]['alt'],
                        'start' : variants_df.iloc[r_index]['start'],
                        'end' : variants_df.iloc[r_index]['end'],
                        'strand' : variants_df.iloc[r_index]['strand'],
                        'var_id' : variants_df.iloc[r_index]['var_id'],
                        'smorf_id': variants_df.iloc[r_index]['smorf_id'],
                        'transcript_id' : 'no_transcript_full_region', 
                        'transcript_type' : '-',
                        'DNA_consequence' : '-',
                        'DNA_seq' :'-',
                        'prot_consequence' : '-',
                        'prot_seq' : '-'
                        }

                    )

                    vars_cons_df = pd.concat([vars_cons_df, consequence_computed])
                    
            

            else: ## var outside the region of interest
                consequence, change, prot_cons, prot_change = 'variant out of region', '-', '-', '-'

                r_index = variants_df.index[variants_df['var_id'] == row.var_id].tolist()

                consequence_computed = pd.DataFrame(
                    {
                    'chrm': variants_df.iloc[r_index]['chrm'],
                    'var_pos' : variants_df.iloc[r_index]['var_pos'],
                    'ref' : variants_df.iloc[r_index]['ref'],
                    'alt' : variants_df.iloc[r_index]['alt'],
                    'start' : variants_df.iloc[r_index]['start'],
                    'end' : variants_df.iloc[r_index]['end'],
                    'strand' : variants_df.iloc[r_index]['strand'],
                    'var_id' : variants_df.iloc[r_index]['var_id'],
                    'transcript_id' : '-', 
                    'transcript_type' : '-',
                    'DNA_consequence' : consequence,
                    'DNA_seq' : change,
                    'prot_consequence' : prot_cons,
                    'prot_seq' : prot_change
                    }

                )
                vars_cons_df = pd.concat([vars_cons_df, consequence_computed])

    ## write_the output
    vars_cons_df.to_csv(outputname, sep='\t', lineterminator='\n', index=False)


def main():
    ## day date
    today = date.today()
    default_outputname = 'output_'+today.strftime("%Y-%m-%d")+'.tsv'
    ##print(default_outputname)

    parser = argparse.ArgumentParser(description='Script to annotate variants within small open reading frames (smORFs)')

    ## define arguments
    parser.add_argument('-r', '--reference_path', metavar='\b', required=True, type=str, help='reference genome dir path')
    parser.add_argument('-t', '--transcripts_filename', metavar='\b', required=True, type=str, help='transcripts database file, .tsv, generated by compute transcripts step')
    parser.add_argument('-i', '--introns_filename', metavar='\b', required=True, type=str, help='introns file, .tsv file, generate by compute introns step')
    parser.add_argument('-s', '--splice_site', metavar='\b', type=int, default=8, help='splice-site size, default = 8 (as VEP)')
    parser.add_argument('-f', '--variants_filename',metavar='\b', required=True, type=str, help='file with the variants and the regions of interest info')
    parser.add_argument('-o', '--output', metavar='\b', type=str, default=default_outputname, help='outputname')

    args = parser.parse_args()


    start_time = time.time()


    ## if the introns file do now exist in the repository, it is computed
    ## takes less than 3 min
    if os.path.exists(args.introns_filename):  
        print('introns file available')

    else: 
        print('introns file NOT available\n')
        print('RUN: smorfinit --transcripts --transc_link <link>\n')
        print('This command generates both transcript and introns files')
        sys.exit(1)
 

    if os.path.exists(args.transcripts_filename):
        print('transcripts coordinates file available')
    else: 
        print('introns not transcript coordinates file NOT available\n')
        print('RUN: smorfinit --transcripts --transc_link <link>\n')
        print('This command generates both transcript and introns files')
        sys.exit(1)

    ## list of files in the reference dir
    ref_dir_files = os.listdir(args.reference_path)
    chrom_list = ['chr1','chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8',\
                'chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16',\
                'chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    matching = [each for each in ref_dir_files if any(ch in each for ch in chrom_list)]
    if len(matching) == 24: ## checks if there are 24 files matching the chromosomes in the reference genome directory
        print('reference genome available')
    else: 
        print('reference genome NOT available\n')
        print('RUN: smorfinit -reference --ref_link <link>')
        sys.exit(1)

    ## run code
    run_smorfep(args.reference_path, args.transcripts_filename,
                args.introns_filename, args.splice_site, 
                args.variants_filename, args.output)

    print('DONE!')

    end_time = (time.time() - start_time)#/ 60.0

    print(end_time, 'seconds.')



if __name__ == '__main__':
    main()
