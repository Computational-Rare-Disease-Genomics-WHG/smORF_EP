#!/usr/bin/python3

## Author: M Fernandes

## Script to test check_smorf_transcript function 
## Apr 2023

import sys
import os
import time

import pandas as pd
import argparse 

from os.path import exists
from datetime import date
from numpy import require

from smorfep.utils.functions import *



def test_check_smorf_transcript(ref_path, transcripts_filename, introns_filename, filename):

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
        list_smorfs = small_df['smorf_id'].unique() ## TODO: check if this is a list or needs to be converted
        print(list_smorfs)

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

            ## Unique IDs from transcripts_small
            t_unique_ids = transcripts_small['transcript_id'].unique()

            if not transcripts_small.empty:
                
                ##4.2 - iterare per transcript
                for index_t, row_t in transcripts_small.iterrows():

                    ## introns per transcript
                    introns_transcript = introns_small.loc[introns_small['transcript_id'] == row_t.transcript_id]

                    check_smorf_transcript(reference_genome[each_chrom], transcripts_small, introns_transcript)

                



def main(): 
    """
    Testing arguments

    """ 

    ref_path = "/Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/"
    transcripts_path = "/Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation.gff3_transcriptCoord_2023-04-14.tsv"
    introns_path = "/Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation.gff3_introns_2023-04-14.tsv" 
    inptuname = "/Users/mariaf/Desktop/GitHub/smORF_EP/smorfep/test/test15_Final_test_introns.tsv"


    test_check_smorf_transcript(ref_path, transcripts_path, introns_path, inptuname)



if __name__ == '__main__':
    main()