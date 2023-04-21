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

        list_smorfs = small_df['smorf_id'].unique() ## TODO: check if this is a list or needs to be converted

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
