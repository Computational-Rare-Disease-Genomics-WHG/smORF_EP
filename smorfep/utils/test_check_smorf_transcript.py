#!/usr/bin/python3

## Author: M Fernandes

## Script to test check_smorf_transcript function 
## Apr 2023

## transcripts file is assumed to do not have duplicates.

import sys
import time
import pandas as pd

from smorfep.utils.functions import *



def test_check_smorf_transcript(ref_path, transcripts_filename, introns_filename, filename):

    ## 1- reads the input file
    variants_df = read_variants_file(filename, '\t', 0)

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
        ## variants/smorfs per chromosome
        small_df = variants_df.loc[variants_df['chrm'] == each_chrom]

        ## transcripts and introns in the chromosome
        transcripts_chr = transcripts_df.loc[transcripts_df['chr'] == 'chr'+str(each_chrom)]
        introns_small = introns_df.loc[introns_df['chr'] == 'chr'+str(each_chrom)]

        ## per smORF
        list_smorfs = small_df['smorf_id'].unique() ## list of unique smorf ID in the chromosome

        for each_smorf in list_smorfs: 
            smorf_df = small_df.loc[small_df['smorf_id'] == each_smorf]
            print(smorf_df)

            for index, row in smorf_df.iterrows():

                smorf_id = row.smorf_id
                smorf_start = row.start
                smorf_end = row.end
                smorf_strand = row.strand

                ## check 1 -- find the transcriptds that cover the full smorf region
                transcripts_small = transcripts_chr.loc[(transcripts_chr.start <= smorf_start) & (transcripts_chr.end >= smorf_end) & (transcripts_chr.strand == smorf_strand)]
                ## transcript needs to cover the full sequence region
                ## transcript in the same strand

                if not transcripts_small.empty:
                    
                    matching_t, unmatching_t, map_gen2transc, map_transc2gen = compatibility_smorf_transcript(reference_genome[each_chrom], transcripts_small, introns_small, smorf_start, smorf_end, smorf_strand)
                    ## this function runs for all transcripts within which the smorf falls within
                
                
            print(smorf_id)
            print(matching_t)
            print(unmatching_t.head)
            # print(map_gen2transc)
            # print(map_transc2gen)

            ##sys.exit(1)
                

def main(): 
    """
    Testing arguments

    """ 

    ref_path = "/Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/"
    transcripts_path = "/Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation.gff3_transcriptCoord_2023-04-14.tsv"
    introns_path = "/Users/mariaf/Desktop/GitHub/smORF_EP/transcripts/gencode.v41.annotation.gff3_introns_2023-04-14.tsv" 
    inputname = "/Users/mariaf/Desktop/GitHub/smORF_EP/smorfep/test/test14_Final_test.tsv"
    ## Ruby's file
    ##inputname = "/Users/mariaf/Desktop/GitHub/smORF_EP/smorfep/test/test_long_smorfinput.tsv"

    test_check_smorf_transcript(ref_path, transcripts_path, introns_path, inputname)


if __name__ == '__main__':

    start_time = time.time()

    main()

    print('DONE!')
    end_time = (time.time() - start_time)#/ 60.0
    print(end_time, 'seconds.')
