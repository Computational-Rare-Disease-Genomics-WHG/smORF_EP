#!/usr/bin/python3

## Contact: M Fernandes



import sys
import os
import time

import pandas as pd
import argparse 

from os.path import exists
from datetime import date
from numpy import require

from smorfep.utils.functions import read_variants_file, read_file, genome2smorf_coords
from smorfep.utils.tool_script_only import *



def run_smorfep_smorfs(ref_path, introns_filename, splice_site, filename, outputname):

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

    ## 1- Import introns 
    introns_df = read_file(introns_filename, '\t', 0) ## computed at the begining of the script
    ##print(introns_df.shape)
    ## controls for duplicates if any - Added on 2024-11-05 - OK
    introns_df = introns_df.drop_duplicates()
    ##print('shape after removing duplicates')
    ##print(introns_df.shape)
    print('introns ready')
    print('')

    ## runs per chromosome
    for each_chrom in all_chromosomes: ## runs per chromosome
        ## TODO: OPTIMIZE --> allow cache freeing after each chromosome -- remove chromosome from the ref_genome dictionary

        ## variants/smorfs per chromosome
        small_df = variants_df.loc[variants_df['chrm'] == each_chrom]
        ## TODO: OPTIMIZE --> After each smORF, we can also remove the smORF from the analysis table -- Write outputfile before

        ## TODO: OPTIMIZE --> allow cache freeing after each chromosome -- remove chromosome from the ref_genome dictionary
        
        ## per smORF
        ## 1 - Collect the list of smORFs
        ## 2 - iterate per smORFs
        ## 3 - collect the transcript info per smORF 
        ## 4 - check_smorf_transcript
        ##      4.1 check trancript info first and exclude the unmatching transcripts first; 
        ##      4.2 run the tool for the matching transcripts only
        ## 5 - Check each variant in the smORF -- effect per transcript

        list_smorfs = small_df['smorf_id'].unique() ## list

        ## itarete per smorf ID
        for smorf_id in list_smorfs:
            print(smorf_id)
            smorf_vars_df = small_df[small_df['smorf_id'] == smorf_id]
            ##print(smorf_vars_df)

            ## smorf info
            smorf_start = smorf_vars_df.iloc[0]['start']
            smorf_end = smorf_vars_df.iloc[0]['end']
            smorf_strand = smorf_vars_df.iloc[0]['strand']
            ##print(smorf_id, smorf_start, smorf_end, smorf_strand)

            ## introns for the smORF
            introns_smORF = introns_df.loc[introns_df['ID'] == smorf_id]

            ## mapping for this smORF
            map_gen2smorf, map_smorf2gen = genome2smorf_coords(smorf_start, smorf_end, smorf_strand, introns_df)

            ## per variant - line
            for index, row in smorf_vars_df.iterrows():
                ##print(index)

                ## row_t is the info about the transctipt
                consequence, change, prot_cons, prot_change = tool(
                    reference_genome[each_chrom], 
                    introns_smORF, 
                    smorf_start,
                    smorf_end,
                    smorf_strand,
                    row.ref, 
                    row.alt, 
                    row.var_pos, 
                    map_gen2smorf, 
                    map_smorf2gen,
                    splice_site)

                ##r_index = smorf_vars_df.index[smorf_vars_df['var_id'] == row.var_id].item()
                ## variant IDs are unique - so we only get one index out of this
                ##print(r_index)
                ##print( smorf_vars_df.index)

                ## adds to the dataframe the output 
                consequence_computed = pd.DataFrame(
                    {
                    'chrm': variants_df.iloc[index]['chrm'],
                    'var_pos' : variants_df.iloc[index]['var_pos'],
                    'ref' : variants_df.iloc[index]['ref'],
                    'alt' : variants_df.iloc[index]['alt'],
                    'start' : variants_df.iloc[index]['start'],
                    'end' : variants_df.iloc[index]['end'],
                    'strand' : variants_df.iloc[index]['strand'],
                    'var_id' : variants_df.iloc[index]['var_id'],
                    'smorf_id': variants_df.iloc[index]['smorf_id'], 
                    'DNA_consequence' : consequence,
                    'DNA_seq' : change,
                    'prot_consequence' : prot_cons,
                    'prot_seq' : prot_change
                    }, index=[index]
                )

                vars_cons_df = pd.concat([vars_cons_df, consequence_computed])

                ## NOTE 1: Removed no transcript -- we report the smORF IDS for smORFs without transcript but don't run the analysis for those

                ## NOTE 2: Removed var outside smorf region -- as we remove the non-matching transcripts from analysis
                ## they are in a separate file


    ##print(vars_cons_df)
            
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
    parser.add_argument('-i', '--introns_filename', metavar='\b', required=True, type=str, help='smORFs introns file. IDs must correspond')
    parser.add_argument('-s', '--splice_site', metavar='\b', type=int, default=8, help='splice-site size, default = 8 (as VEP)')
    parser.add_argument('-f', '--variants_filename',metavar='\b', required=True, type=str, help='file with the variants and the regions of interest info')
    parser.add_argument('-o', '--output', metavar='\b', type=str, default=default_outputname, help='outputname')


    args = parser.parse_args()


    start_time = time.time()


    ## if the introns file do now exist in the repository, it is computed
    ## takes less than 3 min
    if os.path.exists(args.introns_filename):  
        print('smORF introns file available')

    else: 
        print('smORF introns file NOT available\n')
        print('Make use the file is available at the given location.')
        sys.exit(1)
 


    ## list of files in the reference dir
    ref_dir_files = os.listdir(args.reference_path)

    ## NOTE: hard coded -- suffix used by smorfinit step
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
    run_smorfep_smorfs(args.reference_path, 
                args.introns_filename, 
                args.splice_site, 
                args.variants_filename, 
                args.output)

    ## TODO: Add stats on how many smorfs skipped and how many variants successfully annotated (out of XXX)

    print('DONE!')

    end_time = (time.time() - start_time)#/ 60.0

    print(end_time, 'seconds.')



if __name__ == '__main__':
    main()
