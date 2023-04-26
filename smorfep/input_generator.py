#!/usr/bin/python3

"""
Input generator script for smORF_EP
This script generates the formated file for smORF_EP tool. 
Merges the information from the variants and regions of interest (smORFs) into a dataframe and writes the output.

Usage: smorfinput [OPTIONS]

"""

import time

import pandas as pd
import argparse 

from smorfep.utils.functions import *



def bedvcf2intput(ref_path, bedfilename, vcffilename, outputname, bheader, vheader):

    """
    Compiles the input file read by smorfep.

    Input: 
    - ref_path: dir where the reference is stored
    - bedfile: BED file with the smORFs start and end coordinates
    - vcffile: VCF file with the variants 
    - output: output filename
    - bheader: 0 if bed file has header, None otherwise
    - vheader: 0 if vcf file has header, None otherwise
    """
    print("Generating smorfep input file...")

    start_time = time.time()

    ## 1- read bedfile - smorf regions
    smorfs_df = read_file(bedfilename, '\t', bheader)
    ## Keep onlyt the first 6 columns in the BED file -- info we need
    smorfs_df = smorfs_df.iloc[:, :6]
    ## rename the columns
    smorfs_df.columns = ['chrm', 'start', 'end', 'smorf_id', 'score', 'strand']
    ## rm chr prefix from chrm column
    smorfs_df['chrm'] = smorfs_df['chrm'].str.strip('chr')

    ## stats on the smorfs
    total_smorfs = smorfs_df.shape[0] ## number of smORFs in the chromosome
    print('#smORFs: ', total_smorfs) 

    ## 2- Loads reference genome
    ## the chromosomes the smORFs are in
    all_chromosomes = smorfs_df.chrm.unique()

    ## read reference
    ## stores reference per chromosome into a dictionary
    reference_genome = {} 


    ## check files prefix and suffix 
    files_prefix, files_suffix = check_prefix_sufix_ref_files(ref_path)
    
    for chrom_ref in all_chromosomes: 
        r = read_single_fasta(str(chrom_ref), ref_path, files_prefix, files_suffix)

        reference_genome[chrom_ref] = r

    print('reference ready')
    


    ## 3- Process variants
    ## generate a new df to store the new formated vars

    ## read variants file
    variants_df = read_file(vcffilename, '\t', vheader)
    variants_df = variants_df.iloc[:, :7]
    ## name columns
    variants_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']

    variants_df['CHROM'] = variants_df['CHROM'].str.strip('chr')

    ## unique index for variants
    vars_id_index = 1

    ## smorf index 
    smorf_index = 1

    first = True ## to print the header in the file only once
    smORFs_no_vars = 0 ## count the number of smORFs with no variants

    ## iterate per smORF
    for index, row in smorfs_df.iterrows():

        chrm_smorf = row.chrm
        start_smorf = row.start+1 ## +1 for the exact start coordinate as bed has start-1 format
        end_smorf = row.end
        smorfid = row.smorf_id
        strand_smorf = row.strand
        ##print(chrm_smorf, start_smorf, end_smorf, smorfid, strand_smorf)

        ## variants in the smORF
        smorf_variants_df = variants_df[(variants_df['CHROM'] == chrm_smorf) & (variants_df['POS']>= start_smorf) & (variants_df['POS'] <= end_smorf)]

        if not smorf_variants_df.empty:
            for index_var, row_var in smorf_variants_df.iterrows():
                ##print(row_var.POS, row_var.REF, row_var.ALT, row.chrm, row.start+1, row.end, row.strand)

                ## take variant ID from 3rd column in the VCF, if not empty
                if row_var.ID != '':
                    var_id = row_var.ID
                    ##print(var_id)

                else: ## if no ID provided generate one
                    var_id = 'VAR-' + str(vars_id_index)


                if strand_smorf == '+':
                    ## gnomad variants are on the forward strand, no special edits 
                    r = row_var.REF
                    a = row_var.ALT
                    new_var_pos = row_var.POS ## same position as reported

                elif strand_smorf == '-':
                    if len(row_var.REF) == len(row_var.ALT): ##SNV
                        r = complement_seq(row_var.REF)
                        a = complement_seq(row_var.ALT)
                        new_var_pos = row_var.POS ## same position as reported

                    elif len(row_var.REF) > len(row_var.ALT): ## deletion - OK
                        pos_diff = len(row_var.REF) - len(row_var.ALT)

                        a = get_sequence(int(row_var.POS)+pos_diff+1, int(row_var.POS)+pos_diff+1, strand_smorf, reference_genome[chrm_smorf])
                        ref_allele_sufix = reverse_complement_seq(row_var.REF)
                        r = a + ref_allele_sufix[:-1] ## removes the last nt
                        new_var_pos = int(row_var.POS)+pos_diff+1 ## var pos next position after the deletion section
                        

                    elif len(row_var.REF) < len(row_var.ALT): ## insertion 
                        r = get_sequence(int(row_var.POS)+1, int(row_var.POS)+1, strand_smorf, reference_genome[chrm_smorf])
                        alt_allele_sufix = reverse_complement_seq(row_var.ALT)
                        a = r + alt_allele_sufix[:-1] ## removes the last nt
                        new_var_pos = row_var.POS+1 ## same position as reported


                new_var = pd.DataFrame(
                    {'chrm': [chrm_smorf],
                    'var_pos': [new_var_pos],
                    'ref': [r],
                    'alt': [a],
                    'start': [start_smorf],
                    'end': [end_smorf],
                    'strand': [strand_smorf], 
                    'var_id': [var_id], 
                    'smorf_id' : [smorfid]           
                    })

                df = pd.DataFrame(new_var) ## creates a temporary dataframe

                vars_id_index += 1

                if first == True: ## in the first write, write the header line
                    df.to_csv(outputname, sep='\t', lineterminator='\n', index=False, header=True)
                    first = False
                else:
                    # append data frame to CSV file
                    df.to_csv(outputname, mode='a', sep='\t', lineterminator='\n', index=False, header=False)

        else: ## print smORF ID witout variants in it
            #print(smorfid)
            smORFs_no_vars += 1
        
        smorf_index += 1

        if smorf_index % 10000 == 0: 
            print(smorf_index, 'smorfs processed')

    print('#smORFs without vars: ', smORFs_no_vars) 
    print('DONE!')

    end_time = (time.time() - start_time)/ 60.0
    print(end_time, ' minutes.')

    return None



def main():
    """
    Main entry point
    """

    ## by default assumes BED file with no header

    parser = argparse.ArgumentParser(description='Script to convert BED and VCF into smorfep input file')

    parser.add_argument('-r','--refpath', required=True, type=str, help='Path to the reference genome')
    parser.add_argument('-b','--bedfile', required=True, type=str, help='BED file with the smORFs regions')
    parser.add_argument('-v', '--vcffile', required=True, type=str, help='VCF file with the variants')
    parser.add_argument('-o', '--outputfile', required=True, type=str, help='output file name')
    parser.add_argument('--bedheader', help='BED file: first line is the header', action='store_true')
    parser.add_argument('--vcfheader', help='VCF file: first line is the header', action='store_true')


    args = parser.parse_args()

    if args.bedheader:
        bedheader = 0
    else: 
        bedheader = None
    if args.vcfheader: 
        vcfheader = 0
    else: 
        vcfheader = None


    ## generate the input file: var-smorf pairs
    bedvcf2intput(args.refpath, args.bedfile, args.vcffile, args.outputfile, bedheader, vcfheader)


if __name__ == '__main__':
    main()