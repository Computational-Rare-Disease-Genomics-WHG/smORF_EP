"""
Input generator script for smORF_EP
This script generates the formated file for smORF_EP tool. 
Merges the information from the variants and regions of interest (smORFs) into a dataframe and writes the output.

Usage: smorfinput [OPTIONS]

"""

import time
import datetime

import pandas as pd
import argparse 

from smorfep.utils.functions import *



def bedvcf2intput(bedfilename, vcffilename, outputname):
    """
    Compiles the input file read by smorfep.

    Input: 
    - bedfile: BED file with the smORFs start and end coordinates
    - vcffile: VCF file with the variants 
    """
    print("Generating smorfep input file...")

    start_time = time.time()

    ## read bedfile - smorf regions
    smorfs_df = read_file(bedfilename, '\t', None)
    smorfs_df.columns = ['chrm', 'start', 'end', 'smORF_id', 'score', 'strand']
    ## rm chr prefix from chrm column
    smorfs_df['chrm'] = smorfs_df['chrm'].str.strip('chr')


    ## generate a new df to store the new formated vars
    vars_df = pd.DataFrame(data=None, columns=['chrm','var_pos','ref','alt','start','end','strand','var_id', 'smORF_id'])

    chrom_smorf_df = smorfs_df[smorfs_df['chrm'] == c] ## consider only smORFs in the same chromosome ß
    total_smorfs = chrom_smorf_df.shape[0] ## number of smORFs in the chromosome
    print('#smORFs: ', total_smorfs) 
    ##print(chrom_smorf_df)


    ## read variants file
    variants_df = read_file(vcffilename, '\t', 0)
    variants_df = variants_df.drop(['QUAL', 'FILTER', 'AC', 'AN', 'AF'], axis=1)
    variants_df['CHROM'] = variants_df['CHROM'].str.strip('chr')
    ##print(variants_df)

    ## unique index for variants
    vars_id_index = 1

    ## smorf index 
    smorf_index = 1

    ## output name
    now = datetime.datetime.now().strftime('%Y-%m-%d')
    outputname = outputname_prefix + c + '_input_'+ today.strftime("%Y-%m-%d") +'.tsv'

    first = True

    smORFs_no_vars = 0 ## count the number of smORFs with no variants


    ## iterate per smORF
    for index, row in chrom_smorf_df.iterrows():
        chrm_smorf = row.chrm
        start_smorf = row.start+1 ## +1 for the exact start coordinate as bed has start-1 format
        end_smorf = row.end
        smorfid = row.smORF_id
        strand_smorf = row.strand
        ##print(chrm_smorf, start_smorf, end_smorf, smorfid, strand_smorf)

        ## variants in the smORF
        smorf_variants_df = variants_df[(variants_df['CHROM'] == chrm_smorf) & (variants_df['POS']>= start_smorf) & (variants_df['POS'] <= end_smorf)]
        ##print(smorf_variants_df)

        if not smorf_variants_df.empty:
            for index_var, row_var in smorf_variants_df.iterrows():
                var_id = 'XX' + str(vars_id_index)

                ## testing block
                # if row_var.POS == 6228696:
                #     print(6228696)
                #     print(row_var.REF, row_var.ALT)
                #     print(smorfid)
                #     print(strand_smorf)
                #     print(start_smorf, end_smorf)


                if strand_smorf == '+':
                    ## gnomad variants are on the forward strand, no special edits 
                    r = row_var.REF
                    a = row_var.ALT
                    new_var_pos = row_var.POS ## same position as reported

                elif strand_smorf == '-':
                    ##print('here')
                    if len(row_var.REF) == len(row_var.ALT): ##SNV
                        ##print('SNV')
                        r = complement_seq(row_var.REF)
                        a = complement_seq(row_var.ALT)
                        new_var_pos = row_var.POS ## same position as reported

                    elif len(row_var.REF) > len(row_var.ALT): ## deletion - OK
                        ##print('deletion')
                        pos_diff = len(row_var.REF) - len(row_var.ALT)
                        a = get_sequence(str(chrm_smorf), int(row_var.POS)+pos_diff+1, int(row_var.POS)+pos_diff+1, strand_smorf)
                        ref_allele_sufix = reverse_complement_seq(row_var.REF)
                        r = a + ref_allele_sufix[:-1] ## removes the last nt

                        new_var_pos = int(row_var.POS)+pos_diff+1 ## var pos next position after the deletion section
                        

                    elif len(row_var.REF) < len(row_var.ALT): ## insertion 
                        ##print('insertion')
                        r = get_sequence(str(chrm_smorf), int(row_var.POS)+1, int(row_var.POS)+1, strand_smorf)
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
                    'smORF_id' : [smorfid]           
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
            print(smorfid)
            smORFs_no_vars += 1
        
        smorf_index += 1

        if smorf_index % 20 == 0: 
            print(smorf_index)


    print('#smORFs NO vars: ', smORFs_no_vars) 


    print('DONE!')



    end_time = (time.time() - start_time)/ 60.0

    print(end_time, ' minutes.')






    

def main():
    """
    Main entry point
    """

    parser = argparse.ArgumentParser(description='Script to convert BED and VCF into smorfep input file')


    args = parser.parse_args()


    if args.reference: 
        download_ref_genome(args.ref_link)

    elif args.transcripts: 
        download_gencode(args.transc_link)

    elif args.all: 
        download_ref_genome(args.ref_link)
        download_gencode(args.transc_link)



if __name__ == '__main__':
    main()