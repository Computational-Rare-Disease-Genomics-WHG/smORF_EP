#!/usr/bin/python3

"""
Input generator script for smORF_EP
This script generates the formatted file for smORF_EP tool. 
Merges the information from the variants and regions of interest (smORFs) into a dataframe and writes the output.

Usage: smorfinput [OPTIONS]

"""
import time, argparse
import pandas as pd

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
    ## Keep only the first 6 columns in the BED file -- info we need
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

    print('reference ready (new)')

    ## 3- Process variants
    ## generate a new df to store the new formatted vars
    ## read variants file
    variants_df = read_file(vcffilename, '\t', vheader)
    variants_df = variants_df.iloc[:, :7]
    ## name columns
    variants_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    variants_df['CHROM'] = variants_df['CHROM'].str.strip('chr')

    ## unique indices for variants, smorfs w/o variants
    vars_id_index, smorf_index, smORFs_no_vars = [1], [1], []
    df_chrm, df_var_pos, df_ref, df_alt, df_start, df_end, df_strand, df_varid, df_smorfid = [], [], [], [], [], [], [], [], []
    ## iterate per smORF
    
    def parse_smorf(chrm, start, end, smorf_id, strand):
        chrm_smorf, start_smorf, end_smorf, smorfid, strand_smorf = chrm, start+1, end, smorf_id, strand    ## +1 for the exact start coordinate as bed has start-1 format
        smorf_variants_df = variants_df[(variants_df['CHROM'] == chrm_smorf) & (variants_df['POS']>= start_smorf) & (variants_df['POS'] <= end_smorf)]

        if not smorf_variants_df.empty:
            def parse_variants(POS, REF, ALT, ID, chrm_smorf, start_smorf, end_smorf, smorfid, strand_smorf):
                ## take variant ID from 3rd column in the VCF, if not empty
                var_id = ID if ID != '' else f'VAR-{len(vars_id_index)}'

                if strand_smorf == '+':
                    ## gnomad variants are on the forward strand, no special edits 
                    r, a, new_var_pos = REF, ALT, POS   #same position as reported
                elif strand_smorf == '-':
                    if len(REF) == len(ALT): ##SNV
                        r, a, new_var_pos = complement_seq(REF), complement_seq(ALT),  POS ## same position as reported
                    elif len(REF) > len(ALT): ## deletion - OK
                        pos_diff = len(REF) - len(ALT)
                        a = get_sequence(int(POS) + pos_diff + 1, int(POS) + pos_diff + 1, strand_smorf, reference_genome[chrm_smorf])
                        ref_allele_sufix = reverse_complement_seq(REF)
                        r = a + ref_allele_sufix[:-1] ## removes the last nt
                        new_var_pos = int(POS) + pos_diff + 1 ## var pos next position after the deletion section
                    elif len(REF) < len(ALT): ## insertion 
                        r = get_sequence(int(POS) + 1, int(POS) + 1, strand_smorf, reference_genome[chrm_smorf])
                        alt_allele_sufix = reverse_complement_seq(ALT)
                        a = r + alt_allele_sufix[:-1] ## removes the last nt
                        new_var_pos = POS + 1 ## same position as reported

                df_chrm.append(chrm_smorf)
                df_var_pos.append(new_var_pos)
                df_ref.append(r)
                df_alt.append(a)
                df_start.append(start_smorf)
                df_end.append(end_smorf)
                df_strand.append(strand_smorf)
                df_varid.append(var_id)
                df_smorfid.append(smorfid)
                vars_id_index.append(1)

            pd.options.mode.chained_assignment = None  # default='warn'
            smorf_variants_df['chrm_smorf'], smorf_variants_df['start_smorf'], smorf_variants_df['end_smorf'], smorf_variants_df['smorfid'], smorf_variants_df['strand_smorf'] = chrm_smorf, start_smorf, end_smorf, smorfid, strand_smorf
            [parse_variants(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8]) for a in zip(smorf_variants_df['POS'], smorf_variants_df['REF'], smorf_variants_df['ALT'], smorf_variants_df['ID'], \
                                                            smorf_variants_df['chrm_smorf'], smorf_variants_df['start_smorf'], smorf_variants_df['end_smorf'], smorf_variants_df['smorfid'], smorf_variants_df['strand_smorf'])]
        else:
            smORFs_no_vars.append(1)

        smorf_index.append(1)
        if len(smorf_index) % 5000 == 0: 
            print(len(smorf_index), 'smorfs processed')

    [parse_smorf(a[0], a[1], a[2], a[3], a[4]) for a in zip(smorfs_df['chrm'], smorfs_df['start'], smorfs_df['end'], smorfs_df['smorf_id'], smorfs_df['strand'])]
    parse_variants_df = pd.DataFrame(columns=['chrm', 'var_pos', 'ref', 'alt', 'start', 'end', 'strand', 'var_id', 'smorf_id'])
    parse_variants_df['chrm'], parse_variants_df['var_pos'], parse_variants_df['ref'], parse_variants_df['alt'] = df_chrm, df_var_pos, df_ref, df_alt
    parse_variants_df['start'], parse_variants_df['end'], parse_variants_df['strand'] = df_start, df_end, df_strand
    parse_variants_df['var_id'], parse_variants_df['smorf_id'] = df_varid, df_smorfid
    parse_variants_df.to_csv(outputname, sep='\t', lineterminator='\n', index=False, header=True)
    print(f'#smORFs without vars: {len(smORFs_no_vars)}\nDONE!\n{time.time() - start_time} seconds.')
    return None

def main():
    """
    Main entry point
    """
    ## by default assumes BED file with no header
    parser = argparse.ArgumentParser(description='Script to convert BED and VCF into smorfep input file')
    parser.add_argument('-r', '--refpath', required=True, type=str, help='Path to the reference genome')
    parser.add_argument('-b', '--bedfile', required=True, type=str, help='BED file with the smORFs regions')
    parser.add_argument('-v', '--vcffile', required=True, type=str, help='VCF file with the variants')
    parser.add_argument('-o', '--outputfile', required=True, type=str, help='output file name')
    parser.add_argument('--bedheader', help='BED file: first line is the header', action='store_true')
    parser.add_argument('--vcfheader', help='VCF file: first line is the header', action='store_true')

    args = parser.parse_args()

    bedheader = 0 if args.bedheader else None
    vcfheader = 0 if args.vcfheader else None
    ## generate the input file: var-smorf pairs
    bedvcf2intput(args.refpath, args.bedfile, args.vcffile, args.outputfile, bedheader, vcfheader)

if __name__ == '__main__':
    main()