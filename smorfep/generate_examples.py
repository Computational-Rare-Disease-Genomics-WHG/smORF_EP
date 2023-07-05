#!/usr/bin/python3

## Contact: M Fernandes

## Script to generate examples of indels based on an interval and considering a sliding window

import sys

import argparse 
import random

from smorfep.utils.functions import *


### functions 

def choice_excluding(lst, exception):
    possible_choices = [item for item in lst if item != exception]
    return random.choice(possible_choices)



def main():
    """
    Main entry point
    """

    ## by default assumes BED file with no header

    parser = argparse.ArgumentParser(description='Script to generate all the variants within a given region (SNV and indels)')

    parser.add_argument('-c','--chrom', metavar='\b', required=True, type=str, help='chromosome')
    parser.add_argument('-os','--orfStart', metavar='\b', required=True, type=str, help='ORF start coordinate')
    parser.add_argument('-oe','--orfEnd', metavar='\b', required=True, type=str, help='ORF end coordinate')
    parser.add_argument('-s','--strand', metavar='\b', required=True, type=str, help='strand')
    parser.add_argument('-i','--orfID', metavar='\b', required=True, type=str, help='ORF ID')
    parser.add_argument('-ws','--windowStart', metavar='\b', required=True, type=int, help='variants region start coordinate')
    parser.add_argument('-we','--windowEnd', metavar='\b', required=True, type=int, help='variants region end coordinate')
    parser.add_argument('-m','--indelMaxSize', metavar='\b', required=True, type=str, help='indel maximum size')
    parser.add_argument('-r','--refpath', metavar='\b', required=True, type=str, help='Path to the reference genome')
    parser.add_argument('-o', '--outputfile', metavar='\b', required=True, type=str, help='output file name')
    parser.add_argument('-e', '--exonNts', metavar='\b', type=int, default=3, help='output file name')
    parser.add_argument('-j', '--intronNts', metavar='\b', type=int, default=8, help='output file name')

    args = parser.parse_args()

    ##print('python3 generate_examples.py 18 47808957 47930659 - SMAD2 47865134 47868322 2 test_exon-intron_auto_gen.vcf')


    all_nts = ['A', 'T', 'G', 'C']

    ## load the reeference genome chr sequence
    files_prefix, files_suffix = check_prefix_sufix_ref_files(args.ref_path)
    ref = read_single_fasta(args.chrom, args.ref_path, files_prefix, files_suffix)

    ## get sequence
    seq_left = get_sequence(args.windowStart-args.exonNts-1, args.windowStart+args.intronNts-1+args.indelMaxSize, '+', ref)
    print(args.windowStart, args.exonNts, args.intronNts)
    print(args.windowStart-args.exonNts-1, args.windowStart+args.intronNts-1)
    print(seq_left)

    ## output open and header
    out = open(args.outputfile, 'w')
    out.write('chrm\tvar_pos\tref\talt\tstart\tend\tstrand\tvar_id\tsmorf_id')

    ## NOTE: starts 1 postion before the exons nts -- as we use anchor notation
    for var_size in range(args.indelMaxSize+1): ## +1 to include the args.indelMaxSize
        print(var_size)
        diff = args.indelMaxSize-var_size
        iter_seq = seq_left[:len(seq_left)-diff]
        ##print(iter_seq)

        if var_size == 0: ## SNV -- runs only once
            for each_nt in range(len(iter_seq)-var_size):
                var_pos = each_nt + args.windowStart-args.exonNts-1
                ref = iter_seq[each_nt]
                alt = choice_excluding(all_nts, ref)
                var_type = 'SNV'

                ##print(var_pos, ref, alt, var_type)

                ## uncoment next line for SNV
                ##out.write('\n'+chrom+'\t'+str(var_pos)+'\t'+ref+'\t'+alt+'\t'+orf_start+'\t'+orf_end+'\t'+strand+'\t'+var_type+'\t'+orf_id)



        else: ## insertion and deletion  
            for each_nt in range(len(iter_seq)-var_size):
                var_pos = var_pos = each_nt + args.windowStart-args.exonNts-1

                ref_ins = iter_seq[each_nt]
                alt_ins = iter_seq[each_nt] + ''.join(random.choices(all_nts, k=var_size))
                var_type_ins = str(var_size) + 'nt_ins'

                ref_del = iter_seq[each_nt:each_nt+var_size+1]
                alt_del = iter_seq[each_nt]
                var_type_del = str(var_size) + 'nt_del'

                # print(var_pos, ref_ins, alt_ins, var_type_ins)
                # print(var_pos, ref_del, alt_del, var_type_del)
                out.write('\n'+args.chrom+'\t'+str(var_pos)+'\t'+ref_ins+'\t'+alt_ins+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_type_ins+'\t'+args.orfID)
                out.write('\n'+args.chrom+'\t'+str(var_pos)+'\t'+ref_del+'\t'+alt_del+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_type_del+'\t'+args.orfID)


    out.close()





if __name__ == '__main__':
    main()