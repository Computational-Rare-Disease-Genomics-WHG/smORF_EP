#!/usr/bin/python3

## Contact: M Fernandes

## Script to generate examples of indels of a given lenght and considering a sliding window

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
    parser.add_argument('-ws','--intronStart', metavar='\b', required=True, type=int, help='intron start coordinate')
    parser.add_argument('-we','--intronEnd', metavar='\b', required=True, type=int, help='intro end coordinate')
    parser.add_argument('-m','--indelSize', metavar='\b', required=True, type=int, help='indel size')
    parser.add_argument('-r','--refpath', metavar='\b', required=True, type=str, help='Path to the reference genome')
    parser.add_argument('-o', '--outputfile', metavar='\b', required=True, type=str, help='output filename prefix. Default format: .vcf')
    parser.add_argument('-e', '--exonNts', metavar='\b', type=int, default=3, help='output file name')
    parser.add_argument('-j', '--intronNts', metavar='\b', type=int, default=8, help='output file name')

    args = parser.parse_args()

    ##print('python3 generate_examples.py 18 47808957 47930659 - SMAD2 47865134 47868322 2 test_exon-intron_auto_gen.vcf')


    all_nts = ['A', 'T', 'G', 'C']

    ## load the reeference genome chr sequence
    files_prefix, files_suffix = check_prefix_sufix_ref_files(args.refpath)
    ref = read_single_fasta(args.chrom, args.refpath, files_prefix, files_suffix)

    ## get sequence
    seq_left = get_sequence(args.intronStart-args.exonNts-1, args.intronStart+args.intronNts-1+args.indelSize, '+', ref)
    # print(args.windowStart, args.exonNts, args.intronNts)
    # print(args.windowStart-args.exonNts-1, args.windowStart+args.intronNts-1)
    # print(seq_left)

    ## get sequence right end
    seq_right = get_sequence(args.intronEnd-args.intronNts-1, args.intronEnd+args.exonNts-1+args.indelSize, '+', ref)

    ## output open and header
    if args.strand == '+':
        out_left = open(args.outputfile+'_donor.vcf', 'w')
        out_left.write('chrm\tvar_pos\tref\talt\tstart\tend\tstrand\tvar_id\tsmorf_id')
        
        out_right = open(args.outputfile+'_acceptor.vcf', 'w')
        out_right.write('chrm\tvar_pos\tref\talt\tstart\tend\tstrand\tvar_id\tsmorf_id')
    elif args.strand == '-':
        out_right = open(args.outputfile+'_donor.vcf', 'w')
        out_right.write('chrm\tvar_pos\tref\talt\tstart\tend\tstrand\tvar_id\tsmorf_id')
        
        out_left = open(args.outputfile+'_acceptor.vcf', 'w')
        out_left.write('chrm\tvar_pos\tref\talt\tstart\tend\tstrand\tvar_id\tsmorf_id')

    var_index_left = 1
    var_id_left = 'XXX'+str(var_index_left)

    var_index_right = 1
    var_id_right = 'XXX'+str(var_index_right)

    iter_seq_left = seq_left[:len(seq_left)-args.indelSize] ## left side
    iter_seq_right = seq_right[:len(seq_right)-args.indelSize] ## right side


    for each_nt_left in range(len(iter_seq_left)-args.indelSize):
        var_pos_left = each_nt_left + args.intronStart-args.exonNts-1 ## left side


        ref_del = iter_seq_left[each_nt_left:each_nt_left+args.indelSize+1]
        alt_del = iter_seq_left[each_nt_left]
        var_type_del = str(args.indelSize) + 'nt_del'
        out_left.write('\n'+args.chrom+'\t'+str(var_pos_left)+'\t'+ref_del+'\t'+alt_del+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_id_left +'_'+var_type_del+'\t'+args.orfID)
        var_index_left +=1
        var_id_left = 'XXX'+str(var_index_left)


        ref_ins = iter_seq_left[each_nt_left]
        all_alt_ins_left = generate_permutations(all_nts, args.indelSize)
        for each_left_alt in all_alt_ins_left:
            alt_ins = iter_seq_left[each_nt_left] + each_left_alt
 
            var_type_ins = str(args.indelSize) + 'nt_ins'

            out_left.write('\n'+args.chrom+'\t'+str(var_pos_left)+'\t'+ref_ins+'\t'+alt_ins+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_id_left +'_'+var_type_ins+'\t'+args.orfID)
            var_index_left +=1
            var_id_left = 'XXX'+str(var_index_left)


        ## right side of the intron
    for each_nt_right in range(len(iter_seq_right)-args.indelSize):

        var_pos_right = each_nt_right + args.intronEnd-args.intronNts-1 ## right


        ref_del = iter_seq_right[each_nt_right:each_nt_right+args.indelSize+1]
        alt_del = iter_seq_right[each_nt_right]
        var_type_del = str(args.indelSize) + 'nt_del'
        out_right.write('\n'+args.chrom+'\t'+str(var_pos_right)+'\t'+ref_del+'\t'+alt_del+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_id_right +'_'+var_type_del+'\t'+args.orfID)
        var_index_right +=1
        var_id_right = 'XXX'+str(var_index_right)


        ref_ins = iter_seq_right[each_nt_right]
        all_alt_ins_right = generate_permutations(all_nts, args.indelSize)
        for each_right_alt in all_alt_ins_right:
            alt_ins = iter_seq_right[each_nt_right] + each_right_alt
            var_type_ins = str(args.indelSize) + 'nt_ins'

            out_right.write('\n'+args.chrom+'\t'+str(var_pos_right)+'\t'+ref_ins+'\t'+alt_ins+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_id_right +'_'+var_type_ins+'\t'+args.orfID)
            var_index_right +=1
            var_id_right = 'XXX'+str(var_index_right)


    out_left.close()
    out_right.close()




if __name__ == '__main__':
    main()