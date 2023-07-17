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
    parser.add_argument('-ws','--intronStart', metavar='\b', required=True, type=int, help='intron start coordinate')
    parser.add_argument('-we','--intronEnd', metavar='\b', required=True, type=int, help='intro end coordinate')
    parser.add_argument('-m','--indelMaxSize', metavar='\b', required=True, type=int, help='indel maximum size')
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
    seq_left = get_sequence(args.intronStart-args.exonNts-1, args.intronStart+args.intronNts-1+args.indelMaxSize, '+', ref)
    # print(args.windowStart, args.exonNts, args.intronNts)
    # print(args.windowStart-args.exonNts-1, args.windowStart+args.intronNts-1)
    # print(seq_left)

    ## get sequence right end
    seq_right = get_sequence(args.intronEnd-args.intronNts-1, args.intronEnd+args.exonNts-1+args.indelMaxSize, '+', ref)

    ## output open and header
    out_left = open(args.outputfile+'_donor.vcf', 'w')
    out_left.write('chrm\tvar_pos\tref\talt\tstart\tend\tstrand\tvar_id\tsmorf_id')
    
    out_right = open(args.outputfile+'_acceptor.vcf', 'w')
    out_right.write('chrm\tvar_pos\tref\talt\tstart\tend\tstrand\tvar_id\tsmorf_id')

    var_index_left = 1
    var_id_left = 'XXX'+str(var_index_left)

    var_index_right = 1
    var_id_right = 'XXX'+str(var_index_right)


    ## NOTE: starts 1 postion before the exons nts -- as we use anchor notation
    for var_size in range(args.indelMaxSize+1): ## +1 to include the args.indelMaxSize
        print(var_size)
        diff = args.indelMaxSize-var_size
        iter_seq_left = seq_left[:len(seq_left)-diff] ## left side
        iter_seq_right = seq_right[:len(seq_right)-diff] ## right side
        ##print(iter_seq)

        if var_size == 0: ## SNV -- runs only once
            for each_nt_left in range(len(iter_seq_left)-var_size):
                var_pos_left = each_nt_left + args.intronStart-args.exonNts-1 ## left side
                ref = iter_seq_left[each_nt_left]
                alt = choice_excluding(all_nts, ref)
                var_type = 'SNV'

                ##print(var_pos_left, ref, alt, var_type)

                ## uncoment next line for SNV
                ##out_left.write('\n'+chrom+'\t'+str(var_pos_left)+'\t'+ref+'\t'+alt+'\t'+orf_start+'\t'+orf_end+'\t'+strand+'\t'+var_id_left +'\t'+var_type+'\t'+orf_id)
                ##var_index_left +=1
                ##var_id_left = 'XXX'+str(var_index_left)


            for each_nt_right in range(len(iter_seq_right)-var_size):
                var_pos_right = each_nt_right + args.intronEnd-args.intronNts-1 ## right side 
                ref = iter_seq_right[each_nt_right]
                alt = choice_excluding(all_nts, ref)
                var_type = 'SNV'

                ##out_left.write('\n'+chrom+'\t'+str(var_pos_right)+'\t'+ref+'\t'+alt+'\t'+orf_start+'\t'+orf_end+'\t'+strand+'\t'+var_id_right +'\t'+var_type+'\t'+orf_id)
                ##var_index_right +=1
                ##var_id_right = 'XXX'+str(var_index_right)


    ## done until here XXX

        else: ## insertion and deletion  

            ## left side of the intron
            for each_nt_left in range(len(iter_seq_left)-var_size):
                var_pos_left = each_nt_left + args.intronStart-args.exonNts-1 ## left side

                ref_ins = iter_seq_left[each_nt_left]
                alt_ins = iter_seq_left[each_nt_left] + ''.join(random.choices(all_nts, k=var_size))
                var_type_ins = str(var_size) + 'nt_ins'

                ref_del = iter_seq_left[each_nt_left:each_nt_left+var_size+1]
                alt_del = iter_seq_left[each_nt_left]
                var_type_del = str(var_size) + 'nt_del'


                # print(var_pos_left, ref_ins, alt_ins, var_type_ins)
                # print(var_pos_left, ref_del, alt_del, var_type_del)
                out_left.write('\n'+args.chrom+'\t'+str(var_pos_left)+'\t'+ref_ins+'\t'+alt_ins+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_id_left +'_'+var_type_ins+'\t'+args.orfID)
                var_index_left +=1
                var_id_left = 'XXX'+str(var_index_left)
                out_left.write('\n'+args.chrom+'\t'+str(var_pos_left)+'\t'+ref_del+'\t'+alt_del+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_id_left +'_'+var_type_del+'\t'+args.orfID)


                var_index_left +=1
                var_id_left = 'XXX'+str(var_index_left)



             ## right side of the intron
            for each_nt_right in range(len(iter_seq_right)-var_size):

                var_pos_right = each_nt_right + args.intronEnd-args.intronNts-1 ## right
                ref_ins = iter_seq_right[each_nt_right]
                alt_ins = iter_seq_right[each_nt_right] + ''.join(random.choices(all_nts, k=var_size))
                var_type_ins = str(var_size) + 'nt_ins'

                ref_del = iter_seq_right[each_nt_right:each_nt_right+var_size+1]
                alt_del = iter_seq_right[each_nt_right]
                var_type_del = str(var_size) + 'nt_del'


                out_right.write('\n'+args.chrom+'\t'+str(var_pos_right)+'\t'+ref_ins+'\t'+alt_ins+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_id_right +'_'+var_type_ins+'\t'+args.orfID)
                var_index_right +=1
                var_id_right = 'XXX'+str(var_index_right)
                out_right.write('\n'+args.chrom+'\t'+str(var_pos_right)+'\t'+ref_del+'\t'+alt_del+'\t'+args.orfStart+'\t'+args.orfEnd+'\t'+args.strand+'\t'+var_id_right +'_'+var_type_del+'\t'+args.orfID)

                var_index_right +=1
                var_id_right = 'XXX'+str(var_index_right)

    out_left.close()
    out_right.close()




if __name__ == '__main__':
    main()