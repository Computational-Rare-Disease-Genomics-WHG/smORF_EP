#!/usr/bin/python3

## Contact: M Fernandes

## Script to generate examples of indels based on an interval and considering a sliding window

import sys
from smorfep.utils.functions import *
import random

### functions 
def choice_excluding(lst, exception):
    possible_choices = [item for item in lst if item != exception]
    return random.choice(possible_choices)


### --- main

print('python3 generate_examples.py <chr> <orf_start> <orf_end> <strand> <orf_id> <intron_start> <intron_end> <indel_max_size> <outputname>')
print('python3 generate_examples.py 18 47808957 47930659 - SMAD2 47865134 47868322 2 test_exon-intron_auto_gen.vcf')

chrom = sys.argv[1]
orf_start = sys.argv[2]
orf_end = sys.argv[3]
strand = sys.argv[4]
orf_id = sys.argv[5]
intron_start = int(sys.argv[6])
intron_end = int(sys.argv[7])
indel_max_size = int(sys.argv[8])
outputname = sys.argv[9]

exon_nts = 3
intron_nts = 8
all_nts = ['A', 'T', 'G', 'C']

## load the reeference genome chr sequence
ref_path = '/Users/mariaf/Desktop/GitHub/smORF_EP/ref_genome/'

files_prefix, files_suffix = check_prefix_sufix_ref_files(ref_path)
ref = read_single_fasta(str(chrom), ref_path, files_prefix, files_suffix)

## get sequence
seq_left = get_sequence(intron_start-exon_nts-1, intron_start+intron_nts-1+indel_max_size, strand, ref)
# print(intron_start, exon_nts, intron_nts)
# print(intron_start-exon_nts-1, intron_start+intron_nts-1)
##print(seq_left)

## output open and header
out = open(outputname, 'w')
out.write('chrm\tvar_pos\tref\talt\tstart\tend\tstrand\tvar_id\tsmorf_id')

## NOTE: starts 1 postion before the exons nts -- as we use anchor notation
for var_size in range(indel_max_size+1): ## +1 to include the indel_max_size 
    print(var_size)
    diff = indel_max_size-var_size
    iter_seq = seq_left[:len(seq_left)-diff]
    ##print(iter_seq)

    if var_size == 0: ## SNV -- runs only once
        for each_nt in range(len(iter_seq)-var_size):
            var_pos = each_nt + intron_start-exon_nts-1
            ref = iter_seq[each_nt]
            alt = choice_excluding(all_nts, ref)
            var_type = 'SNV'

            ##print(var_pos, ref, alt, var_type)

            ## uncoment next line for SNV
            ##out.write('\n'+chrom+'\t'+str(var_pos)+'\t'+ref+'\t'+alt+'\t'+orf_start+'\t'+orf_end+'\t'+strand+'\t'+var_type+'\t'+orf_id)



    else: ## insertion and deletion  
        for each_nt in range(len(iter_seq)-var_size):
            var_pos = var_pos = each_nt + intron_start-exon_nts-1

            ref_ins = iter_seq[each_nt]
            alt_ins = iter_seq[each_nt] + ''.join(random.choices(all_nts, k=var_size))
            var_type_ins = str(var_size) + 'nt_ins'

            ref_del = iter_seq[each_nt:each_nt+var_size+1]
            alt_del = iter_seq[each_nt]
            var_type_del = str(var_size) + 'nt_del'

            # print(var_pos, ref_ins, alt_ins, var_type_ins)
            # print(var_pos, ref_del, alt_del, var_type_del)
            out.write('\n'+chrom+'\t'+str(var_pos)+'\t'+ref_ins+'\t'+alt_ins+'\t'+orf_start+'\t'+orf_end+'\t'+strand+'\t'+var_type_ins+'\t'+orf_id)
            out.write('\n'+chrom+'\t'+str(var_pos)+'\t'+ref_del+'\t'+alt_del+'\t'+orf_start+'\t'+orf_end+'\t'+strand+'\t'+var_type_del+'\t'+orf_id)


out.close()

## TODO:
# seq_right = get_sequence(intron_end-intron_nts, intron_end+exon_nts+1, strand, ref)
# print(intron_end-intron_nts, intron_end+exon_nts+1)
# print(seq_right)