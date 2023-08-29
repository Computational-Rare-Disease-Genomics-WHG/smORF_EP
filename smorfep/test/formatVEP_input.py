#!/usr/bin/python3

## script to format the input files from the examples files generate with 
## the examplewizard from smORF-EP package

## VEP input format: 3 319780 . GA G . . .


import sys
from smorfep.utils.functions import *

inputname = sys.argv[1]
outputname = sys.argv[2]

file = read_variants_file(inputname, '\t', 0)

##print(file)

file.insert(2, 'id', '.')
file.insert(5, 'score', '.')
file.insert(6, 'other', '.')
file.insert(7, 'last', '.')

file = file.drop(columns=['start', 'end', 'strand', 'var_id', 'smorf_id'])
##print(file)

file.to_csv(outputname, index=False, header=False, sep='\t')