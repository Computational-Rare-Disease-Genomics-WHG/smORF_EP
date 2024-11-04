#!/usr/bin/python3

## Script to get just the main sequence for each chromosome

## usage: ref_chr ref_genome/ filename outputname

import sys
import time


start_time = time.time()

## change filename and output for different names
##filename = 'GCF_000001405.39_GRCh38.p13_genomic.fna'
##output = '_GRCh38.p13_genomic.fna'

path = sys.argv[1]
filename = sys.argv[2]
output = sys.argv[3]


f = open(path + filename, 'r')
out = open(path + 'chr1' + output, 'w') ## used for the cycle to work

for line in f:
    if line[0] == '>': ## header line
        ##print(line)
        if 'scaffold' not in line and 'patch' not in line and 'mitochondrion' not in line:
            ## this gives just the main sequences headers
            ##print(line)
            h = True
            out.close() ## closes previous file
            c = line.split(',')[0].split(' ')[-1]
            out = open(path + 'chr'+ c + output, 'w')
            out.write(line)
        else:
            out.close()
            h = False
            
        ##print(h)
    else:
        if h == True:
            out.write(line)


f.close()
out.close()


print('DONE!')


end = (time.time() - start_time)/ 60.0

print(end, ' minutes.')

