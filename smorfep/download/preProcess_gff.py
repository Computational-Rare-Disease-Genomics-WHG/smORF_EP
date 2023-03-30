#!/usr/bin/python3

## Author: M Fernandes

## Script to remove the head lines on the gff file (# lines )
## adds the fist line as the header -- Column names for the pandas dataframe (For next steps)

## Used for GENCODE 

import sys
import time


##print('python3 0-preProcess_gff.py <inputfilemane> <outputfilename>')
##print('Example: python3 0-preProcess_gff.py gencode/gencode.v41.annotation.gff3 gencode/gencode.v41.annotation_columnNames.gff3')

start_time = time.time()

infilename = sys.argv[1]
outputname = sys.argv[2]

extension = infilename.split('.')[-1]

infile = open(infilename,'r')
out = open(outputname, 'w')



## new header
# line has 9 columns
out.write('chr\tsource\ttype\tstart\tend\tscore\tstrand\tgen_phase\tinfo\n')

if extension == 'gff':
    for line in infile: 
        if line[0] != '#': ## removes header lines
            out.write(line)

elif extension == 'gff3':
    for line in infile: 
        if line[0] != '#' and 'GRCh38' not in line:  ## removes header lines and also the lines that contain just the chromosome info, and transcript separation '###\n'
            out.write(line)



infile.close()
out.close()

end_time = time.time() - start_time

print(end_time, ' seconds')

print('DONE')