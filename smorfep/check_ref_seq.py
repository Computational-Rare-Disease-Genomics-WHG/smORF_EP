#!/usr/bin/python3

"""
Downloader script for smORF_EP
This script downloads the reference genome and the transcriptome from NCBI and Gencode respectively

Usage: checkseq [OPTIONS]

"""

import sys
from smorfep.utils.functions import get_sequence, read_single_fasta


def main():

    parser = argparse.ArgumentParser(description='Script to convert BED and VCF into smorfep input file')

    parser.add_argument('-r','--refpath', required=True, type=str, help='Path to the reference genome')
    parser.add_argument('-c', '--chromosome', required=True, type=str, help='chromosome')
    parser.add_argument('-s', '--start', required=True, type=str, help='start coordinate')
    parser.add_argument('-e', '--end', required=True, type=str, help='end coordinate')
    parser.add_argument('-p', '--strand', required=True, type=str, help='strand')

    args = parser.parse_args()


    ref = read_single_fasta(args.chromosome, args.refpath).upper() ## upper required as the sequence has capital and lower letters 

    seq = get_sequence(args.start, args.end, args.strand, ref)

    print(seq)


if __name__ == '__main__':
    main()
    

