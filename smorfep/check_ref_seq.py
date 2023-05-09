#!/usr/bin/python3

"""
Downloader script for smORF_EP
This script downloads the reference genome and the transcriptome from NCBI and Gencode respectively

Usage: checkseq [OPTIONS]

"""

import sys
import smorfep.utils.functions import get_sequence, read_single_fasta


def main():

    ref_path = sys.argv[1]
    chrom = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    strand = sys.argv[5]

    ref = read_single_fasta(chrom, ref_path).upper() ## upper required as the sequence has capital and lower letters 

    seq = get_sequence(start, end, strand, ref)

    print(seq)


if __name__ == '__main__':
    main()
    

