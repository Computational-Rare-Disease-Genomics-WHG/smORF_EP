"""
Downloader script for smORF_EP
This script downloads the reference genome and the transcriptome from NCBI and Gencode respectively

Usage: smorfinit [OPTIONS]

"""


import os 
import sys
import urllib.request
import gzip
import datetime
import argparse 

def bedvcf2intput():
    """
    Compiles the input file read by smorfep.

    Input: 
    - bedfile: BED file with the smORFs start and end coordinates
    - vcffile: VCF file with the variants 
    """
    print("Generating smorfep input file...")
    

def main():
    """
    Main entry point
    """

    parser = argparse.ArgumentParser(description='Script to convert BED and VCF into smorfep input file')


    args = parser.parse_args()


    if args.reference: 
        download_ref_genome(args.ref_link)

    elif args.transcripts: 
        download_gencode(args.transc_link)

    elif args.all: 
        download_ref_genome(args.ref_link)
        download_gencode(args.transc_link)



if __name__ == '__main__':
    main()