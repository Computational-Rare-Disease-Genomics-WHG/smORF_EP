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

def download_ref_genome(ref_link):
    """
    Downloads the reference genome from NCBI and splits it per chromosome

    Input: 
    - ref_link: Full link for the reference file to be downloaded.
    """
    print("Downloading reference genome...")
    # 1. Creates a new dir for the reference genome
    # os.makedirs('ref_genome', exist_ok=True)

    # # 2. Download reference from NCBI to the ref_genome/repository -- For GRC38.p13
    # ##url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109.20211119/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz'
    
    # ## get the filename from path (last bit): 
    outputname = ref_link.split('/')[-1] 
    outputname_unconpress = outputname.strip('.gz')
    # urllib.request.urlretrieve(ref_link, 'ref_genome/' + outputname)

    # # 3. Uncompress reference genome and delete compressed version
    # with gzip.open('ref_genome/' + outputname, 'rb') as f_in:
    #     with open('ref_genome/'+ outputname_unconpress, 'wb') as f_out:
    #         f_out.write(f_in.read())
    # os.remove('ref_genome/' + outputname)

    # 4. Split reference per chromosome
    os.system('ref_per_chr.py ref_genome/ ' + outputname_unconpress + ' _GRCh38.p13_genomic.fna')
    # Excludes 'scaffold', 'patch', 'mitochondrion' sequences 
    # If others to exclude add condition to line 24 
    # outputname format: chr<num>_<suffix_user_given> (line above: chr<num>_GRCh38.p13_genomic.fna)

    print("Done")



def download_gencode(transc_link):
    """
    Downloads Gencode gff3 format (extensions used are .gff3 specific)

    Input: 
    -transc_link: Full link for the transcripts gff3 file (we used GENCODE)
    """
    print("Downloading transcripts...")
    # Script to obtain transcript information from GENCODE
    # Adapt to other sources

    # Uses urllib.request

    now = datetime.datetime.now().strftime('%Y-%m-%d')
    # print(now)

    # 1- Creates a new dir for the reference genome
    os.makedirs('transcripts', exist_ok=True)

    # 2- Download GENCODE
    ##url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gff3.gz'
    outputname = transc_link.split('/')[-1] 
    outputname_unconpress = outputname.strip('.gz')

    urllib.request.urlretrieve(transc_link, 'transcripts/'+ outputname)

    # 3- Uncompress reference genome and delete compressed version
    with gzip.open('transcripts/' + outputname, 'rb') as f_in:
        with open('transcripts/'+ outputname_unconpress, 'wb') as f_out:
            f_out.write(f_in.read())
    os.remove('transcripts/'+ outputname)

    # 4- Precomputations on genecode file
    # Pre-process gff3 file - single header ## TODO: make the names 
    prefix = outputname_unconpress.strip('.gff3')
    os.system('preProcess_gff.py transcripts/gencode.v41.annotation.gff3 transcripts/gencode.v41.annotation_columnNames.gff3')

    # transcript coordinates (15.seconds M1 ship)
    os.system(f'compute_transcripts_GENCODE.py transcripts/gencode.v41.annotation_columnNames.gff3 transcripts/gencode.v41.annotation_transcriptCoord_{now}.tsv')

    # intron coordinates
    os.system(f'compute_introns_GENCODE_perTransc.py transcripts/gencode.v41.annotation_columnNames.gff3 transcripts/gencode.v41.annotation_introns_{now}.tsv')

    print("Done")


def main():
    """
    Main entry point
    """

    ## TODO: Add arguments to download ref and transcripts to allow different versions
    parser = argparse.ArgumentParser(description='Script download the reference and/or the transcript information')
    ##Add the arguments as in the run.py script

    ## define arguments
    ## arguments mutually exclusive -- ro, to, all (reference only, trancripts only, both)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--reference', help='Run reference download', action='store_true')
    group.add_argument('--transcripts', help='Run transripts download', action='store_true')
    group.add_argument('--all', help='Run reference and transcripts download', action='store_true')


    ## define the mandatory other arguments deoending on the mutually excluding arguments
    if '-reference' in sys.argv: 
        parser.add_argument('--ref_link', required=True, type=str, help='reference genome link')

    elif '-transcripts' in sys.argv: 
        parser.add_argument('--transc_link', required=True, type=str, help='transcripts link')
        
    elif '-all' in sys.argv: 
        parser.add_argument('--ref_link', required=True, type=str, help='reference genome link')
        parser.add_argument('--transc_link', required=True, type=str, help='transcripts link')


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