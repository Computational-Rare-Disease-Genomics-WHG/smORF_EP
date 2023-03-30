"""
Downloader script for smORF_EP
This script downloads the reference genome and the transcriptome from NCBI and Gencode respectively

Usage: smorfinit [OPTIONS]

"""


import os 
import urllib.request
import gzip
import datetime

def download_ref_genome():
    """
    Downloads the reference genome from NCBI and splits it per chromosome
    """
    print("Downloading reference genome...")
    # 1. Creates a new dir for the reference genome
    os.makedirs('ref_genome', exist_ok=True)

    # 2. Download reference from NCBI to the ref_genome/repository -- For GRC38.p13
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109.20211119/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz'
    urllib.request.urlretrieve(url, 'ref_genome/GCF_000001405.39_GRCh38.p13_genomic.fna.gz')

    # 3. Uncompress reference genome and delete compressed version
    with gzip.open('ref_genome/GCF_000001405.39_GRCh38.p13_genomic.fna.gz', 'rb') as f_in:
        with open('ref_genome/GCF_000001405.39_GRCh38.p13_genomic.fna', 'wb') as f_out:
            f_out.write(f_in.read())
    os.remove('ref_genome/GCF_000001405.39_GRCh38.p13_genomic.fna.gz')

    # 4. Split reference per chromosome
    os.system('ref_per_chr.py GCF_000001405.39_GRCh38.p13_genomic.fna _GRCh38.p13_genomic.fna')
    # Excludes 'scaffold', 'patch', 'mitochondrion' sequences 
    # If others to exclude add condition to line 24 
    # outputname format: chr<num>_<suffix_user_given> (line above: chr<num>_GRCh38.p13_genomic.fna)

    print("Done")



def download_gencode():
    """
    Downloads Gencode
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
    url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gff3.gz'
    urllib.request.urlretrieve(url, 'transcripts/gencode.v41.annotation.gff3.gz')

    # 3- Uncompress reference genome and delete compressed version
    with gzip.open('transcripts/gencode.v41.annotation.gff3.gz', 'rb') as f_in:
        with open('transcripts/gencode.v41.annotation.gff3', 'wb') as f_out:
            f_out.write(f_in.read())
    os.remove('transcripts/gencode.v41.annotation.gff3.gz')

    # 4- Precomputations on genecode file
    # Pre-process gff3 file - single header
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
    download_ref_genome()
    download_gencode()

if __name__ == '__main__':
    main()