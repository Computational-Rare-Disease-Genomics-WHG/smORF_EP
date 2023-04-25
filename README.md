# smORF-EP

smORF-EP: predicting the effect of variants in small open reading frames.


- [Requirements](#requirements)
- [Install smORF-EP](#install-smorf-ep)
- [Download Reference and transcripts](#donwload-reference-and-transcripts)
- [Generate input](#generate-input)
- [Run smORF-EP](#run-smorf-ep)
- [Input and output](#intput-and-output)
- [Annotations description](#annotations-description)





# Requirements

## Python libraries used

requirements.txt file provides the dependencies for smORF-EP.


<!-- 
- pandas (might require installation)
- re (might require installation)
- datetime
- os
- sys
- time -->



# Install smORF-EP

## Installation from repository

1. Download smORF-EP repository
```
git clone https://github.com/Computational-Rare-Disease-Genomics-WHG/smORF-EP/
```
or

Through zip dowlnoad: 
![alt text](https://github.com/Computational-Rare-Disease-Genomics-WHG/smORF_EP/blob/main/support/download_git_repo.png?raw=true)


2. Move inside smORF-EP repository
```
cd smORF-EP
```

3. Install through pip
```
pip3 install .

python setup.py install
```


<!--
```

```
-->

<!--
Command when the package is in pip

pip install smORF-EP
-->


# Download reference and GENCODE transcripts
<!--Follow the [instructions.](https://github.com/Computational-Rare-Disease-Genomics-WHG/smORF-EP/blob/main/data)-->

## Download reference and GENCODE transcripts 

Example (download reference GCF_000001405.39_GRCh38.p13 and GENCODE v41): 

```
smorfinit --all --ref_link https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109.20211119/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz --transc_link https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gff3.gz
```


## Download refernce genome only

Example (download reference version GCF_000001405.39_GRCh38.p13): 

```
smorfinit -reference -ref_link https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109.20211119/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
```


## Download GENCODE transcripts only

Example (download GENCODE v41):

```
smorfinit --transcripts --transc_link https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gff3.gz
```


# Generate input

```
smorfinput -b <BED file> -v <VCF file> -o <outputname>
```

Add --bedheader if the first line in the BED file is the file header. 

Add --vcfheader if the first line in the VCF file is the file header.

Example for both input files with header: 
```
smorfinput -b <BED file> --bedheader -v <VCF file> --vcfheader -o <outputname>
```


# Run smORF-EP

Example:
```
smorfep -r <path_to_ref_genome/>ref_genome/ -t <path_to_transcriptCoord_file/transcriptCoord.tsv> -i <path_to_transcriptCoord_file/introns.tsv> -f <variants_file>.tsv -o <outputname.tsv>
```
**-r**: reference genome repository, where the sequences per chromosome are stored.
**-t**: file with the transcripts coordiantes (start and end).
**-i**: file with the introns coordinates for the transcripts in the study.
**-f**: input file with the variants and respective smORF region (see Input(#input)).
**-o**: output file (see Output(#output)). 

Note: 
Add path before input and output filename for scpecific locations.


# Input and output

## Input

## Output




# Annotations description


<br><hr>
[:arrow_up: Back to top](#smorf-ep)
