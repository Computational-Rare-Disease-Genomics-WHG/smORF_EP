# smORF-EP

smORF-EP: predicting the effect of variants in small open reading frames.


- [How to cite](#how-to-cite)
- [Requirements](#requirements)
- [Install smORF-EP](#install-smorf-ep)
- [Download Reference and transcripts](#donwload-reference-and-transcripts)
- [Generate input](#generate-input)
- [Run smORF-EP](#run-smorf-ep)
- [Input and output](#intput-and-output)
- [Annotations description](#annotations-description)




# How to cite

<!-- 
'smORF-EP: predicting the effect of variants in small open reading frames'
-->

<br><hr>
[:arrow_up: Back to top](#smorf-ep)


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

<br><hr>
[:arrow_up: Back to top](#smorf-ep)


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

python3 setup.py install
```

<br><hr>
[:arrow_up: Back to top](#smorf-ep)

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
smorfinit --reference --ref_link https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109.20211119/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
```


## Download GENCODE transcripts only

Example (download GENCODE v41):

```
smorfinit --transcripts --transc_link https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gff3.gz
```

<br><hr>
[:arrow_up: Back to top](#smorf-ep)


# Generate input

```
smorfinput -r <path_to_ref_genome/>ref_genome/ -b <BED file> -v <VCF file> -o <outputname>
```

Add --bedheader if the first line in the BED file is the file header. 

Add --vcfheader if the first line in the VCF file is the file header.

Example for both input files with header: 
```
smorfinput -r <path_to_ref_genome/>ref_genome/ -b <BED file> --bedheader -v <VCF file> --vcfheader -o <outputname>
```

<br><hr>
[:arrow_up: Back to top](#smorf-ep)


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
Add path before input and output filename for specific locations.

<br><hr>
[:arrow_up: Back to top](#smorf-ep)


# Output

smorfep command generates three outputfiles:
- **consequence file**: File with the variants annotation per transcript
- **excluded file**: File listing the transcripts excluded per smORF and respective flag
- **no transcript file**: contains the smORF IDs for smORFs where no transcript was found, therefore the tool doesn't run on them.

## Output examples:


### Consequence file



### Excluded transcript file

The file maps per transcript ID the transcripts that are excluded from analysis and the respective exclusion criteria (type). Additionally, the last column reports the length of the sequence, used in case of exclusion by non-3nt periodicity. 

![alt text](https://github.com/Computational-Rare-Disease-Genomics-WHG/smORF_EP/blob/main/support/excluded_file_example.png?raw=true)

### No transcript file


<br><hr>
[:arrow_up: Back to top](#smorf-ep)


# Annotations description


<br><hr>
[:arrow_up: Back to top](#smorf-ep)
