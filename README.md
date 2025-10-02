
[![DOI](https://sandbox.zenodo.org/badge/605157988.svg)](https://handle.test.datacite.org/10.5072/zenodo.342187)


![logo](https://github.com/Computational-Rare-Disease-Genomics-WHG/smORF_EP/blob/main/support/smorfEPlogo_resized3.png) 
# smORF-EP

smORF-EP: predicting the effect of variants in small open reading frames.


- [smROF-EP Overview](#smorf-ep-overview)
- [How to cite](#how-to-cite)
- [Docker image](#docker-image)
- [Requirements](#requirements)
- [Install smORF-EP](#install-smorf-ep)
- [Download Reference and transcripts](#donwload-reference-and-transcripts)
- [Generate input](#generate-input)
- [Run smORF-EP](#run-smorf-ep)
- [Input and output](#intput-and-output)
- [Annotations description](#annotations-description)
- [Issues and new features](#issues-and-new-features)


# smORF-EP overview


<br><hr>
[:arrow_up: Back to top](#smorf-ep)


# How to cite

To cite smORF-EP please use the repository DOI: [![DOI](https://sandbox.zenodo.org/badge/605157988.svg)](https://handle.test.datacite.org/10.5072/zenodo.342187)



<br><hr>
[:arrow_up: Back to top](#smorf-ep)


# Docker image

The simplest way to install and run smORF-EP is using out Docker image. 

## Build Docker image
```
docker build -t "smorfepdev:Dockerfile" 
```
Where 'smorfepdev' is the docker image name.

## Check the Docker image built sucessfully
```
docker images
```
<br><hr>
[:arrow_up: Back to top](#smorf-ep)


# Requirements

## Python libraries

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


## Download reference genome only

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
smorfinput -b <BED file> --bedheader -v <VCF file> --vcfheader -o <outputname>
```

<br><hr>
[:arrow_up: Back to top](#smorf-ep)


# Run smORF-EP

Example:
```
smorfep -r <path_to_ref_genome/>ref_genome/ -t <path_to_transcriptCoord_file/transcriptCoord.tsv> -i <path_to_transcriptCoord_file/introns.tsv> -f <variants_file>.tsv -o <outputname.tsv>
```
**-r**: reference genome repository, where the sequences per chromosome are stored. <br />
**-t**: file with the transcripts coordiantes (start and end). <br />
**-i**: file with the introns coordinates for the transcripts in the study. <br />
**-f**: input file with the variants and respective smORF region (see Input(#input)). <br />
**-o**: output file (see Output(#output)). <br />


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

Annotations provided by smORF-EP follow the same standard of VEP.


| Annotation        | Description| 
| :---------------- | :--------- | 
| missense_variant   | Variant that changes one or more bases and results in a change in the amino acid sequence, however, the lentght is preserved  |
| inframe_insertion  | Inframe non synonymous variant that inserts bases  |
| inframe_deletion  | Inframe non synonymous variant that deletes bases  |
| frameshift_variant  | Variant which disrupts the translational reading frame, as the number of inserted or deleted nucleotides is not multiple of three  |
| stop_lost   | Variant changing at least one base of the terminatior codon (stop), leading to transcript elongation   |
| stop_gained   | Variant changing at least one base of the terminatior codon (stop), leading to transcript shortening   |
| stop_retained_variant  | Variant changing at least one base of the terminatior codon (stop), but the terminator remains  |
| start_lost   | Variant that affects the start codon and the reference sequence has a canonical start 'ATG'   |
| start_variant   | Variant that affects the start codon and the reference sequence has a non-canonical start   |
| start_retained_variant  | Variant changing at least one base of the start codon, but the start remains  |
| splice_donor_variant  | Variant that changes the 2 base region at the 5'end of the intron   |
| splice_donor_5th_base_variant  | Variant that changes the 5th base after the start of the intron  |
| splice_donor_region_variant   | Variant that changes between 3-6th bases after the splice junction on the 5'end    |
| splice_acceptor_variant     | Variant that changes the 2 base region at the 3'end of the intron  |
| splice_region_variant       | Variant that changes one of these regions: (i) 1-3 bases of the exon next to the splice donor/acceptor region, or (ii) 3-8 bases of the intron   |
| intron_variant  | Variant occurring within an intron  |
| splice_polypyrimidine_tract_variant   | Variant falling in the polypyrimidine tract at 3' end of intron between 17 and 3 bases from the end of the intron |
| synonymous_variant   | Variant that does not change the encoding amino acid  |

<!--
|   |   |
-->


<br><hr>
[:arrow_up: Back to top](#smorf-ep)

# Issues and new features

You can report any issue through the "Bug report" form or contact the authors.

There is also a form for "Feature request". 


<br><hr>
[:arrow_up: Back to top](#smorf-ep)

<!--

# Caveats
- Tested on 10 replicas for variants within exon-intron boundaries (forawrd and reverse strand, donor and acceptor sites)
-Teste on insertions crosisng exon-intron boundaries: all permutations of 6nts (forawrd and reverse strand, donor and acceptor sites)

-->
