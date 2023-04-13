# smORF-EP

smORF-EP: predicting the effect of variants in small open reading frames.


- [Requirements](#requirements)
- [Install smORF-EP](#install-smorf-ep)
- [Download Ref and transcripts](#donwload-reference-and-transcripts)
- [Run smORF-EP](#run-smorf-ep)
- [Annotations description](#annotations-description)
- [Sequence quality filters](#sequence-quality-filtes)





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



# Download reference and GENCODE transcripts
<!--Follow the [instructions.](https://github.com/Computational-Rare-Disease-Genomics-WHG/smORF-EP/blob/main/data)-->


# Install smORF-EP

## Installation from repository

1. Download smORF-EP repository
```
git clone https://github.com/Computational-Rare-Disease-Genomics-WHG/smORF-EP/
```

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



# Run smORF-EP

Example:
```
python3 smORF-EP.py -r ref_genome/ -t transcriptCoord.tsv -i introns.tsv -f variants_file.tsv -o outputname.tsv
```
**-r**: reference genome repository, where the sequences per chromosome are stored.
**-t**: file with the transcripts coordiantes (start and end).
**-i**: file with the introns coordinates for the transcripts in the study.
**-f**: input file with the variants and respective smORF region (see Input(#input)).
**-o**: output file (see Output(#output))



# Annotations description


# Input and output description

## Input

## Output

