# smORF-EP

smORF-EP: predicting the effect of variants in small open reading frames.


- [Requirements](#requirements)
  * [Python libraries used](#python-libraries-used)
- [Usage](#usage)
- [Annotations description](#annotations-description)
- [Squence quality filters](#sequence-quality-filtes)





# Requirements

## Python libraries used

- pandas (might require installation)
- re (might require installation)
- datetime
- os
- sys
- time 


## Files download and decompression

To obtain the reference genome and transcript information we used *wget* and for decompression we used *gzip*. 
The scripts provided are using those, and pre-installation is required.


# Downlaod human reference and GENCODE transcripts



# Running smORF-EP

Example:
```
python3 smORF-EP.py -r ref_genome/ -t transcriptCoord.tsv -i introns.tsv -f variants_file.tsv -o outputname.tsv
```
**-r**: reference genome repository, where the sequences per chromosome are stored.
**-t**: file with the transcripts coordiantes (start and end).
**-i**: file with the introns coordinates for the transcripts in the study.
**-f**: input file with the variants and respective smORF region (see Input(#input)).
**-o**: output file (see Output(#output))

<!-- # Installation

pip install smORF-EP -->



# Annotations description


# Input and output description

## Input

## Output


# Sequence quality filters

- 3nt periodicity
- multiple stops
- last sequence trio is a STOP codon
- start/end not within an intron (from known transcripts)
