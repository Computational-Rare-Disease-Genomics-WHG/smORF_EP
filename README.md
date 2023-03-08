# smORF-EP
smORF-EP: predicting the effect of variants in small open reading frames


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
The scripts provided are using those, and pre-installation is required

# Installation

pip install smORF-EP


# Runnng smORF-EP


# Annotations description



# Sequence quality filters

- 3nt periodicity
- multiple stops
- last sequence trio is a STOP codon
- start/end not within an intron (from known transcripts)
