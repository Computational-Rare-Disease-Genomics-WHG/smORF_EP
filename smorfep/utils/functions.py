#!/usr/bin/python3

## functions used by main script smORF-EP

## VEP annotation matching to VEP v101 -- July 2023


from smorfep.utils.genetic_code import * 
from typing import List, Optional, Dict, Union
from itertools import product
import pandas as pd
import re
import os



def generate_permutations(characters, length):
    """
        Function to generate all the permutations of a given length, 
        given a set of possible characters.
        Input: 
        - characters: list with the possible values for permutations
        - lenght of the sequence

        Returns a list with all the possible permutations

        Note: Used for the examples all permutations for insertions case -- generate_examples.py
        We could use the fixed anchor function, or use this one and add the anchor before -- Used

    """
    return [''.join(p) for p in product(characters, repeat=length)]


# # Example usage:
# characters = ['A', 'T', 'G', 'C']
# length = 3

# permutations = generate_permutations(characters, length)
# for perm in permutations:
#     print(perm)


def generate_permutations_anchorFixed(characters, length, starting_char):
    """
        Function to generate the permutations of a given length, starting with a specific character  
        for a given a set of possible characters.

        Input:
        - characters: list with the possible values for permutations
        - lenght of the sequence
        - starting_char: used to filter and keep only permutations starting with a given character 

        Returns a list with all the possible permutations that start with starting_char

        Note: Useful if we want to consider a first fixed nucleotide (anchor nucleotide)

    """

    return [''.join(p) for p in product(characters, repeat=length) if p[0] == starting_char]

## Testing Block:
# characters = ['A', 'T', 'G', 'C']
# length = 3
# starting_letter = 'A'

# permutations = generate_permutations_anchorFixed(characters, length, starting_letter)
# for perm in permutations:
#     print(perm)


def read_hierarchy(filename) :
    """
        Function that reads a file and inport the content to a ordered list. 

        Returns the ordered list.
    """

    olist = []
    f = open(filename, 'r')

    for line in f: 
        line = line.strip('\n')
        olist.append(line)

    f.close()

    return olist


def high_var(list_flags, hierarchy):
    """
        Function to return the higher flag in a list, given a users priority list.

        Return the higher priority flag
    """

    for higher in hierarchy:
        if higher in list_flags:
            return higher


def check_prefix_sufix_ref_files(filepath):
    """
    Function to take the prefix and sufix of reference, 
    based on the files in the directory. 
    
    Assumes the per chromosome files are in the same format, 
    this means one common prefix and sufix. 

    Returns the file prefix and suffix.
    """

    files_list = os.listdir(filepath)
    ## takes first file
    f = files_list[0]
    ## find number in the filename -- variable part
    num = re.findall(r'\d+', f)[0]  ## [0] as the result is a list with a single element
    files_prefix = f.split(num)[0]
    files_suffix = f.split(num)[1]

    return files_prefix, files_suffix


def check_substring(substring, slist): 
    """
        Function to search a substring in a list of words.

        Returns the list of works where the substring is in. 
    """
    found_elements = [each for each in slist if substring in each]

    return found_elements


def read_single_fasta(chromosome, path, files_prefix, files_suffix):
    """
    Function to open fasta files.
    Input:
    - chromosome
    - path: where the reference sequence files are located
    - files_prefix: everything before the chromosome number on the filename
    - files_suffix: everything after the chromosome number on the filename

    Note: Focusing the study in a chromosome or a small set, the list option
        saves memory.

    Returns a dictionary with the chrmosome as key and the seqeunce as value.
    """

    fasta = '' ## dictionary with ref sequence per chromosome

    filename = path + files_prefix + chromosome + files_suffix
    f = open(filename, 'r')

    for line in f:
        if line[0] == '>':
            line = line.strip('\n')
            header = line
            seq = ''

        else:
            line = line.strip('\n')
            seq += line

    ## as each file has one chromosome, then at the end saves it in the dictionary
    fasta = seq.upper()
    ## upper required as the sequence has capital and lower letters 
    
    f.close()

    return fasta



def reverse_seq (seq):
    """
        Funtion that reads from the end of the sequence (inverts seqeunce)
    """
    return seq[::-1]



def complement_seq(seq):
    """
        Function to complement the sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    new_seq = ''.join([complement[base] for base in seq])
    
    return new_seq



def reverse_complement_seq(seq):
    """
        Function that applies both, reverse and complement of a DNA sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    new_seq = ''.join([complement[base] for base in seq[::-1]])
    
    return new_seq


def get_sequence(start, end, strand, ref):

    """
        Function to get the genomic sequence from the reference genome
        given the start position, the end position and the strand.

        It requires the specification of the location of the reference
        genome. 

        ## NOTE: 
        The chrososome does not need to be specified as it is assumed 
        that the reference is a single chromosome sequence. 

        Input: 
        - sequence start coordinate - ref genome
        - sequence end coordinate - ref genome
        - strand: for + strand it pulls the reference, for - strand it reverse complement the sequence

        Output: 
        genomci sequence for the given start, end and strand

    """

    ## .upper needed so the sequence is all capital letters. Ref contains lowercase letters 
    if strand == '+':
        seq = ref[start-1:end].upper() 

    elif strand == '-':
        seq = reverse_complement_seq(ref[start-1:end].upper())
    ## upper to have all capital letters needed for protein sequence
    
    return seq



def get_trios(seq):
    """
        Function that devides a sequence into trios.
        It can be addapted to other subsequence length.

        Input:
        - DNA sequence

        Returns a list with the trios --> used later for protein sequence check
    """

    n= 3 ## for trios
    trios_list = []

    for i in range(0, len(seq),n):
        trios_list.append(seq[i:i+n])

    return trios_list



def get_protein(sequence):
    """
        Giving a sequence converts a genomic sequence to proteing sequence.

        Input:
        - DNA sequence

        Returns the protein sequence
    """

    trios = get_trios(sequence)
    prot = ''

    for each in trios:
        amino = genetic_code_gen[each]

        if amino != 'STOP':    
           prot += amino
        elif amino == 'STOP':
            return prot

    return prot




## for sequences without introns
def add_variant(sequence, start, end,  ref, alt, position, strand):
    """
        Function to add a variant given:
        - sequence
        - start -- start genomic position of the sequence
        - end -- end genomic position of the sequence
        - ref -- reference allele
        - alt -- alternative allele
        - position -- genomic position of the variant
        - strand -- strand 

        Returns the new genomic sequence
    """

    if strand == '+':
        variant_index = position - start
    elif strand == '-':
        variant_index = end - position

    if sequence[variant_index:variant_index+len(ref)] == ref:
        prefix = sequence[:variant_index] ## excludes the variant position
        suffix = sequence[variant_index+1:] ## excludes the variant position

        if len(ref) > len(alt): ## deletions
            del_size = len(ref) - len(alt)
            new_alt = '-'* del_size

            ## for deletions the first postion on the ref is mantained, the remaining is removed
            ## so add 1 to prefix to keep this base
            prefix = sequence[:variant_index+1] 
            suffix = sequence[variant_index+len(ref):]

            new_seq = prefix + new_alt + suffix

        elif len(ref) <= len(alt): ## insertions and SNPs
            new_seq = prefix + alt + suffix
 
        
    else:
        print(position, ref, alt, strand, start, end)
        print('ref allele does not correspond!')
        print('reference genome: ', sequence[variant_index:variant_index+len(ref)])
        print('ref input: ', ref)
        print(position) 
        print(sequence)
        return None, sequence[variant_index:variant_index+len(ref)], ref

    return new_seq, None, None


## map genomic coordinates to transcript coordinates
def genome2transcript_coords(start, end, strand, introns_df):
    """
    Function that maps genomic coordinates and transcript coordinates.
    Takes into account introns if they exist. 

    - start: start coordinate (genomic)
    - end: end coordinate (genomic)
    - strand: strand (if reverse strand the counting starts from the end)
    - introns_df: df with introns coordinates, one intron per line

    Output: 
    (1)dictionary where keys are genomic coordinates and values the 
    corresponding transcript cooridnates.
    (2)dictionary where keys are transcript coordinates and values the 
    corresponding genomic cooridnates.
    """

    map_genome2transcript = {} ## dictionary genCoord --> transcCoord
    map_transcript2genome = {} ## dictionary transcCoord  --> genCoord

    if introns_df.empty: ## no introns from smORF start until transcript end
        if strand == '+': ## Checked - OK
            pos_transc = 0
            for each in range(start, end+1):  ## +1 to include the last position in the seq
                map_genome2transcript[each] = pos_transc  
                map_transcript2genome[pos_transc] = each

                pos_transc += 1

        elif strand == '-': ## Checked -- OK
            pos_transc = 0
            for i in range(end, start-1, -1):
                map_genome2transcript[i] = pos_transc
                map_transcript2genome[pos_transc] = i

                pos_transc += 1

    else: ## introns

        if strand == '+': ## checked - OK
            introns_df = introns_df.sort_values(by=['start']) ## introns in crescent order

            intron_num = 1
            pos_transc = 0 ## index starts in 0
            for index, row in introns_df.iterrows(): 
                if intron_num == 1: ## first intron
                    for val in range(start, row['start']): ## the actual position is row['start']-1 which is the postion before the intron start
                        ## we set row['start'] as range is end exclusive, so it stops at row['start']-1
                        map_genome2transcript[val] = pos_transc
                        map_transcript2genome[pos_transc] = val
                        
                        pos_transc += 1

                    next_start = row['end'] +1

                else: 
                    for v in range(next_start, row['start']):  ## Checked
                        map_genome2transcript[v] = pos_transc
                        map_transcript2genome[pos_transc] = v

                        pos_transc += 1

                    next_start = row['end'] +1

                intron_num += 1
            
            ## last bit
            for val_end in range(next_start, end+1): ## +1 to include the end position 
                map_genome2transcript[val_end] = pos_transc
                map_transcript2genome[pos_transc] = val_end
                pos_transc += 1

        elif strand == '-':
            introns_df = introns_df.sort_values(by=['start'], ascending=False) ## introns in crestcent order

            intron_num = 1
            pos_transc = 0 ## index starts in 0
            for index, row in introns_df.iterrows(): 
                if intron_num == 1: ## first intron

                    for val in range(end, row['end'],-1): ## genomic coord would be row['end']+1 -- python index starts at 0
                        map_genome2transcript[val] = pos_transc
                        map_transcript2genome[pos_transc] = val
                        
                        pos_transc += 1
                        
                    next_start = row['start'] -1

                else: 
                    for v in range(next_start, row['end'], -1): ## Checked - OK
                        map_genome2transcript[v] = pos_transc
                        map_transcript2genome[pos_transc] = v
                        
                        pos_transc += 1

                    next_start = row['start'] -1
                
                intron_num += 1
            
            ## last bit
            for val_end in range(next_start, start-1, -1): 
                map_genome2transcript[val_end] = pos_transc
                map_transcript2genome[pos_transc] = val_end
                pos_transc += 1


    # print(introns_df)
    # print(introns_df.start)
    # print(introns_df.end)
    # print(introns_df.transcript_id)
    # print(start, end)
    # print(map_transcript2genome)

    return map_genome2transcript, map_transcript2genome
    


## for sequences with introns
def add_variant_transcriptSeq(sequence, start, end, ref, alt, position, map_coordinates):
    """
        Function to add a variant given:
        - sequence
        - start -- start genomic position of the sequence
        - end -- end genomic position of the sequence
        - ref -- reference allele
        - alt -- alternative allele
        - position -- genomic position of the variant
        - strand -- strand 

        NOTE: It removes the introns before adding the variant. 
        Not used for variants that cross exon-intron sections.

        Returns the new genomic sequence
    """

    ## position in the sequence
    variant_index = map_coordinates[position]

    if sequence[variant_index:variant_index+len(ref)] == ref:
        prefix = sequence[:variant_index] ## excludes the variant position
        suffix = sequence[variant_index+1:] ## excludes the variant position

        if len(ref) > len(alt): ## deletions
            del_size = len(ref) - len(alt)
            new_alt = '-'* del_size

            prefix = sequence[:variant_index+1] 
            suffix = sequence[variant_index+len(ref):]

            new_seq = prefix + new_alt + suffix

        elif len(ref) <= len(alt): ## insertions and SNPs
            new_seq = prefix + alt + suffix
 
        
    else:
        print('add_variant_transcriptSeq function')
        print(position, ref, alt, start, end)
        print('ref allele does not correspond!')
        print('reference genome: ', sequence[variant_index:variant_index+len(ref)])
        print('ref input: ', ref) 
        ##print(map_coordinates)
        ##print(position) 
        print(sequence)
        return None, sequence[variant_index:variant_index+len(ref)], ref

    return new_seq, None, None






def nt2aaMAP(sequence): 
    """
        This function creates the index map from nucleotides sequence to aminoacid seqeunce. 
        So, nucleotides 1,2,3 have corresponds to the first aminoacid so all have index 0.
    """
    nt_aa = {} 

    index_aa = 0

    for i in range(0,len(sequence)):
        real_nt = i+1 ## because index starts at 0 -- nt 3 will be index 2, nt 6 -> index 5, ...
        
        if real_nt % 3 == 0 and i != 0: ## per 3 nts we pass for next aminoacid  
            nt_aa[i] = index_aa
            index_aa += 1 ## next aa
        else:
            nt_aa[i] = index_aa

    return nt_aa


def read_file(filename, delimiter, header):
    """ Function to open a file with a specified delimiter and return the
        information  it constains.

        Input:
        - filename
        - delimiter: '\t', ',', ';', ':', '.'
        - header: 0 or None ## 1 for the first line as header
        
        Returns the information on the file in a dictionary with the
        first column as key, and a columns_index list with the info
        for each record in the file.

        Usage example: openFile(<filename>, '\t', 1) ## first line is the header
    """
    
    df = pd.read_csv(filename, sep=delimiter, header=header, lineterminator='\n')

    return df



def read_variants_file(filename, delimiter, header):
    """ Function to open a file with a specified delimiter and return the
        information  it constains.

        Input:
        - filename
        - delimiter: '\t', ',', ';', ':', '.'
        - header: 0 or None ## 1 for the first line as header
        
        Returns the information on the file in a dictionary with the
        first column as key, and a columns_index list with the info
        for each record in the file.

        Usage example: openFile(<filename>, '\t', 1) ## first line is the header

    """
    
    df = pd.read_csv(filename, sep=delimiter, header=header)

    return df


def read_vep_annotations(filename):
    """
        Function to read and format the VEP annotations into a pandas dataframe
        for further processing. 

        NOTE: Columns format used here match the output from running VEP on the HPC 

        Returns the pandas dataframe
    """

    vep_df = pd.read_csv(filename, sep='\t', header=None, comment='#')
    vep_df.columns = ['var_id', 'location',	'Allele', 'Gene', 'Feature', \
        'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', \
        'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', \
        'IMPACT', 'DISTANCE', 'STRAND', 'FLAGS', 'existing_InFrame_oORFs', \
        'existing_OutOfFrame_oORFs', 'existing_uORFs', 'five_prime_UTR_variant_annotation', 'five_prime_UTR_variant_consequence']

    ## split location column to chromosome and postion
    vep_df[['chrm', 'var_pos']] = vep_df['location'].str.split(':',expand=True)

    ## split var_pos and keep first only -- some entries have more than one position
    vep_df['var_pos'] = vep_df['var_pos'].str.split('-').str[0]

    vep_df = vep_df.drop(columns=['location', 'Existing_variation', 'DISTANCE', 'FLAGS', 'existing_InFrame_oORFs', 'existing_OutOfFrame_oORFs', 'existing_uORFs'])

    ## re-order columns in the df
    vep_df = vep_df[['var_id','chrm', 'var_pos', 'Allele', 'STRAND', 'Gene', 'Feature', \
        'Feature_type', 'Consequence', 'Protein_position', 'Amino_acids', 'Codons', \
        'IMPACT', 'five_prime_UTR_variant_annotation', 'five_prime_UTR_variant_consequence']]

    ## replace nomenclature 1/-1 by +/- for strand
    vep_df['STRAND'] = vep_df['STRAND'].replace(1,'+')
    vep_df['STRAND'] = vep_df['STRAND'].replace(-1,'-')


    return vep_df



def read_vep_web(filename):
    """
        Function to read and format the VEP annotations into a pandas dataframe
        for further processing. 

        NOTE: Columns format from file obtained using VEP web API

        Returns the pandas dataframe
    """

    vep_df = pd.read_csv(filename, sep ='\t', header=None, comment='#')
    vep_df.columns = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']

    ## Create df with a consequence per transcript
    vep_transc_annot_df = pd.DataFrame(data=None, columns=['CHROM','POS','ID','REF','ALT', 'TRANSC_ID', 'SO','IMPACT'])

    index_gen = 0

    for index, row in vep_df.iterrows():
        c = row.CHROM
        pos = row.POS
        var_id = row.ID
        r = row.REF
        a = row.ALT

        info = row.INFO.split(',')
        for each_entry in info: 
            each_entry = each_entry.split('|')

            transc_id = each_entry[6]
            so = each_entry[1]
            impact = each_entry[2]

            new_line = pd.DataFrame(
                {
                'CHROM': c,
                'POS' : pos,
                'ID' : var_id,
                'REF' : r,
                'ALT' : a,
                'TRANSC_ID' : transc_id,
                'SO' : so,
                'IMPACT' : impact
                }, index=[index_gen]
            )

            vep_transc_annot_df = pd.concat([vep_transc_annot_df, new_line])

    return vep_transc_annot_df
    


def protein_consequence(seq, new_seq, var_pos, start, end, strand):
    """
        Function that checks the difference between protein sequences, 
        given two DNA sequences: initial and after mutation

        Input: 
        - seq = sequence without mutation
        - new_seq = sequence with mutation

    """

    ##  Get protein sequence from Ref_sequence
    prot_seq = get_protein(seq) 

    ##  New protein sequence
    new_prot = get_protein(new_seq)

    ## synonymous variant 
    if prot_seq == new_prot:
        consequence = 'synonymous_variant'
        change = 0
        
        return consequence, change

        

    elif len(prot_seq) == len(new_prot): ## single aminoacid change 

        consequence = 'missense_variant'
        ## since is a single nucleotide that changes, sets diff can be computed

        ## For the aminoacid change report - change
        nt_aa = nt2aaMAP(seq)

        if strand == '+': 
            change_index =(var_pos-start) 
        elif strand == '-':
            change_index = (end - var_pos) 
        aa_index = nt_aa[change_index]

        ref_aa = prot_seq[aa_index]
        alt_aa = new_prot[aa_index]

        change = ref_aa + '>' + alt_aa

        return consequence, change

    else: ## insertions/delitions -- The protein sequence is fixed, so frameshift does not apply at this level
    
        change = len(new_prot) - len(prot_seq)
        if change > 0: 
            consequence = 'protein_elongation'
        elif change < 0: 
            consequence = 'protein_truncation'

        return consequence, change



def protein_consequence_transcript(seq, new_seq, var_pos, map_coordinates):
    """
        Function that checks the difference between protein sequences, 
        given two DNA sequences: initial and after mutation

        Input: 
        - seq = sequence without mutation
        - new_seq = sequence with mutation
        - var_pos = variant position
        - map_coordinates = dictionary with genome to transcript coordinates mapping

        Output: function returns the consequence of the variant at the protein level

    """

    ##  Get protein sequence from Ref_sequence
    prot_seq = get_protein(seq) 

    ##  New protein sequence
    new_prot = get_protein(new_seq)

    ## synonymous variant 
    if prot_seq == new_prot:
        consequence = 'synonymous_variant'
        change = 0
        
        return consequence, change

        

    elif len(prot_seq) == len(new_prot): ## single aminoacid change 
        consequence = 'missense_variant'
        ## since is a single nucleotide that changes, sets diff can be computed

        ## For the aminoacid change report - change
        nt_aa = nt2aaMAP(seq)

        change_index =map_coordinates[var_pos] 

        aa_index = nt_aa[change_index]

        ref_aa = prot_seq[aa_index]
        alt_aa = new_prot[aa_index]

        change = ref_aa + '>' + alt_aa

        return consequence, change

    else: ## insertions/deletions -- The protein sequence is fixed, so frameshift does not apply at this level
        change = len(new_prot) - len(prot_seq)
        if change > 0: 
            consequence = 'protein_elongation'
        elif change < 0: 
            consequence = 'protein_truncation'

        return consequence, change




def frameshift(seq, transcript_extension, map_coordinates): 
    """
        Function that extends the sequence and searches for the new stop codon. 
        Used to access the consequence on the protein when a frameshift mutation happens.

        Works on genomic sequence: ATGC

        Input: 
        - seq = sequence after mutation
        - transcript_extension = Seqeunce within which we search for a new stop codon - if off this 
            range we just flag as out of transcript bound
        - map_coordinates = dictionary with the correspondence of transcript to genomic coordinates 

        Note: as we provide the searcheable range we don't need strand neither coordinates 
            as in the previous version of this function


        Output: 
        - new_seq = sequence until the first stop codon found, given the frameshift
    """

    seq_len = len(seq)
    difference = seq_len%3 ## -> how many nucleotides it shifts: 0 (inframe), 1 or 2
    print(difference)

    if difference == 1: ## we need to add 2 base to get new frame
        corrected_seq = seq + transcript_extension[:2]
        transcript_extension = transcript_extension[2:] ## removes the 2 added characters from the extension

    elif difference == 2: ## we need to add 1 base to get new frame 
        corrected_seq = seq + transcript_extension[:1]
        transcript_extension = transcript_extension[1:]

    new_seq, new_stop = stop_transcript_search(corrected_seq, transcript_extension, map_coordinates)

    if new_stop == None: ## stop not found
        return None

    else: ## stop found (1-on the corrected sequence 2-on the extended sequence (until end of the transcript))
        return new_seq



def find_stop_inframe(seq, map_coordinates):
    """
        Function to find inframe stop codons.
        Used to find regions of interest with a start and end codon.

        Works on genomic coordinates and do not take into account introns.

        Input: seq - sequence to be search
               map_coordinates - mapping between transcript and genomic coordinates, used to report the genomic 
                        coordinate for the new stop codon inframe

        output: Corrdinate of the last nucleotide of the new stop codon, and index of the new stop codon in the sequence.
                Or 'None, None' if any inframe stop codon is found. 

    """

    seq_trios = get_trios(seq)

    stop_trios = []
    if 'TAG' in seq_trios:
        indices = [i for i, x in enumerate(seq_trios) if x == 'TAG']
        stop_trios.extend(indices)
    if 'TAA' in seq_trios:
        indices = [i for i, x in enumerate(seq_trios) if x == 'TAA']
        stop_trios.extend(indices)
    if 'TGA' in seq_trios:
        indices = [i for i, x in enumerate(seq_trios) if x == 'TGA']
        stop_trios.extend(indices)
    
    if stop_trios == []: ## no stop  in frame
        return None, None
    else: 
        ## as we work in the sequence we want alway the first stop position -> first seq_index on the stops list 
        stop_trios.sort()
        new_stop_index = stop_trios[0]*3 ## gives the index of the first letter on the stop codon

        ##print(map_coordinates)
        ##print(new_stop_index, new_stop_index+3)

        return map_coordinates[new_stop_index+3], new_stop_index+3
    
    
def find_stop_inframe_sequence(seq):
    """
        Function to find inframe stop codons.
        Used to check if an insertion adds a stop codon. 
        check_introns function used this function.

        Works on the input sequence only, computes trios checks if there is a 
        stop codon in the tios list. 

        NOTE: Does not matter if the last trio has missing nts -- will not influence 
        the search as it will not match the stop codons. No need to remove it.

        Input: seq - sequence to be search

        output: True if there is a stop codon in the sequence  

    """
    
    seq_trios = get_trios(seq)

    stop_trios = []
    if 'TAG' in seq_trios:
        indices = [i for i, x in enumerate(seq_trios) if x == 'TAG']
        stop_trios.extend(indices)
    if 'TAA' in seq_trios:
        indices = [i for i, x in enumerate(seq_trios) if x == 'TAA']
        stop_trios.extend(indices)
    if 'TGA' in seq_trios:
        indices = [i for i, x in enumerate(seq_trios) if x == 'TGA']
        stop_trios.extend(indices)
    
    if stop_trios == []: ## no stop  in frame
        found_stop = False
    else: 
        ## as we work in the sequence we want alway the first stop position -> first seq_index on the stops list 
        stop_trios.sort()
        found_stop = True

    return found_stop
    

def stop_transcript_search(seq, transcript_extension, map_coordinates):
    """
        Function to search a new stop when a stop lost variant occurs. 

        Searches within the sequence until the end of a transcript. 

        Takes into account the introns, as we give the correspondence between transcript and genomic coordinates.

        Input:
        - seq: 
        - start: 
        - strand: 
        - transcript_extension: 
        - map_coordinates:

        Returns the new sequence and the new stop coordinate
    """

    ## 1- first search in the corrected sequence 
    new_stop, new_stop_index = find_stop_inframe(seq, map_coordinates)

    if new_stop != None: 
        new_seq = seq[:new_stop_index]

    
    else: 
        ## 2 - if not in the correected sequence, search until the end of the transcript
        extended_sequence = seq + transcript_extension
        new_stop, new_stop_index = find_stop_inframe(extended_sequence, map_coordinates)  ## genomic coordinate of the new stop (including stop codon)

        if new_stop == None:
            ## stop not within the transcript range
            return None, None


        else: ## stop found within the transcript
            ## get sequence index of the new stop
            new_seq = extended_sequence[:new_stop_index]
            ## independent of the strand, as we work on the transcript sequence start --> stop

    return new_seq, new_stop





def remove_introns(introns_df, start, end, strand, ref):
    """
    Function that removes introns from a sequence and also output the new length.

    Input: 
    - introns_df: dataframe with introns coordinates, one per line
    - start: start coordinate of the sequence in study
    - end: end coordinate of the sequence in study
    - strand: which strand we are taking
    - ref: reference genome sequence

    Output: 
    sequence without introns and new length.
    """

    seq_without_introns = '' 
    
    introns_checked = 0

    if strand == '+': ## sort crescent order
        introns_df = introns_df.sort_values(by='start', ascending=True)
    elif strand == '-': 
        introns_df = introns_df.sort_values(by='start', ascending=False)

    for index, row in introns_df.iterrows(): ## compute the sum of all introns len

        if introns_checked == 0: ## first exon
            ## start is the region start (start)
            if strand == '+':
                exon_end = row['start'] - 1
                seq_without_introns = get_sequence(start, exon_end, strand, ref)
                next_start = row['end'] + 1
            elif strand == '-':
                exon_end = row['end'] + 1
                seq_without_introns = get_sequence(exon_end, end, strand, ref)
                next_end = row['start'] - 1
            ## first bit of the sequence
            
        elif introns_checked != introns_df.shape[0]: ## not the last line - exons in the middle

            if strand == '+': ## add sequence after
                exon_start = next_start
                exon_end = row['start'] - 1
                seq_without_introns = seq_without_introns + get_sequence(exon_start, exon_end, strand, ref)
                next_start = row['end'] + 1
            
            elif strand == '-': ## add sequence before
                exon_start = row['end'] + 1
                exon_end = next_end
                seq_without_introns = seq_without_introns + get_sequence(exon_start, exon_end, strand, ref)
                next_end = row['start'] - 1
 
        
        introns_checked += 1 
        
    ## last exon
    ## end is the end of the region
    if strand == '+':
        exon_start = next_start 
        seq_without_introns = seq_without_introns + get_sequence(exon_start, end, strand, ref)
    elif strand == '-':
        exon_end = next_end
        seq_without_introns = seq_without_introns + get_sequence(start, exon_end, strand, ref)

        
    new_len = len(seq_without_introns)

    
    return seq_without_introns, new_len



def check_start(seq, new_sequence, start, end, variant_pos, strand): 
    """
        Function to check if the variant affects the start codon.
        Input: 
        - sequence: before variant
        - new_sequence: sequence with the variant
        - start: start coordinate of the sequence (genomic coord)
        - end: end coordinate of the sequence (genomic coord)
        - variant_pos: variant position (genomic coord)
        - strand

        Note: Start retain is would be reported as None and other type of effect as frameshift will be checked

        Returns: start_lost if the variant affects the start codon, or None otherwise.

    """


    if strand == '+':
        len_change = '-'
        prot_cons = '-'
        change_prot = '-'

        if variant_pos >= start and variant_pos <= start +2 and seq[:3] != new_sequence[:3]: 
        ## last condition is for indels right after the start codon, as the coordinate will be the last nt of the start codon, which is kept the same
            len_change = seq[:3] + '->' + new_sequence[:3]
            return 'start_lost', len_change, prot_cons, change_prot
        
        else:
            return None, len_change, prot_cons, change_prot

    elif strand == '-': 
        len_change = '-'
        prot_cons = '-'
        change_prot = '-'

        if variant_pos <= end and variant_pos >= end -2 and seq[:3] != new_sequence[:3]: ## last condition is for indels right after the start codon, as the coordinate will be the last nt of the start codon, which is kept the same
            len_change = seq[:3] + '->' + new_sequence[:3]
            return 'start_lost', len_change, prot_cons, change_prot
    
        else: 

            return None, len_change, prot_cons, change_prot




def check_start_transcript(seq, new_sequence, variant_pos, map_coordinates): 
    """
        Function to check if the variant affects the start codon.
        Input: 
        - sequence: before variant
        - new_sequence: sequence with the variant
        - variant_pos: variant position (genomic coord)
        - map_coordinates: dictionary with the genomic to transcript positions mapping
   

        Returns: start_lost if the variant affects the start codon, or None otherwise.

    """

    transcript_var_pos = map_coordinates[variant_pos]

    len_change = '-'
    prot_cons = '-'
    change_prot = '-'

    ## independent of the strand, as the trancript is always from start codon to stop codon
    ## mapping is addapted: transcript starts from start if forward strand, and end for the reverse strand
    if transcript_var_pos in [0,1,2] and seq[:3] != new_sequence[:3]: 
        if new_sequence[:3] == 'ATG': ## not considering change to non-canonical start codon.
            len_change = seq[:3] + '->' + new_sequence[:3]
            return 'start_retained_variant', len_change, prot_cons, change_prot
        else: 
            len_change = seq[:3] + '->' + new_sequence[:3]
            return 'start_lost', len_change, prot_cons, change_prot

    else: ## not a start affecting variant
        return None, len_change, prot_cons, change_prot




def check_stop(seq, new_sequence, start, end, variant_pos, strand, transcript_info, ref_sequence, map_coordinates):
    """
        Function to check if the variant affects the stop codon.
        Input: 
        - sequence: before variant
        - new_sequence: sequence with the variant
        - start: start coordinate of the sequence (genomic coord)
        - end: end coordinate of the sequence (genomic coord)
        - variant_pos: variant position (genomic coord)
        - strand
        - transcript_info: used to check the transcript end coordinate
        - ref_sequence need in case the stop is lost and there is an extension
        - map_cordinates: mapping between transcript and genomic coordinates

        Returns: stop_lost or stop_retain if the variant affects the stop codon, or None otherwise.

    """

    stop_codons = ['TAG', 'TAA', 'TGA']

    if strand == '+': 
        if variant_pos <= end and variant_pos >= end -2: 
            ## new_sequence is the sequence with the variant
            ## new end is not a stop codon

            if end == transcript_info.iloc[0].end: ## smorf stops in the last position of the transcript
                if len(seq) == len(new_sequence):
                    len_change = seq[len(seq)-3:] + '->' + new_sequence[len(new_sequence)-3:]

                    if new_sequence[len(new_sequence)-3:] in stop_codons:
                        return 'stop_retained_variant', len_change, '-', '-'
                else:
                    return None, 'off_transcript_stop', '-', '-'

            else:
                seq2transcEnd = get_sequence(end+1, transcript_info.iloc[0].end, transcript_info.iloc[0].strand, ref_sequence)

                new_seq, new_stop = stop_transcript_search(new_sequence, seq2transcEnd, map_coordinates)

                if new_seq == None: 
                    len_change = 'off_transcript_stop'
                    prot_cons = '-'
                    change_prot = '-'

                else: 
                    len_change = len(new_seq) - len(seq)

                    if len_change == 0: 
                        len_change = seq[len(seq)-3:] + '->' + new_seq[len(new_seq)-3:]
                        return 'stop_retained_variant', len_change, '-', '-'
                    else:
                        prot_cons, change_prot = protein_consequence(seq, new_seq, variant_pos, start, end, strand)
                
                return 'stop_lost', len_change, prot_cons, change_prot

        return None, '-', '-', '-'
        
    if strand == '-':

        if variant_pos >= start and variant_pos <= start +2: 

            if start == transcript_info.iloc[0].start: ## smorf stops in the last position of the transcript
                if len(seq) == len(new_sequence):
                    len_change = seq[len(seq)-3:] + '->' + new_sequence[len(new_sequence)-3:]
                    if new_sequence[len(new_sequence)-3:] in stop_codons:
                        return 'stop_retained_variant', len_change, '-', '-'
                else:
                    return None, 'off_transcript_stop', '-', '-'

            else:
                seq2transcEnd = get_sequence(transcript_info.iloc[0].start, start, transcript_info.iloc[0].strand, ref_sequence)

                new_seq, new_stop = stop_transcript_search(new_sequence, seq2transcEnd, map_coordinates)
                
                if new_seq == None: 
                    len_change = 'off_transcript_stop'
                    prot_cons = '-'
                    change_prot = '-'

                else: 
                    len_change = len(new_seq) - len(seq)

                    if len_change == 0: 
                        len_change = seq[len(seq)-3:] + '->' + new_seq[len(new_seq)-3:]

                        return 'stop_retained_variant', len_change, '-', '-'

                    else:
                        prot_cons, change_prot = protein_consequence(seq, new_seq, variant_pos, start, end, strand)
            
                    return 'stop_lost', len_change, prot_cons, change_prot

        return None, '-', '-', '-'



def check_stop_transcript(seq, new_sequence, start, end, variant_pos, strand, map_coordinates_g, map_coordinates_t, extension_seq):
    """
        Function to check if the variant affects the stop codon.
        
        Input: 
        - sequence: before variant
        - new_sequence: sequence with the variant
        - start: start coordinate of the sequence (genomic coord)
        - end: end coordinate of the sequence (genomic coord)
        - variant_pos: variant position (genomic coord)
        - strand
        - map_coordinates_g: dictionary with genomic to transcript coordinates mapping
        - map_coordinates_t: dictionary with transcript to genomic coordinates mapping
        - ref_sequence need in case the stop is lost and there is an extension

        Returns: stop_lost or stop_retain if the variant affects the start codon, or None otherwise.

    """

    stop_codons = ['TAG', 'TAA', 'TGA']

    ## need this because the mapping is made until the end of the transcript, not the end of the region of interest
    if strand == '+':
        last_position = map_coordinates_g[end]
        stop_positions = [map_coordinates_t[last_position-2], map_coordinates_t[last_position-1], map_coordinates_t[last_position]]
        ## below does not work if there are introns in the stop codon
        stop_positions = [end-2, end-1, end]
    elif strand == '-':
        last_position = map_coordinates_g[start]

        stop_positions = [map_coordinates_t[last_position-2], map_coordinates_t[last_position-1], map_coordinates_t[last_position]]
        ## below does not work if there are introns in the stop codon
        stop_positions = [start+2, start+1, start]

    ## mapping is addapted: transcript starts from start if forward strand, and end for the reverse strand
    if variant_pos in stop_positions: 
        ## new_sequence is the sequence with the variant
        ## new end is not a stop codon
        new_seq, new_stop = stop_transcript_search(new_sequence, extension_seq, map_coordinates_t)

        if new_seq == None: 
            len_change = 'off_transcript_stop'
            prot_cons = '-'
            change_prot = '-'

        else: 
            len_change = len(new_seq) - len(seq)

            if len_change == 0: 
                len_change = seq[len(seq)-3:] + '->' + new_seq[len(new_seq)-3:]

                return 'stop_retained_variant', len_change, '-', '-'
            
            else: 
                prot_cons, change_prot = protein_consequence(seq, new_seq, variant_pos, start, end, strand)

            return 'stop_lost', len_change, prot_cons, change_prot
    
    return None, '-', '-', '-'
        

def compatibility_smorf_transcript(ref_sequence, transcript_info, introns_df, smorf_id, smorf_start, smorf_end, strand):
    """
    Function to check the compatibility of smORF and transcript coordinates. 

    Compares one smorf with all the transcripts it falls completely within.

    Input:
    - ref_sequence:  reference sequence for the chromosome in analysis
    - transcript_info: dataframe structure line with the information about the transcript (start, end, id, etc)
    - introns_df: intronic regions dataframe for the chromosome in analysis
    - smorf_id (for final output format)
    - smorf_start
    - smorfs_end
    - strand 

    returns a list of matching transcripts, a dataframe with unmatching transcript and the respective reason for exclusion
    """

    matching_transcripts = []
    unmatching_trancripts = pd.DataFrame(columns=['smorf_id','transcript_id', 'type', 'length'])

    transcripts_mapping_dictionary = {} ## per transcript ID stores map_gen2transc, map_transc2gen

    stop_codons = ['TAG', 'TAA', 'TGA']

    for index, row in transcript_info.iterrows(): ## transcript coordinates, each line assumed to be a different transcript

        t_id = row.transcript_id
        t_start = row.start
        t_end = row.end
        t_strand = row.strand 
        ##print(t_id)

        ## introns for the transcript
        introns_transcript = introns_df[introns_df['transcript_id'] == t_id]
        ## from this point on we use introns_transcript to filter to work on the transcript introns only
        ##print(introns_transcript)

        if strand == '+':
            introns_transcript = introns_transcript[(introns_transcript['start']>= smorf_start) & (introns_transcript['end']<= t_end)]

        elif strand == '-':
            introns_transcript = introns_transcript[(introns_transcript['start']>= t_start) & (introns_transcript['end']<= smorf_end)]

        ## smorf without introns
        if introns_transcript.empty:
            smorf_seq = get_sequence(smorf_start, smorf_end, strand, ref_sequence)
            smorf_len = len(smorf_seq)

            ## compute coordinates map 
            if strand == '+':
                map_gen2transc, map_transc2gen = genome2transcript_coords(smorf_start, t_end, strand, introns_transcript)
            elif strand == '-':
                map_gen2transc, map_transc2gen = genome2transcript_coords(t_start, smorf_end, strand, introns_transcript)


            ## Checks the number of stop codons in the sequence
            s, s_index = find_stop_inframe(smorf_seq[:len(smorf_seq)-3], map_transc2gen) ## removes last codon and searches for stop codons inframe
                

        ## smorf with introns
        elif not introns_transcript.empty:
            smorf_seq, smorf_len = remove_introns(introns_transcript, smorf_start, smorf_end, strand, ref_sequence)

            ## compute coordinates map 
            if strand == '+':
                map_gen2transc, map_transc2gen = genome2transcript_coords(smorf_start, t_end, strand, introns_transcript)
            elif strand == '-':
                map_gen2transc, map_transc2gen = genome2transcript_coords(t_start, smorf_end, strand, introns_transcript)

            ## Check the number of stop codons in the sequence
            s, s_index = find_stop_inframe(smorf_seq[:len(smorf_seq)-3], map_transc2gen) ## removes last codon and searches for stop codons inframe


            ## 1- check introns
            if strand == '+':
                ## check if start is within an intron
                start_intron = introns_transcript[(introns_transcript['start']<= smorf_start) & (introns_transcript['end']>= smorf_start)]
                ## check if end is within an intron
                end_intron = introns_transcript[(introns_transcript['start']<= smorf_end) & (introns_transcript['end']>= smorf_end)]

                if not start_intron.empty:
                    new_row = pd.DataFrame({
                        'smorf_id':[smorf_id],
                        'transcript_id': [t_id], 
                        'type': ['start_within_intron'], 
                        'length': ['-']})
                    unmatching_trancripts = pd.concat([unmatching_trancripts, new_row], ignore_index=True)
                    continue

                elif not end_intron.empty:
                    new_row = pd.DataFrame({
                        'smorf_id':[smorf_id],
                        'transcript_id': [t_id], 
                        'type': ['end_within_intron'], 
                        'length': ['-']})
                    unmatching_trancripts = pd.concat([unmatching_trancripts, new_row], ignore_index=True)
                    continue

            elif strand == '-': 
                ## check if start is within an intron
                start_intron = introns_transcript[(introns_transcript['start']<= smorf_end) & (introns_transcript['end']>= smorf_end)]
                ## check if end is within an intron
                end_intron = introns_transcript[(introns_transcript['start']<= smorf_start) & (introns_transcript['end']>= smorf_start)]
                
                if not start_intron.empty:
                    new_row = pd.DataFrame({
                        'smorf_id':[smorf_id],
                        'transcript_id': [t_id], 
                        'type':['start_within_intron'], 
                        'length': ['-']})
                    unmatching_trancripts = pd.concat([unmatching_trancripts, new_row], ignore_index=True)
                    continue
                
                elif not end_intron.empty:
                    new_row = pd.DataFrame({
                        'smorf_id':[smorf_id],
                        'transcript_id': [t_id], 
                        'type': ['end_within_intron'], 
                        'length': ['-']})
                    unmatching_trancripts = pd.concat([unmatching_trancripts, new_row], ignore_index=True)
                    continue
        

        ## 2- check 3nt periodicity
        elif smorf_len % 3 != 0: 
            new_row = pd.DataFrame({
                'smorf_id':[smorf_id],
                'transcript_id':[t_id], 
                'type': ['not_multiple_of_3'], 
                'length': [smorf_len]})
            unmatching_trancripts = pd.concat([unmatching_trancripts, new_row], ignore_index=True)
            continue


        ## 3- Last trio is not a stop
        elif smorf_seq[len(smorf_seq)-3:len(smorf_seq)+1] not in stop_codons:
            last_codon = smorf_seq[len(smorf_seq)-3:len(smorf_seq)+1]
            new_row = pd.DataFrame({
                'smorf_id':[smorf_id],
                'transcript_id':[t_id], 
                'type': ['last_trio_not_a_stop'], 
                'length': [last_codon]})
            unmatching_trancripts = pd.concat([unmatching_trancripts, new_row], ignore_index=True)
            continue

        
        ## 4- check multiple stop codons
        elif s != None: ## Multiple stop codons in the sequence 
            new_row = pd.DataFrame({
                'smorf_id':[smorf_id],
                'transcript_id': [t_id], 
                'type': ['more_than_one_stop'], 
                'length': ['-'] })
            unmatching_trancripts = pd.concat([unmatching_trancripts, new_row], ignore_index=True)
            continue
        
        ## transcript matches 
        matching_transcripts.append(t_id)
        transcripts_mapping_dictionary[t_id] = [map_gen2transc, map_transc2gen]


    return matching_transcripts, unmatching_trancripts, transcripts_mapping_dictionary



def check_var_type(ref, alt):
    """
        Function to check if the variant is an indel.

        Considers checks also the VCF format (anchor/no_anchor nt)
        Note: nomenclature is different, also variant position changes in one case or other

        Insertions do not have the no-anchor format as the variant position corresponds to a position in the reference genome.

        Input: 
        - ref: reference allele
        - alt: alternative allele

        Output: 
        var_type: snv (single nucleotide variant), indel (insertion-deletion), del (deletion), ins (insetion)
        vcf_format: anchor or no_anchor 
 
    """

    prefix_del = ref.startswith(alt)

    if alt in ['','.', '*'] and prefix_del == False: ## del without anchor nt
        vcf_format = 'no_anchor'
    else: 
        vcf_format = 'anchor'

    return vcf_format



def add_anchor_nt(var_pos, ref, alt, ref_genome):
    """
        Function to add the anchor nt to no-anchor format. 
        
        Note: smorfep functions were designed to process anchor-based variants. 
    """
    
    new_var_pos = var_pos - 1
    anchor_nt = get_sequence(new_var_pos, new_var_pos, '+',ref_genome)

    if alt == '*': ## deletion 
        new_ref = anchor_nt + ref
        new_alt = anchor_nt
        
    return new_var_pos, new_ref, new_alt 



def find_position(dct, position):
    """
        Function to quickly check is a position is in the dicionary of mapped positions.
        Iterated on the keys.

        NOTE: for the smorfep if the position is in the dictionary keys list -- means the position is in an exonic region.

        Single position search.

        Input:
        - dct: dictionary of position to search 
        - position: position of interest

        Returns a boolean result
    """

    if position in dct.keys():
        return True
    else: 
        return False



def within_exon(start, end, mapgen2transc):
    """
        Function to count how many nucleotides are in the exon region. 

        Used to determine if the indel is inframe or frameshift

        Input: 
        - start position
        - end position
        - dictionary with the exon genomic coordinates as keys

        Returns the number of nucleotides that fall wihin exonic region(s).
    """

    # exonnts = 0 ## count the number of 
    # for position in range(start, end):
    #     if position in mapgen2transc.keys():
    #         print('exon pos:', position)
    #         exonnts += 1

    positions_exon = [x for x in range(start, end+1) if x in mapgen2transc.keys()]

    return len(positions_exon)


## reported on the forward strand
def map_splice_regions(introns_df, splice_size, intron_exon_size=3, splice_da_size=2):

    """
        Function to create a dataframe with the splice regions considering the introns within a smORF. 
        Maps the region as follow: 
        (donor)
        - splice_region_start = intron_start - into_exon_size
        - splice_region_end = intron_start + splice-size
        (acceptor)
        - splice_region_start = intron_end - splice_size
        - splice_region_end = intron_start + into_exon_size

        Input: 
        - intron_coordinates: dataframe with the coordinates of intron regions
        - splice_size: user defined or default (8bps as VEP) splice region size 
        - intron_exon_size: number of nucleotides in the exon precede the intron (by default 2bps as VEP)
        - splice_da_size: size of the donor/acceptor region (by default 2bps as VEP)

        Returns a dataframe with the start and end coordinate of splice regions -- 2 per intron. 

        NOTE 1: donor and acceptor reported on the forward strand. For annotations on the reverse strand donor and acceptor should be reversed. 
                
        NOTE 2: by default the 3 last/first nucleotides in the exon (resp. donor/acceptor for forward strand, otherwise for reverse strand) are considered part of the splice_region (same as VEP).
        To change this, change into_exon_size variable in this function.

    """ 

    ##into_exon_size = 3

    ## final dataframe - start empty
    splice_regions_df = pd.DataFrame(data=None, columns=['chr', 'start', 'end', 'splice_region', 'ID', 'da_start', 'da_end'])
    
    new_line_index = 0
    for index, row in introns_df.iterrows():
        splice_region_donor_start = row.start - intron_exon_size
        splice_region_donor_end = row.start + splice_size - 1 ## -1 as row.start is the first position of the intron
        splice_region_acceptor_start = row.end - splice_size + 1 ## -1 as row.end is the last position of the intron
        splice_region_acceptor_end = row.end + intron_exon_size 

        num = row.intron_number
        t_id = row.transcript_id
        c = row.chr

        new_line_donor = pd.DataFrame(
            {
            'chr': c, 
            'start': splice_region_donor_start, 
            'end': splice_region_donor_end, 
            'splice_region': 'donor_sr_intron_'+str(num), 
            'ID': t_id,
            'da_start': row.start, 
            'da_end': row.start + splice_da_size -1
            }, index=[new_line_index]
        )

        new_line_index += 1
        splice_regions_df = pd.concat([splice_regions_df, new_line_donor])


        new_line_acceptor = pd.DataFrame(
            {
            'chr': c, 
            'start': splice_region_acceptor_start, 
            'end': splice_region_acceptor_end, 
            'splice_region': 'acceptor_sr_intron_'+str(num), 
            'ID': t_id,
            'da_start': row.end - splice_da_size +1, 
            'da_end': row.end
            }, index=[new_line_index]
        )

        new_line_index += 1        
        splice_regions_df = pd.concat([splice_regions_df, new_line_acceptor])
    
    return splice_regions_df


def check_introns(seq, start_orf, end_orf, var_pos, ref, alt, strand, map_gen2transc, splice_regions_df, splice_size =8, intron_exon_size=3, splice_da_size=2):
    """ 
        Function to check if a variant crosses exon-intron boundaries.
        # Special case, only required for indels.

        Assumes ref and alt in with-anchor-nt format.

        Input: 
        - var_pos: variant position
        - ref: reference allele
        - alt: alternative allele
        - strand: strand is required to set in which direction we need to search
        - map_gen2transc: mapping between genomic and transcript coordinates (used to obtain the exon and intron coordinates)
        - splice_regions_df: dataframe with splice regions 
        
        Note1: Intron coordinates are not included in the mapping, as introns are not present in the transcript sequence.

        Note2: donor is the interval on the left of the intron; acceptor is the end on the right end of the intron (strand independent).

        Note3: Protein consequence is not assessed for this kind of variants. the only case there is a protein consequence annotation
        is if the variant falls within the exonic, splice region nucleotides.

        Output: 
        Consequnce(s) of the variant for this case. 
        Or 'None' if the full length of the variant is within the exon --> run the normal analysis after.
    """

    prot_cons = '-' ## Protein consequence is not assessed for this kind of variants 


    var_type = ''
    ## NOTE: var_next_pos -- used to check if the indel is in the in between the last nt of the exon and the first of the intron -- Special case

    ## Get the postions affected by the variant
    ## range is ordered - crescent way
    if len(ref) == len(alt): ## SNV
        var_type = 'SNV'
        all_var_pos = None
    elif len(ref) < len(alt): ## insertion
        var_type = 'insertion'
        if strand == '+':
            all_var_pos = [v for v in range(var_pos, var_pos + len(alt))] ## no need for -1, as range is non-inclusive of end pos
            var_next_pos = var_pos+1
        elif strand == '-':
            all_var_pos = [v for v in range(var_pos-(len(alt)-1), var_pos+1)]
            var_next_pos = var_pos-1
    elif len(ref) > len(alt): ## deletion
        var_type = 'deletion'
        if strand == '+':
            all_var_pos = [v for v in range(var_pos, var_pos + len(ref))] ## no need for -1, as range is non-inclusive of end pos
            var_next_pos = var_pos+1
        elif strand == '-':
            all_var_pos = [v for v in range(var_pos - (len(ref)-1), var_pos+1)]
            var_next_pos = var_pos-1

    ##print(var_type)
    
    # print('all_var_positions:')
    # print(all_var_pos)


    ## Filter donor and acceptor regions -- should contain the variant position
    ## NOTE: donor and acceptor are mutually exclusive
    filtered_rows = splice_regions_df[splice_regions_df['splice_region'].str.contains('donor_sr_intron')] ## only donor
    ##filtered_donor = filtered_donor[(filtered_donor['start'] <= var_pos) & (filtered_donor['end'] >= var_pos)]
    ## gets the donor site of interest

    filtered_donor = pd.DataFrame(columns=splice_regions_df.columns)
    if var_type != 'SNV': ## to include exon-intron crossing variants
        for index, row in filtered_rows.iterrows():
            for pos in all_var_pos:
                if row['start'] <= pos <= row['end']:
                    filtered_donor = pd.DataFrame([row])
                    break

    else:
        filtered_donor = filtered_rows[(filtered_rows['start'] <= var_pos) & (filtered_rows['end'] >= var_pos)]


    filtered_rows_ac = splice_regions_df[splice_regions_df['splice_region'].str.contains('acceptor_sr_intron')] ## only acceptor
    ##filtered_acceptor = filtered_acceptor[(filtered_acceptor['start'] <= var_pos) & (filtered_acceptor['end'] >= var_pos)]  
    filtered_acceptor = pd.DataFrame(columns=splice_regions_df.columns)
    if var_type != 'SNV':
        for index, row_a in filtered_rows_ac.iterrows():
            for pos in all_var_pos:
                if row_a['start'] <= pos <= row_a['end']:
                    filtered_acceptor = pd.DataFrame([row_a])
                    break
    else: 
        filtered_acceptor = filtered_rows_ac[(filtered_rows_ac['start'] <= var_pos) & (filtered_rows_ac['end'] >= var_pos)]  

    ## save donor and acceptor positions
    ## for all the donor and acceptor sites per transcript 
    ## for donor/acceptor variants
    donor_acceptor_positions = []

    ## splice_regions (1-3 bases in the exon flanking the intron start + 7-XX base within the intron (XX for VEP is 8, but user can define other lenght))
    ## for splice_region_variant
    ## NOTE: left -- donor; right -- acceptor (foraward strand, reverse otherwise) 
    splice_region = []
    splice_region_exon_nts = [] ## used for the SNV case, as if a SNV falls into the exon bases of splice region there there is also missense/synonymous annotation

    ## Positions between 3-6 bases (VEP default) within the intron
    ## for splice_donor/acceptor_region_variant
    splice_donor_acceptor_region = []

    ## Positions between 3-17 bases from the acceptor end of the intron
    ## Acceptor ONLY
    polypirimidine_region = []

    ## 5th position within the intron donor side 
    ## Donor ONLY
    fifthbase = None

    intron_end_region = None

    ##XXtestXX
    ##print('donor empty, acceptor empty', filtered_donor.empty, filtered_acceptor.empty)
    if not filtered_donor.empty: 
        row = filtered_donor.iloc[0] ## left side
        if strand == '+':
            ##print('+ donor end')
            intron_end_region = 'donor_end'

            donor_acceptor_positions.extend([i for i in range(row.da_start, row.da_end+1)])

            splice_region.extend([m for m in range(row.start, row.start+intron_exon_size)]) 
            splice_region.extend([m for m in range(row.end-(splice_size-6)+1, row.end+1)]) ## VEP uses 6th base up to splice region upper range (8bps) as for splice region

            splice_region_exon_nts.extend([m for m in range(row.start, row.start+intron_exon_size)]) 
            #print(splice_region_exon_nts)

            fifthbase = row.da_start + 4 ## +4 as da_start is the first base of the intron

            splice_donor_acceptor_region.extend([t for t in range(row.da_start+splice_da_size, fifthbase)])
            splice_donor_acceptor_region.extend([fifthbase+1]) ## adds 6th base -- default VEP 

        elif strand == '-':
            ##print('- acceptor end')
            intron_end_region = 'acceptor_end'

            donor_acceptor_positions.extend([i for i in range(row.da_start, row.da_end+1)])
            ##print(donor_acceptor_positions)

            polypirimidine_region.extend([j for j in range(row.da_end+1,row.da_end+16)])
            ##print(polypirimidine_region)

            splice_region.extend([m for m in range(row.start, row.start+intron_exon_size)]) 
            splice_region.extend([m for m in range(row.end-(splice_size-6)+1, row.end+1)]) ## VEP uses 6th base up to splice region upper range (8bps) as for splice region

            splice_region_exon_nts.extend([m for m in range(row.start, row.start+intron_exon_size)]) 

            splice_donor_acceptor_region.extend([t for t in range(row.da_start+splice_da_size, row.end-(splice_size-6)+1)])

    elif not filtered_acceptor.empty:
        row_a = filtered_acceptor.iloc[0] ## right side

        if strand == '+':
            ##print('+ acceptor end')

            intron_end_region = 'acceptor_end'

            donor_acceptor_positions.extend([g for g in range(row_a.da_start, row_a.da_end+1)])
            ##print(donor_acceptor_positions)

            polypirimidine_region.extend([j for j in range(row_a.da_start-17, row_a.da_start)])
            ##print(polypirimidine_region)

            splice_region.extend([j for j in range(row_a.start, row_a.start+(splice_size-6))]) ## VEP uses 6th base up to splice region upper range (8bps) as for splice region
            splice_region.extend([j for j in range(row_a.end-intron_exon_size+1, row_a.end+1)]) 
            ##print(splice_region) 

            splice_region_exon_nts.extend([j for j in range(row_a.end-intron_exon_size+1, row_a.end+1)])

            splice_donor_acceptor_region.extend([t for t in range(row_a.da_start-(6-splice_da_size), row_a.da_start)]) ## as range is exclusive of end -- reports up to row_a.da_start-1 position
            ## OK
        
        elif strand == '-':
            ##print('- donor end')
            intron_end_region = 'donor_end'

            donor_acceptor_positions.extend([g for g in range(row_a.da_start, row_a.da_end+1)])

            splice_region.extend([j for j in range(row_a.start, row_a.start+(splice_size-6))]) ## VEP uses 6th base up to splice region upper range (8bps) as for splice region
            splice_region.extend([j for j in range(row_a.end-intron_exon_size+1, row_a.end+1)]) 

            splice_region_exon_nts.extend([j for j in range(row_a.end-intron_exon_size+1, row_a.end+1)])

            fifthbase = row_a.da_end - 4 ## -4 as da_end is the first base of the intron

            splice_donor_acceptor_region.extend([t for t in range(fifthbase+1, row_a.da_end-splice_da_size+1)])
            splice_donor_acceptor_region.extend([fifthbase-1]) ## adds 6th base -- default VEP 

    ##XXtestXX
    ##print(intron_end_region)

    if var_type != 'SNV': ## computations required only for indels
        ## variant start and end coordinates
        ##print(all_var_pos)
        var_start = all_var_pos[0]
        var_end = all_var_pos[-1]
        ##print('var_start', var_start)
        ##print('var_end', var_end)

        ## check start and end within exon: 
        var_pos_check = find_position(map_gen2transc, var_pos)
        var_start_check = find_position(map_gen2transc, var_start)
        var_end_check = find_position(map_gen2transc, var_end)
        ##print('var_pos', var_pos)
        ##print('var_pos_check', var_pos_check)
        ##print('var_end_check', var_end_check)

        var_next_pos_check = find_position(map_gen2transc, var_next_pos)
        ##print('var_next_pos_check', var_next_pos_check)

        ## don't need swap for reverse strand -- var_start is always < than var_end
        exon_nts = within_exon(var_start, var_end, map_gen2transc)
        ##print('exon_nts', exon_nts)

        if strand == '+': 
            check_no_anchor = all_var_pos[1:] ## first position in the range is the anchor base
        elif strand == '-':
            check_no_anchor = all_var_pos[:len(all_var_pos)-1] ## last position in the range is the anchor base

        
        ##XXtestXX
        ##print('positions no anchor', check_no_anchor)

    ## 1- donor splice region variants
    if intron_end_region == 'donor_end':

        if var_type == 'SNV': ## single nucleotide variant
            ##print('SNV')
            if var_pos in donor_acceptor_positions: ## if it is a deletion and affects the splice site is donor 
                dna_cons = 'splice_donor_variant'
                prot_cons = '-'
            
            elif var_pos == fifthbase: 
                dna_cons = 'splice_donor_5th_base_variant&intron_variant'
                prot_cons = '-'

            elif var_pos in splice_donor_acceptor_region:
                dna_cons = 'splice_donor_region_variant&intron_variant'
                prot_cons = '-'
        
            elif var_pos in splice_region: ## three possibilities
                ## VEP reports: 
                # missense_variant&splice_region_variant - if in the exon and changes the aminoacid
                # splice_region_variant&synonymous_variant - if in the exon and do NOT change the aminoacid
                # splice_region_variant&intron_variant - if in the intron 
                
                if var_pos in splice_region_exon_nts:

                    new_sequence, ref_original, ref_inFile = add_variant_transcriptSeq(seq, start_orf, end_orf, ref, alt, var_pos, map_gen2transc)
                    prot_cons, change_prot = protein_consequence_transcript(seq, new_sequence, var_pos, map_gen2transc)

                    ## should only report those two annotations
                    if prot_cons == 'missense_variant':
                        dna_cons = 'missense_variant&splice_region_variant'
                    elif prot_cons == 'synonymous_variant':
                        dna_cons = 'splice_region_variant&synonymous_variant'
                    elif prot_cons == 'protein_truncation':
                        dna_cons = 'stop_gained&splice_region_variant'
                    
                else:
                    print('var_position in splice region and not in the exon nucleotides region')
                    dna_cons = 'splice_region_variant&intron_variant'

            
            elif var_pos not in map_gen2transc.keys():
                dna_cons = 'intron_variant'
            
            else: ## run the deep intron and exon analysis 
                dna_cons = 'Not_intronic'



        ## insertions
        elif var_type == 'insertion': 
            ##print('insertion')

            insertion_size = len(alt) -1 ## -1 to remove anchor base

            if var_start_check == True and var_pos not in splice_region_exon_nts and var_next_pos in splice_region_exon_nts and strand == '+': ## anchor position is the position before the splice region exonic nts
                ##XXtestXX
                ## print('first condition - anchor before splice region')
                dna_cons = 'Not_intronic'

            elif var_start_check == True and var_end_check == True and [x for x in check_no_anchor if x in splice_region_exon_nts] != []: ## variant within the exon, but on the splice region -- last 3 nt of the exon (VEP default)

                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion

            elif var_next_pos in splice_region_exon_nts and strand == '+':  ## forward strand working
                ##XXtestXX
                ##print('var_next_position within splice region exon nts')

                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion

            elif var_pos_check == True and var_next_pos_check == False and [x for x in check_no_anchor if x in splice_region_exon_nts] == [] and strand == '+': ## insertion after the last exon nt ## forward strand
                ##XXtestXX
                ##print('var_pos check True, var_next_positio_check False, check_no_anchor within splice region exon')

                if insertion_size % 3 == 0:

                    ## Check if there is a stop inframe for all the insertion length -- not limited to the first codon as a stop
                    seq_new = seq[:map_gen2transc[splice_region_exon_nts[-1]]+1] ## from begining until the intron
                    seq_insertion = seq_new + alt[1:] ## adds all the alternative without the anchor nt
                    find_stop = find_stop_inframe_sequence(seq_insertion)
                    
                    if find_stop == True: ## if the insertion introduces a stop codon
                        dna_cons = 'stop_gained&splice_region_variant'
                    else: 
                        dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion

            ## NOTE: this and next condition to get the forward and reverse strand cases
            elif exon_nts >= 1 and var_end_check == False and var_start_check == True: ## insertion after the last nt in the exon

                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion

            elif exon_nts >= 1 and var_end_check == True and var_start_check == False: ## insertion after the last nt in the exon

                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion

            elif var_pos_check == True and var_next_pos_check == False: 

                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion
            
            elif var_pos_check == False and var_next_pos == True: 

                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion


            elif [x for x in check_no_anchor if x in donor_acceptor_positions] != []: ## if it is a insertion and affects the splice site is donor 
                dna_cons = 'splice_donor_variant'
            
            elif var_pos in donor_acceptor_positions and var_next_pos not in donor_acceptor_positions: ## added for the forward strand -- Also used for the reverse -- removed strand == '+'
                ##XXtestXX
                ##print('var_pos in donor acceptor positions, and var_nex_pos not in donor acceptor positions and strand +')
                dna_cons = 'splice_region_variant&intron_variant'


            ## insertion inframe with fifthbase special case 
            elif fifthbase in check_no_anchor and insertion_size%3 != 0: 
                ## For insertions that cross the 5th base VEP annotates with splice_donor_region_variant&intron_variant -- We match 
                ##dna_cons = 'splice_donor_5th_base_variant&intron_variant'
                
                ##XXtestXX
                ##print('fifthbase insertion not multiple of 3')
                dna_cons = 'splice_donor_region_variant&intron_variant' 

            elif [x for x in all_var_pos if x in donor_acceptor_positions] != [] and [x for x in all_var_pos if x in splice_donor_acceptor_region] != []: ## if the insertion happens between the donor main site and the splice_donor_region '-- insertion on the 3rd base within intron
                ##XXtestXX
                ##print('all position ins donor acceptor and splice fonor acceptor region too and strand +')
                dna_cons = 'splice_region_variant&intron_variant'


            elif [x for x in check_no_anchor if x in splice_donor_acceptor_region] != []:
                ##XXtestXX
                ##print('second donor region condition')
                dna_cons = 'splice_donor_region_variant&intron_variant'

        
            elif [x for x in check_no_anchor if x in splice_region] != []:
                ##XXtestXX
                ##print('last splice region condition')
                dna_cons = 'splice_region_variant&intron_variant'

            
            elif var_start_check == True and var_end_check == True: ## ## variant fully in the exon run the deep intron and exon analysis 
                dna_cons = 'Not_intronic'

            else: 
                if [x for x in check_no_anchor if x in map_gen2transc.keys()] == []:
                    dna_cons = 'intron_variant'
                
                else:
                    dna_cons = 'Not_intronic'
                    

        ## deletions
        elif var_type == 'deletion':
            ##print('deletion')

            if [x for x in check_no_anchor if x in donor_acceptor_positions] != []: ## if it is a deletion and affects the splice site is donor 
                dna_cons = 'splice_donor_variant'

            elif var_start_check == True and var_end_check == True and [x for x in check_no_anchor if x in splice_region] != []: ## variant within the exon, but on the splice region -- last 3 nt of the exon (VEP default)
                deletion_size = len(ref) -1 ## -1 to remove anchor base

                if deletion_size % 3 == 0: 
                    dna_cons = 'inframe_deletion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_deletion

            ## NOTE: This and next condition to consider forward and reverse strand
            elif exon_nts >= 1 and var_end_check == False and var_start_check == True: ## deletion after the last nt in the exon
                deletion_size = len(ref) -1 ## -1 to remove anchor base

                if deletion_size % 3 == 0: 
                    dna_cons = 'inframe_deletion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_deletion

            elif exon_nts >= 1 and var_end_check == True and var_start_check == False: ## insertion after the last nt in the exon
                deletion_size = len(ref) -1 ## -1 to remove anchor base

                if deletion_size % 3 == 0: 
                    dna_cons = 'inframe_deletion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_deletion
            

            elif fifthbase in check_no_anchor: 
                dna_cons = 'splice_donor_5th_base_variant&intron_variant'

            elif [x for x in check_no_anchor if x in splice_donor_acceptor_region] != []:
                dna_cons = 'splice_donor_region_variant&intron_variant'
            
            elif [x for x in check_no_anchor if x in splice_region] != []:
                dna_cons = 'splice_region_variant&intron_variant'

            elif var_start_check == True and var_end_check == True: ## ## variant fully in the exon run the deep intron and exon analysis 
                dna_cons = 'Not_intronic'
            else: 
                if [x for x in check_no_anchor if x in map_gen2transc.keys()] == []:
                    dna_cons = 'intron_variant'
                
                else:
                    dna_cons = 'Not_intronic'

 

    ## 2- acceptor splice region variants
    elif intron_end_region == 'acceptor_end':
        if var_type == 'SNV': 
            ##print('SNV acceptor')
            if var_pos in donor_acceptor_positions: ## if it is a deletion and affects the splice site is donor 
                dna_cons = 'splice_acceptor_variant'


            elif var_pos in splice_region: ## three possibilities
                ## VEP reports: 
                # missense_variant&splice_region_variant - if in the exon and changes the aminoacid
                # splice_region_variant&synonymous_variant - if in the exon and do NOT change the aminoacid
                # splice_region_variant&intron_variant - if in the intron 
                
                if var_pos in splice_region_exon_nts:

                    new_sequence, ref_original, ref_inFile = add_variant_transcriptSeq(seq, start_orf, end_orf, ref, alt, var_pos, map_gen2transc)
                    prot_cons, change_prot = protein_consequence_transcript(seq, new_sequence, var_pos, map_gen2transc)

                    ## should only report those two annotations
                    if prot_cons == 'missense_variant':
                        dna_cons = 'missense_variant&splice_region_variant'
                    elif prot_cons == 'synonymous_variant':
                        dna_cons = 'splice_region_variant&synonymous_variant'
                    elif prot_cons == 'protein_truncation':
                        dna_cons = 'stop_gained&splice_region_variant'

                else:
                    ##print('splice_region, not exon nts')
                    dna_cons = 'splice_region_variant&intron_variant'

            
            elif var_pos not in map_gen2transc.keys():
                dna_cons = 'intron_variant'
            
            else: ## run the deep intron and exon analysis 
                dna_cons = 'Not_intronic'
            

        elif var_type == 'insertion': 
            insertion_size = len(alt) -1 ## -1 to remove anchor base

            ##print('acceptor insertion')

            ##print('check no anchor in splice_region', [x for x in check_no_anchor if x in splice_region])
            ##print('check no anchor in splice_region_exon_nts', [x for x in check_no_anchor if x in splice_region_exon_nts])
            ##print('splice_region exon nts', splice_region_exon_nts)

            if len([x for x in all_var_pos if x in donor_acceptor_positions]) == len(donor_acceptor_positions) and var_pos in donor_acceptor_positions: ## insertion in the middle of the acceptor main site
                ##XXtestXX
                ##print('insertion in the middle of the main splice acceptor site')
                dna_cons = 'splice_acceptor_variant'


            elif var_start_check == True and var_end_check == True and [x for x in check_no_anchor if x in splice_region_exon_nts] != [] and strand == '+': ## variant within the exon, but on the splice region -- last 3 nt of the exon (VEP default)
                ##XXtestXX
                ##print('first protein altering condition')

                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion

            
            elif var_start_check == True and var_end_check == True and [x for x in check_no_anchor if x in splice_region_exon_nts] != [] and strand == '-':
                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion


            elif [x for x in check_no_anchor if x in donor_acceptor_positions] != []: ## if it is an insertion and affects the splice site is acceptor

                if len([x for x in all_var_pos if x in donor_acceptor_positions]) != len(all_var_pos): ## insertion overlaps, but is before the acceptor splice site - there is at least one nt not within the acceptor main splice site
                    print('variant overlaps splice region ')
                    if insertion_size % 3 == 0: ## multiple of 3 conditions added on 15-08-2023 working for strand +
                        print('multiple of 3')
                        if len([x for x in all_var_pos if x in donor_acceptor_positions]) > 1: ## condition need for splice_polypyrimidine_tract_variant&intron_variant annotation 15-08-2023 working for strand +
                            dna_cons = 'splice_polypyrimidine_tract_variant&intron_variant'
                        else:
                            dna_cons = 'splice_region_variant&splice_polypyrimidine_tract_variant&intron_variant'
                    else:
                        print('not multiple of 3')
                        dna_cons = 'splice_polypyrimidine_tract_variant&intron_variant'

                else: 
                    print('splice acceptor variant insertion affect the splice site but occurs before')
                    dna_cons = 'splice_acceptor_variant'
            

            elif exon_nts >= 1 and var_end_check == False and var_start_check == True and strand == '-': ## insertion after the last nt in the exon
                ##print('end false, start true, - strand')

                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion


            elif exon_nts >= 1 and var_end_check == True and var_start_check == False and strand == '-': ## insertion after the last nt in the exon -- reverse strand
                ##print('end true, start False, - strand')
                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion


            elif var_pos_check == True and var_next_pos_check == False: 
                ##print('vas_pos true, var_next false')
                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion

            
            elif var_pos_check == False and var_next_pos_check == True and strand == '+': ## insertion btween last nt in the intron and first of the exon -- forward strand
                print('insertion on the last intron base and chack protein change, + strand')
                
                if insertion_size % 3 == 0:
                    dna_cons = 'inframe_insertion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_insertion



            elif [x for x in all_var_pos if x in donor_acceptor_positions] != [] and [x for x in all_var_pos if x in splice_donor_acceptor_region] != []: ## if the insertion happens between the donor main site and the splice_donor_region '-- insertion on the 3rd base within intron
                print('all_pos -- splice donor_acceptor and splice donor_acceptor_region')
                dna_cons = 'splice_region_variant&intron_variant'


            elif [x for x in check_no_anchor if x in splice_region] != [] and strand == '+': ## added 15-08-2023
                print('no anchor var pos -- splice region and strand +')
                ##dna_cons = 'splice_region_variant&intron_variant'
                dna_cons = 'splice_polypyrimidine_tract_variant&intron_variant'


            elif [x for x in check_no_anchor if x in splice_region] != [] and strand == '-': ## edited on the 23-08-2023
                print('no anchor var pos -- splice region and strand - ')
                ##dna_cons = 'splice_region_variant&intron_variant'
                dna_cons = 'splice_polypyrimidine_tract_variant&intron_variant'

            
            elif var_start_check == True and var_end_check == True: ## ## variant fully in the exon run the deep intron and exon analysis 
                dna_cons = 'Not_intronic'

            else: 
                if [x for x in check_no_anchor if x in map_gen2transc.keys()] == []:
                    dna_cons = 'intron_variant'
                
                else:
                    dna_cons = 'Not_intronic'


        elif var_type == 'deletion':
            ##print('acceptor deletion')

            if [x for x in check_no_anchor if x in donor_acceptor_positions] != []: ## if it is a deletion and affects the splice site is donor 
                dna_cons = 'splice_acceptor_variant'

            elif var_start_check == True and var_end_check == True and [x for x in check_no_anchor if x in splice_region] != []: ## variant within the exon, but on the splice region -- last 3 nt of the exon (VEP default)
                deletion_size = len(ref) -1 ## -1 to remove anchor base

                if deletion_size % 3 == 0: 
                    dna_cons = 'inframe_deletion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_deletion


            ## NOTE: This and next condition to consider forward and reverse strand
            elif exon_nts >= 1 and var_end_check == False and var_start_check == True: ## deletion after the last nt in the exon
                deletion_size = len(ref) -1 ## -1 to remove anchor base

                if deletion_size % 3 == 0: 
                    dna_cons = 'inframe_deletion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_deletion


            elif exon_nts >= 1 and var_end_check == True and var_start_check == False: ## insertion after the last nt in the exon
                deletion_size = len(ref) -1 ## -1 to remove anchor base

                if deletion_size % 3 == 0: 
                    dna_cons = 'inframe_deletion&splice_region_variant'
                else: 
                    dna_cons = 'frameshift_variant&splice_region_variant' ## frameshift_deletion

            
            
            elif [x for x in check_no_anchor if x in splice_region] != []:
                print('no anchor check splice region')
                dna_cons = 'splice_region_variant&intron_variant'


            elif var_start_check == True and var_end_check == True: ## ## variant fully in the exon run the deep intron and exon analysis 
                dna_cons = 'Not_intronic'


            else: 
                if [x for x in check_no_anchor if x in map_gen2transc.keys()] == []:
                    dna_cons = 'intron_variant'

                else:
                    dna_cons = 'Not_intronic'
    
    else: ## variant not in the donor or acceptor splice region
        ## check intron 
        if var_type == 'SNV': 
            if var_pos not in map_gen2transc.keys():
                dna_cons = 'intron_variant'
            else:
                dna_cons = 'Not_intronic'

        else: 
            if [x for x in all_var_pos if x in map_gen2transc] == []:
                dna_cons = 'intron_variant'
            elif var_pos not in map_gen2transc.keys():
                dna_cons = 'intron_variant'
            elif var_next_pos_check == False: ## position after anchor is within the intron
                dna_cons = 'intron_variant'
            
            else:
                dna_cons = 'Not_intronic'
    
    print(dna_cons, '-', prot_cons, '-', all_var_pos)

    return dna_cons, '-', prot_cons, '-', all_var_pos




# def compute_start_end_coordinate(transcript_annotations): 
#     """
#     Function to compute the start and end coordinates of a combination of intervals.
#     Used for computing the start and end coordinate of genomic regions (e.g., 5'UTR, CDS, 3'UTR)
#     from .gff format files (as GENCODE).

#     Input: 
#     - transcript_annotations- pandas dataframe with the annotations for a given transcript

#     Output: 
#     (1) When annotations are available, the output has 3 ranges: 5'UTR, CDS, and 3'UTR.
#     (2) When annotations are NOT available, the output returns 'NA' in the undefined types.
#     """

#     ## select the section we want in the dataframe
#     ## five_prime_UTR
#     ## CDS
#     ## three_prime_UTR

#     ## if one is empty - Nothing
#     ## else
#     ## select stat minimum of start column and max from end column
#     ## output intervals






#     #return fiveprime_coord, cds_coord, threeprime_cood



def check_frame(smorf_start_coord, cds_start_coord):
    """
    Computes the difference between the smORF start and CDS sctart to determine if the 
    smORF start is in frame with the CDS start. 
    Input:
    - smorf_start_coord
    - CDS_start_coord

    Output: 
    - "inframe" - if the smorf and CDS start are in the same frame
    - "+1" - if smORF and CDS start are shifter of one position
    - "+2" - if smORF and CDS start are shifter of two positions
    """


    starts_difference = smorf_start_coord - cds_start_coord
    frame = starts_difference%3 ## -> how many nucleotides it shifts: 0 (inframe), 1 or 2

    if frame == 0: 
        return "inframe"
    elif frame in [-1, 1]:
        return "+1"
    elif frame in [-2, 2]:
        return "+2"
