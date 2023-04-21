#!/usr/bin/python3

## Author: M Fernandes

## functions used by main script smORF-EP


from smorfep.utils.genetic_code import * 
from typing import List, Optional, Dict, Union
import pandas as pd
import re
import os

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
    fasta = seq

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
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    new_seq = ''.join([complement[base] for base in seq])
    
    return new_seq



def reverse_complement_seq(seq):
    """
        Function that applies both, reverse and complement of a DNA sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

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
        print(position, ref, alt, strand , start, end)
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
    map_tranccript2genome = {} ## dictionary transcCoord  --> genCoord

    if introns_df.empty: ## no introns case
        if strand == '+':
            pos_transc = 0
            for each in range(start, end): 
                map_genome2transcript[each] = pos_transc  
                map_tranccript2genome[pos_transc] = each

                pos_transc += 1

        elif strand == '-':
            pos_transc = 0
            for i in range(end, start, -1):
                map_genome2transcript[i] = pos_transc
                map_tranccript2genome[pos_transc] = i

                pos_transc += 1

    else: ## introns

        if strand == '+':
            introns_df = introns_df.sort_values(by=['start']) ## introns in crescent order

            intron_num = 1
            pos_transc = 0 ## index starts in 0
            for index, row in introns_df.iterrows(): 
                if intron_num == 1: ## first intron
                    for val in range(start, row['start']): ## the actual position is row['start']-1 which is the postion before the intron start
                        ## we set row['start'] as range is end exclusive, so it stops at row['start']-1
                        map_genome2transcript[val] = pos_transc
                        map_tranccript2genome[pos_transc] = val
                        
                        pos_transc += 1

                    next_start = row['end'] +1

                else: 
                    for v in range(next_start, row['start']):  ## Checked
                        map_genome2transcript[v] = pos_transc
                        map_tranccript2genome[pos_transc] = v

                        pos_transc += 1

                    next_start = row['end'] +1

                intron_num += 1
            
            ## last bit
            for val_end in range(next_start, end+1): ## +1 to include the end position 
                map_genome2transcript[val_end] = pos_transc
                map_tranccript2genome[pos_transc] = val_end
                pos_transc += 1

        elif strand == '-':
            introns_df = introns_df.sort_values(by=['start'], ascending=False) ## introns in crestcent order

            intron_num = 1
            pos_transc = 0 ## index starts in 0
            for index, row in introns_df.iterrows(): 
                if intron_num == 1: ## first intron

                    for val in range(end, row['end'],-1): ## genomic coord would be row['end']+1 -- python index starts at 0
                        map_genome2transcript[val] = pos_transc
                        map_tranccript2genome[pos_transc] = val
                        
                        pos_transc += 1
                        
                    next_start = row['start'] -1

                else: 
                    for v in range(next_start, row['end'], -1):  ## TODO: XXX to check 
                        map_genome2transcript[v] = pos_transc
                        map_tranccript2genome[pos_transc] = v
                        
                        pos_transc += 1

                    next_start = row['start'] -1
                
                intron_num += 1
            
            ## last bit
            for val_end in range(next_start, start-1, -1): 
                map_genome2transcript[val_end] = pos_transc
                map_tranccript2genome[pos_transc] = val_end
                pos_transc += 1


    return map_genome2transcript, map_tranccript2genome
    


## for sequences with introns
def add_variant_transcriptSeq(sequence, start, end,  ref, alt, position, map_coordinates):
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
        print(position, ref, alt, start, end)
        print('ref allele does not correspond!')
        print('reference genome: ', sequence[variant_index:variant_index+len(ref)])
        print('ref input: ', ref) 
        print(map_coordinates)
        print(position) 
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




def search_introns(introns, var_pos, strand, splice_site = 8): 
    """
        Function to search if the variants are within an intronic region. 
        Input: 
        - introns = pandas dataframe with the information of the introns formated as gff MANE files
        Contains: 
              .chr   .source  .type    .start   .end
              .col5  .strand  .col7    .ID  .Parent
              .gene_id   .transcript_id   .gene_type   .gene_name   .transcript_type 
              .transcript_name   .exon_number .exon_id .tag .protein_id  .Dbxref
        - strand 
        - var_pos = position of the variant
        - splice_site = range considered splice-site within the intron region, both ends
                    by defult this value is defines as 8bps as used by VEP.

        Assumes introns are sequential and do not overlap -- MANE case. Code needs changes otherwise.

        Return: 
        - result = None, intron var, splice-site var

        splice-site var if fall within 8bases at most from intron start or end 
    """

    result = ""


    ## check if the variant falls into an intron
    intron_select_start = introns[(introns.start <= var_pos) & (introns.end >= var_pos)]['start'].tolist()
    intron_select_start.sort() ## sorts all the starts for the intervals/introns the variant falls in 


    if len(intron_select_start) == 0: 
        result = 'Not intronic' 
    
    else: 
        ## if the variant is in the intron, but not within the splice-site
        result = 'intron_variant'

        ## check if affects the splice-sites
        ## assumes that it max gets 1 intron as correspondence - NO overlap

        if strand == '+':
            s = intron_select_start[0]
            e = introns[introns.start == s]['end'].item()

            ## if it gets within splice-site size(8bps default) from start of intron: splice donor
            if var_pos >= s+2 and var_pos <= s+ (splice_site-1):  ## splice_site -1 as the first base is s
                result = 'splice_region_variant'
            elif var_pos >= s and var_pos <= s+ 2: ## variant within the 2 first bases of the intron
                result = 'splice_donor_variant'
            
            ## if it gets within splice-site size(8bps default) from end of intron: splice acceptor
            elif var_pos <= e -2 and var_pos >= e - (splice_site-1):
                result = 'splice_region_variant'
            elif var_pos <= e and var_pos >= e - 2:
                result = 'splice_acceptor_variant'
        

        elif strand == '-': 
            s = intron_select_start[0] ## end on the intron in the reverse strand
            e = introns[introns.start == s]['end'].item() ## start of the intron in the reverse strand

            ## if it gets within splice-site size(8bps default) from start of intron: splice donor
            if var_pos >= s+2 and var_pos <= s+ (splice_site-1):  ## splice_site -1 as the first base is s
                result = 'splice_region_variant'
            elif var_pos >= s and var_pos <= s+ 2: ## variant within the 2 first bases of the intron
                result = 'splice_acceptor_variant'
            
            ## if it gets within splice-site size(8bps default) from end of intron: splice acceptor
            elif var_pos <= e -2 and var_pos >= e - (splice_site-1):
                result = 'splice_region_variant'
            elif var_pos <= e and var_pos >= e - 2:
                result = 'splice_donor_variant'

    return result


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

    else: ## insertions/delitions -- The protein sequence is fixed, so frameshift does not apply at this level
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
               start_pos - sequence start postion -- for strand '-' is the end coordinate
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

        return map_coordinates[new_stop_index+3], new_stop_index+3
    
    

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
            return 'start_lost', len_change, prot_cons, change_prot
        
        else:
            return None, len_change, prot_cons, change_prot

    elif strand == '-': 
        len_change = '-'
        prot_cons = '-'
        change_prot = '-'

        if variant_pos <= end and variant_pos >= end -2 and seq[:3] != new_sequence[:3]: ## last condition is for indels right after the start codon, as the coordinate will be the last nt of the start codon, which is kept the same

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
        return 'start_lost', len_change, prot_cons, change_prot
    
    else: 
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
            seq2transcEnd = get_sequence(end+1, transcript_info.end, transcript_info.strand, ref_sequence)

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
            seq2transcEnd = get_sequence(transcript_info.start, start, transcript_info.strand, ref_sequence)

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
        

def check_smorf_transcript(ref_sequence, transcript_info, introns_df, smorf_start, smorf_end, strand):
    """
    Function to check the compatibility of smORF and transcript coordinates. 

    Compares one smorf with all the transcripts it falls completely within.

    Input:
    - ref_sequence:  reference sequence for the chromosome in analysis
    - transcript_info: dataframe structure line with the information about the transcript (start, end, id, etc)
    - introns_df: intronic regions dataframe for the chromosome in analysis
    - start: start of the region of interest
    - end: end of the region of interest
    - strand 

    returns a list of matching transcripts, a dataframe with unmatching transcript and the respective reason for exclusion
    """

    matching_transcripts = []
    unmatching_trancripts = pd.DataFrame(columns=['transcript_id','flag', 'type', 'length'])

    for index, row in transcript_info.iterrows():

        ## 1- check smorf start/end within transcript


        ## 2- check introns
        if strand == '+':
            ## check if start is within an intron
            start_intron = introns_df[(introns_df['start']<= smorf_start) & (introns_df['end']>= smorf_start)]
            ## check if end is within an intron
            end_intron = introns_df[(introns_df['start']<= smorf_end) & (introns_df['end']>= smorf_end)]

            if not start_intron.empty:
                return 'wrong_sequence', 'start within intron'

            elif not end_intron.empty:
                return 'wrong_sequence', 'end within intron'

        elif strand == '-': 
            ## check if start is within an intron
            start_intron = introns_df[(introns_df['start']<= smorf_end) & (introns_df['end']>= smorf_end)]
            ## check if end is within an intron
            end_intron = introns_df[(introns_df['start']<= smorf_start) & (introns_df['end']>= smorf_start)]
            
            if not start_intron.empty:
                new_row = {'transcript_id': , 'flag': 'wrong_sequence', 'type':'start within intron' , 'length': '-' }
                unmatching_trancripts = unmatching_trancripts.append(new_row, ignore_index=True)
            
            elif not end_intron.empty:
                new_row = {'transcript_id': , 'flag': 'wrong_sequence', 'type': 'end within intron', 'length': '-' }
                unmatching_trancripts = unmatching_trancripts.append(new_row, ignore_index=True)

        ## 3- check periodicity


        ## 4- check multiple stop codons
        
        

