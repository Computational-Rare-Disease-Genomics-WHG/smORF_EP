#!/usr/bin/python3

## Contact: M Fernandes

## main tool function 
## Performs the variants annotation 

from smorfep.utils.functions import *


def tool(ref_sequence, transcript_info, transcript_introns_df, start, end, strand, ref, alt, variant_pos, map_gen2transc, map_transc2gen, splice_site=8, intron_exon_size=3, donor_acceptor_size = 2):

    """
        Function that runs the variant consequence check.

        Input:
        - ref_sequence: reference sequence for the chromosome in analysis
        - transcript_info: dataframe structure line with the information about the transcript (start, end, id, etc)
        - transcript_introns_df: intronic regions dataframe for the chromosome in analysis
        - start: start of the region of interest
        - end: end of the region of interest
        - strand 
        - ref: reference allele
        - alt: alternative allele
        - variant_pos: variant genomic position
        - splice_site: size from the start/end of exonic region to be considered splice site. 
                        ## 1-2 bases are assumed as donor and acceptor
                        By default, this is set to 8 (VEP used range)

        returns the variant consequence
    """

    stop_codons = ['TAG', 'TAA', 'TGA']

    ## transcript information - from input 1 line dataframe transcript_info
    t_id = transcript_info.iloc[0].transcript_id
    t_start = transcript_info.iloc[0].start
    t_end = transcript_info.iloc[0].end
    t_strand = transcript_info.iloc[0].strand 

    print('transcript start and end:')
    print(t_start, t_end)

    ## testing print -- variant header
    # print('\n')
    # print(variant_pos, ref, alt)

    ## 1 - Get sequence from Ref genome
    ##print(strand)
    ##print(start, end)
    seq = get_sequence(start, end, strand, ref_sequence)
    print('')
    print('sequence')
    print(seq)
    # print('seq len before variant ', len(seq))
    # print('diff end - start ', end-start)


    ## check if the format is deletion without anchor: 
    check_anchor_nt = check_var_type(ref, alt)
    if check_anchor_nt == 'no_anchor':
        ##print('no anchor')
        ## convert to with-anchor-nt format
        variant_pos, ref, alt = add_anchor_nt(variant_pos, ref, alt, ref_sequence)
        ## Checked: varainat reported on the forward strand -- OK

    ## introns on the extension -- pre-compute
    if strand == '+':
        ## smorf start until transcript end
        transcript_introns_df = transcript_introns_df[(transcript_introns_df['start']>= start) & (transcript_introns_df['end']<= t_end)]
        transcript_introns_df_extension = transcript_introns_df[(transcript_introns_df['start']>= end+1) & (transcript_introns_df['end']<=t_end)]
    

    elif strand == '-':
        ## transcript start until smORF end
        transcript_introns_df =  transcript_introns_df[(transcript_introns_df['start']>= t_start) & (transcript_introns_df['end']<= end)]
        transcript_introns_df_extension = transcript_introns_df[(transcript_introns_df['start']>= t_start) & (transcript_introns_df['end']<=start-1)]
        

    ## introns of the smorf region
    introns_smorf = transcript_introns_df[(transcript_introns_df['start']>= start) & (transcript_introns_df['end']<=end)]
    ##sort introns
    introns_smorf = introns_smorf.sort_values(by=['start'])
    ##print(introns_smorf[['start','end']])

    # Collects all the coordinates from the smORF start until the transcript end
    # Collect also the extension, from end of region until end of transcript
    ## considers two cases, with and without introns
    if strand == '+':
        ##map_gen2transc, map_transc2gen = genome2transcript_coords(start, t_end, strand, transcript_introns_df)
        
        if not transcript_introns_df_extension.empty:
            extension_seq, ext_len = remove_introns(transcript_introns_df_extension, end+1, t_end, strand, ref_sequence)

        else: 
            extension_seq = get_sequence(end+1, t_end, strand, ref_sequence)
            ext_len = len(extension_seq)
    
    elif strand == '-':
        ##map_gen2transc, map_transc2gen = genome2transcript_coords(t_start, end, strand, transcript_introns_df)

        if not transcript_introns_df_extension.empty:
            extension_seq, ext_len = remove_introns(transcript_introns_df_extension, t_start, start, strand, ref_sequence)

        else: 
            extension_seq = get_sequence(t_start, start, strand, ref_sequence)
            ext_len = len(extension_seq)
        ## Extension correct - 2025-01-08
        ##print(extension_seq) 

        ## Variant conversion to reverse strand
        ## pre-process variant on the reverse strand
        ref = reverse_complement_seq(ref)
        alt = reverse_complement_seq(alt)
\
    ## 1.1 check reference allele matching
    reference_allele = get_sequence(variant_pos-1, variant_pos-1+len(ref), strand, ref_sequence) ## the same for forward and reverse strand: OK 2025-01-08
    ##TEST reference allele is the right varaint position: OK 2025-01-08
    ##print('test_reference:')
    ##print(reference_allele, ref)
    #print(get_sequence(variant_pos-2, variant_pos+len(ref), strand, ref_sequence)) ## check the neighboor nucleotides to confirm the positon is correct
    if ref != reference_allele:
        return 'Reference_mismatch', 'ref_genome: '+ str(reference_allele), 'reference_given: ' + str(ref), '-' 


    ## 2 - Processing presence of introns in the smORF
    if not introns_smorf.empty:  ## if there are introns in the smORF range
        print('smORF with introns')

        ## map the splice regions in the smorf
        splice_regions_df = map_splice_regions(introns_smorf, splice_site)
        ##print(splice_regions_df.head)

        ## 2.1- Get sequence without introns
        seq, new_len = remove_introns(introns_smorf, start, end, strand, ref_sequence)
        ## region sequence without introns
        print(introns_smorf)
        print(introns_smorf['start'])
        print(introns_smorf['end'])
        print(seq)


        ## 2.2 - Check if variant falls into an intron region     
        ## notes:
        ## - Use 8bp (VEP standard) for splice-site affecting var
        ## - for the analysis in GEL - Update AggV2; use all intron length to spot variants (also for denovo)  
        
        ## Check: exon-intron crossing variants
        try: 

            dna_c, dna_seq_c, prot_c, prot_seq_c, all_var_pos = check_introns(seq, start, end, variant_pos, ref, alt, strand, map_gen2transc, splice_regions_df)

        except UnboundLocalError:
            print('no dna_cons')
            print(variant_pos, ref, alt, strand)
        
        ##print(' ')
        ##print(variant_pos, ref, alt)
        ## print(dna_c, dna_seq_c, prot_c, prot_seq_c)

        ## 2.3 -- not intron related variants
        if dna_c == 'Not_intronic':  ## runs if the variant is not intronic

            ## Computing splice region coordinates
            if strand == '+':
                donor_splice_region_exon = []
                ## compute the donor regions in the sequence range
                for index,row in splice_regions_df.iterrows():
                    donor_splice_region_exon.extend([m for m in range(row.start, row.start+intron_exon_size)]) 

            elif strand == '-':
                donor_splice_region_exon = []
                for index,row in splice_regions_df.iterrows():
                    donor_splice_region_exon.extend([j for j in range(row.end-intron_exon_size+1, row.end+1)])


            ## 2.3.1- introduce the variant 
            new_sequence, ref_original, ref_inFile = add_variant_transcriptSeq(seq, start, end, ref, alt, variant_pos, map_gen2transc)
            if new_sequence == None: 
                print('add_variant_transcript')
                return 'Reference_mismatch', 'ref_genome: '+ str(ref_original), 'reference_given: ' + str(ref_inFile), '-'


            if len(ref) > len(alt): ## deletion
                new_sequence = new_sequence.replace('-','') # we need to take the dash out for seq processing 

            ## 2.3.2- check start and stop affecting variants
            
            ## 2.5.2.1 - start related
            ## working with positions allows non-canonical starts
            start_var, len_change, prot_cons, change_prot = check_start_transcript(seq, new_sequence, variant_pos, map_gen2transc)
            if start_var != None:        
                ref_seq_ref_allele = get_sequence(variant_pos, variant_pos+len(ref), strand, ref_sequence)
                return start_var, len_change, prot_cons, change_prot
            
            ## 2.5.2.2 - stop related
            stop_var, len_change, prot_cons, change_prot = check_stop_transcript(seq, new_sequence, start, end, variant_pos, strand, map_gen2transc, map_transc2gen, extension_seq)
            if stop_var != None:  
                ref_seq_ref_allele = get_sequence(variant_pos, variant_pos+len(ref), strand, ref_sequence)
                return stop_var, len_change, prot_cons, change_prot


            ## 2.3.3 - All other types of variants

            ## 2.3.3.1 - stop gain
            seq_trios = get_trios(new_sequence)
            seq_trios.pop()

            ## search if there are stop codons in the sequence, after removing the last codon = region stop
            matching = [s for s in seq_trios if any(xs in s for xs in stop_codons)] 


            transcript_var_positon = map_gen2transc[variant_pos]
            if len(alt) > len(ref): ## insertions 
                var_trio = int(transcript_var_positon//3.0 +1)
            else: 
                var_trio = int(transcript_var_positon//3.0)

            if matching != [] and seq_trios[var_trio] in stop_codons:
                first_stop_found = matching[0]
                index_stop = seq_trios.index(first_stop_found)

                new_seq = new_sequence[:3*(index_stop+1)]

                len_change = len(new_seq) - len(seq)

                prot_cons, change_prot = protein_consequence_transcript(seq, new_seq, variant_pos, map_gen2transc)

                return 'stop_gained', len_change, prot_cons, change_prot


            ## 2.3.3.2 - Insertions
            if len(alt) > len(ref):
                ##print('insertion')

                len_change = len(alt) - len(ref)

                ## inframe
                if len(new_sequence) % 3 == 0 and len_change % 3 == 0: 

                    prot_cons, prot_change = protein_consequence_transcript(seq, new_sequence, variant_pos, map_gen2transc)
                    
                    ## all_var_pos compiled all the positions the variant includes
                    if strand == '+':
                        no_anchor_all_var_pos = all_var_pos[1:]
                        ##print('excluded anchor all_var_pos', no_anchor_all_var_pos)
                    elif strand == '-':
                        no_anchor_all_var_pos = all_var_pos[:len(all_var_pos)]
                        ##print('excluded anchor all_var_pos', no_anchor_all_var_pos)

                    ##print('splice_region_exon_nts', donor_splice_region_exon) ## splice_region_exon_nts is empty if the var does not fall into it
                
                    ##print('inframe_insertion')
                    return 'inframe_insertion', len_change, prot_cons, prot_change


                ## frameshit insertion
                else: 
                    new_seq = frameshift(new_sequence, extension_seq, map_transc2gen)

                    if new_seq != None: ## stop found within the transcript
                        len_change = len(new_seq) - len(seq)
                        ## protein consequence
                        prot_cons, change_prot = protein_consequence(seq, new_seq,variant_pos, start, end, strand)
                    else: 
                        len_change = 'off_transcript_stop'
                        prot_cons = '-'
                        change_prot = '-'

                    return 'frameshift_variant', len_change, prot_cons, change_prot


            ## 2.3.3.3 - Deletions
            elif len(ref) > len(alt):

                ## inframe
                if len(new_sequence) % 3 == 0: 
                    len_change = len(new_sequence) - len(seq)

                    prot_cons, prot_change = protein_consequence_transcript(seq, new_sequence, variant_pos, map_gen2transc)
                    
                    return 'inframe_deletion', len_change, prot_cons, change_prot

                ## frameshift deletion
                else: 
                    new_seq = frameshift(new_sequence, extension_seq, map_transc2gen)

                    if new_seq != None: ## stop found within the transcript
                        len_change = len(new_seq) - len(seq)
                        prot_cons, change_prot = protein_consequence(seq, new_seq, variant_pos, start, end, strand)
                    else: 
                        len_change = 'off_transcript_stop'
                        prot_cons = '-'
                        change_prot = '-'
                    
                    return 'frameshift_variant', len_change, prot_cons, change_prot

            prot_cons, change_prot = protein_consequence_transcript(seq, new_sequence, variant_pos, map_gen2transc)

            if prot_cons == 'missense_variant':
                return 'missense_variant', 0, prot_cons, change_prot
            elif prot_cons == 'synonymous_variant':
                return 'synonymous_variant', 0, prot_cons, change_prot

        else:
            return dna_c, dna_seq_c, prot_c, prot_seq_c
        ## ------- end of sequence with introns check  -------------


    ## 3 -- processing for smORFs without introns
    else: ##  NOTE: if extension are needed we still need to check introns on the extension

        # ## 3.1 - Introduce the variant
        new_sequence, ref_original, ref_inFile = add_variant(seq, start, end, ref, alt, variant_pos, strand)
        if new_sequence == None: 
            print('add_variant')
            return 'Reference_mismatch', 'ref_genome: '+ str(ref_original), 'reference_given: ' + str(ref_inFile), '-'


        if len(ref) > len(alt): ## deletion
            new_sequence = new_sequence.replace('-','') # we need to take the dash out for seq processing 
            
        ## 3.2 Start and stop variants

        if transcript_introns_df_extension.empty: ## no introns on the extension

            ## 3.4.1 - affect start
            start_var, len_change, prot_cons, change_prot = check_start(seq, new_sequence, start, end, variant_pos, strand)
            if start_var != None:    
                return start_var, len_change, prot_cons, change_prot

            ## 3.2.2 - affect stop
            stop_var, len_change, prot_cons, change_prot = check_stop(seq, new_sequence, start, end, variant_pos, strand, transcript_info, ref_sequence, map_transc2gen)
            if stop_var != None: 
                return stop_var, len_change, prot_cons, change_prot
        else: 
            ## 3.2.3 - start related
            # print(variant_pos)
            # print(variant_pos in map_gen2transc.keys())
            # print(start in map_gen2transc.keys())
            # print(end in map_gen2transc.keys())
            # print(map_gen2transc.keys())

            start_var, len_change, prot_cons, change_prot = check_start_transcript(seq, new_sequence, variant_pos, map_gen2transc)
            if start_var != None:     
                return start_var, len_change, prot_cons, change_prot
            
            ## 3.2.4 - stop related
            stop_var, len_change, prot_cons, change_prot = check_stop_transcript(seq, new_sequence, start, end, variant_pos, strand, map_gen2transc, map_transc2gen, extension_seq)
            if stop_var != None:  
                return stop_var, len_change, prot_cons, change_prot

  
        ## 3.3 - all other types of variants
        
        ## 3.3.1 Stop gain - new stop codon before the one in the original seq
        seq_trios = get_trios(new_sequence)
        
        seq_trios.pop() ## remove the last trio -- conventional stop codon
        
        matching = [s for s in seq_trios if any(xs in s for xs in stop_codons)] 

        ## check in which codon the variant happens, so if does not create a stop here is not considered a stop gain
        if strand == '+':
            if len(alt) > 3 and len(ref) == 1: ## deletions
                var_trio = int((variant_pos - start)//3.0 +1)

            else:
                var_trio = int((variant_pos - start)//3.0)

        elif strand == '-':
            if len(alt) > len(ref): ## insertions
                var_trio = int((end - variant_pos)//3.0 +1)

            else:
                var_trio = int((end - variant_pos)//3.0) 

        if matching != [] and seq_trios[var_trio] in ['TAG', 'TAA', 'TGA']: ## if we find a stop codon it was created by the variant
        ## second part of condition -- the codon where the variant happens is a stop codon

            first_stop_found = matching[0]
            index_stop = seq_trios.index(first_stop_found)

            new_seq = new_sequence[:3*(index_stop+1)] ## +1 as the index count starts in 0

            len_change = len(new_seq) - len(seq)

            prot_cons, change_prot = protein_consequence(seq, new_seq, variant_pos, start, end, strand)
            

            return 'stop_gained', len_change, prot_cons, change_prot


        ## 3.4 - insetions
        if len(ref) < len(alt): 
            len_change = len(alt) - len(ref)

            ## 3.4.1 inframe
            if len(new_sequence) % 3 == 0 and len_change % 3 == 0: 

                prot_cons, prot_change = protein_consequence(seq, new_sequence, variant_pos, start, end, strand)

                return 'inframe_insertion', len_change, prot_cons, prot_change

            ## 3.4.2 frameshift insertion
            else:
                new_seq = frameshift(new_sequence, extension_seq, map_transc2gen)

                if new_seq != None: 
                    len_change = len(new_seq) - len(seq)
                    ## protein consequence
                    prot_cons, change_prot = protein_consequence(seq, new_seq,variant_pos, start, end, strand)
                else: 
                    len_change = 'off_transcript_stop'
                    prot_cons = '-'
                    change_prot = '-'

                return 'frameshift_variant', len_change, prot_cons, change_prot

        ## 3.5 - deletions
        elif len(ref) > len(alt): 

            ## 3.5.1 inframe
            if len(new_sequence) % 3 == 0: 

                len_change = len(new_sequence) - len(seq)

                prot_cons, change_prot = protein_consequence(seq, new_sequence, variant_pos, start, end, strand)

                return 'inframe_deletion', len_change, prot_cons, change_prot

            ## 3.5.2 frameshift deletion -- WORKING 2023-01-24
            else:
                new_seq = frameshift(new_sequence, extension_seq, map_transc2gen)

                if new_seq != None: ## stop found within the transcript
                    len_change = len(new_seq) - len(seq)

                    ##protein consequence
                    prot_cons, change_prot = protein_consequence(seq, new_seq, variant_pos, start, end, strand)
                    ##print(prot_cons, change_prot)

                else: ## no new stop within the transcript 
                    len_change = 'off_transcript_stop'
                    prot_cons = '-'
                    change_prot = '-'

                return 'frameshift_variant', len_change, prot_cons, change_prot


        ## 3.6 - single nucleotide change
        prot_cons, change_prot = protein_consequence(seq, new_sequence, variant_pos, start, end, strand)

        if prot_cons == 'missense_variant':
            return 'missense_variant', 0, prot_cons, change_prot
        elif prot_cons == 'synonymous_variant':
            return 'synonymous_variant', 0, prot_cons, change_prot
        
