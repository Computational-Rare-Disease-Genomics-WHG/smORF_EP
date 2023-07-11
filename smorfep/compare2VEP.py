#!/usr/bin/python3

## Contact: M Fernandes

## Script to comapre the annotations from smORF-EP and VEP

## VEP testing files were obtained by runing VEP on a VCF input file with the variants using the online ENSEMBL tool

## NOTE: working with the output from VEP web
## TODO: check how it does with the VEP output from the cluster

import argparse

from smorfep.utils.functions import *



def main():
    """
    Main entry point
    """

    parser = argparse.ArgumentParser(description='Script to compare smORF-EP output annotations with VEP output annotarions')

    parser.add_argument('-s', '--smorfepfilename', metavar='\b', required=True, type=str, help='smORF-EP annotations filename')
    parser.add_argument('-v','--vepfilename', metavar='\b', required=True, type=str, help='VEP annotations filename')
    parser.add_argument('-o', '--output', metavar='\b', default='comparison_output.vcf', type=str, help='outputname')

    args = parser.parse_args()

    ## read smorfep file
    smorfep_df = read_variants_file(args.smorfepfilename, '\t', 0)
    ##print(smorfep_df)

    ## read vep file
    ## command works when computed the VEP annotaitons on BMRC
    ##vep_df = read_vep_annotations(args.vepfilename)

    ## VEP annotations obtained VEP website API
    vep_df = read_vep_web(args.vepfilename)
    ##print(vep_df)
    
    
    comparison_df = pd.DataFrame(data=None, columns=['chrom', 'var_pos', 'ref', 'alt', 'var_id', 'transc_id', 'smorfep_so', 'vep_so', 'match'])


    index_gen = 0
    ## iterate per line in the smorfep file -- filters per transcript that cover the smORFs in analysis -- smaller possible df
    for index, row in smorfep_df.iterrows():

        consequence = row.DNA_consequence.split(', ')
        ##print(consequence)
        
        ## get the vep so for the same transcript 
        vep_annotation = vep_df[(vep_df.CHROM == row.chrm) \
                        & (vep_df.POS == row.var_pos) \
                        & (vep_df.REF == row.ref) \
                        & (vep_df.ALT == row.alt) \
                        & (vep_df.TRANSC_ID == row.transcript_id)]

        vep_so = vep_annotation.SO.item().split('&')
        ##print(vep_so)


        compare = [x for x in consequence if x in vep_so]
        ##print('compare', compare)

        if compare == consequence: 
            match = 'Same annotation'
        else: 
            match = 'Diff annotation'

        new_line = pd.DataFrame(
            {
            'chrom': row.chrm,
            'var_pos' : row.var_pos,
            'ref' : row.ref,
            'alt' : row.alt,
            'var_id' : row.var_id,
            'transc_id' : row.transcript_id,
            'smorfep_so' : row.DNA_consequence, 
            'vep_so' : vep_annotation.SO.item(), 
            'match' : match
            }, index=[index_gen]
        )

        comparison_df = pd.concat([comparison_df, new_line])
        
        index_gen += 1

    ## write the output to a file
    comparison_df.to_csv(args.output, sep='\t', lineterminator='\n', index=False)

    ## print the unmatching entries -- If empty all good
    missed_annotations = comparison_df[comparison_df.match == 'Diff annotation']
    print('Wrong annotations: ')
    print(missed_annotations) 



if __name__ == '__main__':
    main()