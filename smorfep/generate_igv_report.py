#!/usr/bin/python3

## Contact: M Fernandes

## script to generate the viewing files 

import argparse 
import datetime

from smorfep.utils.functions import *

def generate_igv_files(smorf_vars_filename, all_smorfs_coordinates_filename, var_filename = 'variants.vcf',smorf_filename='main_smorf.bed', overlap_filename='overlap.bed', clinvar_filename= None):
    """
        Function that runs generates the three input files used by igv-reports from 2 mandatory and o
        
        Input:
        - smorf_vars_filename = output files rerported by smORF-EP 
        - all_smorfs_coordinates_filename = File with all the smORFs in the dataset analysed. 
            This file is used to check the smORFs overlapping the target smORF
        - var_filename = outputname of the variants VCF file. By default: 'variants.vcf'.
        - smorf_filename = outputname of the target smORF BED file. By default: 'main_smorf.bed'.
        - overlap_filename = outputname of the other smORFs that overlapt the target smORf region (+/- 1000 bases). 
            By default: 'overlap.bed'.
        - clinvar_filename = clinvar annotaitons file. In case the variant has annotated significance, it will be presented in the report.
            Two Clinvar information fields are collected: CLNSIG and CLNREVSTAT.

        Output: 
        - var_file (var_filename): File with the variants information to be shown in the igv report: 
            Chromosome, position, ref_allele, alt_allele, var_ID, predicted effect in the target smORF, target smORF, 
            clinical significance (CLNSIG), and clinical significance revision status (CNLREVSTAT). 
            These two last details are from ClinVar (if the variant is not found in ClinVar, these fields will contain '-';
            if the ClinVar file is not provided these fields will be filled with 'NA')
        - smorf_file (smorf_filename): File with the coordinates of the target smORF in BED format.
        - overlap_file (overlap_filename): File with the coordinates of other smORFs that overlap the target smORF region

    """

    


    return None



def run_igv_reports(smorf_vars_filename, all_smorfs_coordinates_filename, outputname = None, clinvar_filename= None):
    """
        Function that runs the igv-reports generation.

        Input: 
        - smorf_vars_filename = output files rerported by smORF-EP 
        - all_smorfs_coordinates_filename = File with all the smORFs in the dataset analysed. 
            This file is used to check the smORFs overlapping the target smORF
        - outputname = user defined name for the output file (if not defined it is set to 'report_<toady_date>.html')
        - clinvar_filename = clinvar annotaitons file. In case the variant has annotated significance, it will be presented in the report.
            Two Clinvar information fields are collected: CLNSIG and CLNREVSTAT.

        output: 
        html file with interactive report for a given smORF

    """

    ## 0 - Variables setup
    var_filename = 'variants.vcf'
    smorf_filename = 'main_smorf.bed'
    overlap_filename = 'overlap.bed'
    flanking_size = 1000
    if outputname == None: ## default name
        now = datetime.datetime.now().strftime('%Y-%m-%d')
        outputfilename = 'report_'+now+'.html'
    else: 
        outputfilename = outputname + '.html'
    

    ## 1 - generate the input files used by igv-ºreports library


    os.system('create_report '+ var_filename + ' --flanking '+str(flanking_size)+' --info-columns PRED_EFFECT SMORF CLNSIG CLNREVSTAT --tracks '+smorf_filename+' '+var_filename+' '+overlap_filename+' --output ' +outputfilename)

    print('report generated')


def main():
    """
    Main entry point
    """

    parser = argparse.ArgumentParser(description='Script to generate an interactive report for specific smORFs')


    parser.add_argument('-s','--smorfFile', metavar='\b', required=True, type=str, help='smorf-variants filename')
    parser.add_argument('-a','--allsmORFsFile', metavar='\b', required=True, type=int, help='all smorfs to compare filename')
    parser.add_argument('-o', '--outputfile', metavar='\b', type=str, help='output filename prefix. Default format: .html')

    args = parser.parse_args()




if __name__ == '__main__':
    main()