#!/usr/bin/python3

## Contact: M Fernandes

## script to generate the viewing files 

import argparse 
import datetime

from smorfep.utils.functions import *

def generate_igv_files(smorf_vars_filename, all_smorfs_coordinates_filename, outputname = None, clinvar_filename= None):
    """
        Function that runs the 
    """



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