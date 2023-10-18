#!/usr/bin/python3

## Contact: M Fernandes

## script to generate the viewing files 

## NOTE 1: input filenames for igc-reports are manually set in the code: Lines 
## NOTE 2: flanking size used to compute the overlaping smORFs is manually defined in Line 


import argparse 
import datetime
import numpy

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

    ## 1 - Compile input files for igv-reports
    ## read variants to a pandas dataframe
    variants_df = read_variants_file(smorf_vars_filename, '\t', 0)

    print(variants_df)

    ## create VCF file df - CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    column_vcf = ['CHRM','POS', 'ID', 'REF','ALT','QUAL','FILTER','INFO']
    vcf_output_df = pd.DataFrame(data=None, columns=column_vcf)
    ##print(vcf_output_df)

    vcf_output_df['CHRM'] = 'chr'+variants_df['chrm'].astype(str)
    vcf_output_df['POS'] = variants_df['var_pos']
    vcf_output_df['ID'] = variants_df['var_id']
    vcf_output_df['REF'] = variants_df['ref']
    vcf_output_df['ALT'] = variants_df['alt']
    vcf_output_df['QUAL'] = '.'
    vcf_output_df['FILTER'] = '.'
    vcf_output_df['INFO'] = 'PRED_EFFECT='+variants_df['DNA_consequence']+';SMORF='+variants_df['smorf_id']+';'

    print(vcf_output_df)

    ## TODO: Add ClinVar annotations
    ## Add clinvar annotations to INFO column
    ##CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter



    ## 2 - compile the BED file with the target smORF info
    ## NOTE: considers the start and end positions in the vairants-smORF file
    ## BED file with no header
    ##     18	47808957	47930659	SMAD2	2	-
    ## TODO: Check to include introns in the smORF track

    ##check first line in the variants-smorf df
    first_line = variants_df.iloc[0,:]
    ##print(first_line)

    smorf_chrom = first_line['chrm']
    smorf_start = int(first_line['start'])
    smorf_end = int(first_line['end'])
    smorf_id = first_line['smorf_id']
    smorf_score = '-'
    smorf_strand = first_line['strand']

    smorf_bed = open(smorf_filename,'w')
    smorf_bed.write(str(smorf_chrom)+'\t'+str(smorf_start)+'\t'+str(smorf_end)+'\t'+smorf_id+'\t'+str(smorf_score)+'\t'+smorf_strand)
    smorf_bed.close()
    ## STEP -2 OK


    ## 3 - compile the overlapping smorfs file
    ## NOTE: To add a flanking region edit next line
    extended_overlap_size = 0 ##1000 ## NOTE: This means 500 bases each end are considered for the extended overlap

    extended_min = smorf_start - extended_overlap_size/2 
    extended_max = smorf_end - extended_overlap_size/2

    ## open other smORFs coordinates files 
    smorfs_set_df = read_file(all_smorfs_coordinates_filename, '\t', None)
    smorfs_set_df.columns = ['chr','start','end','id','score','strand']
    ##print(smorfs_set_df)
    
    ## filtering
    ## TODO: allow with and without prefix crm column
    smorfs_set_df = smorfs_set_df[smorfs_set_df['chr'] == 'chr'+str(smorf_chrom)] ## filter per chromosome
    overlap_left_df = smorfs_set_df[(smorfs_set_df['start'] <= extended_min) & (smorfs_set_df['end'] >= extended_min) & (smorfs_set_df['end'] <= extended_max)] ## overlap a region on the left end
    overlap_right_df = smorfs_set_df[(smorfs_set_df['start'] >= extended_min) & (smorfs_set_df['start'] <= extended_max) & (smorfs_set_df['end'] <= extended_max)] ## overlap a region on the right end
    overlap_full_within_df = smorfs_set_df[(smorfs_set_df['start'] >= extended_min) & (smorfs_set_df['end'] <= extended_max)] ## overlapping smORFs fully within the range of main smorf
    overlap_full_over_df  = smorfs_set_df[(smorfs_set_df['start'] <= extended_min) & (smorfs_set_df['end'] >= extended_max)] ## overlapping smORFs are larger than the main smorf
    
    overlap_final_df = pd.concat([overlap_left_df,overlap_right_df,overlap_full_within_df,overlap_full_over_df])
    
    print(overlap_final_df)


    return None



def run_igv_reports(smorf_vars_filename, all_smorfs_coordinates_filename, clinvar_filename= None, outputname = None):
    """
        Function that runs the igv-reports generation.

        Input: 
        - smorf_vars_filename = output files rerported by smORF-EP 
        - all_smorfs_coordinates_filename = File with all the smORFs in the dataset analysed. 
            This file is used to check the smORFs overlapping the target smORF
        - clinvar_filename = clinvar annotaitons file. In case the variant has annotated significance, it will be presented in the report.
            Two Clinvar information fields are collected: CLNSIG and CLNREVSTAT.
        - outputname = user defined name for the output file (if not defined it is set to 'report_<toady_date>.html')
        
        output: 
        html file with interactive report for a given smORF

    """

    ## 0 - Variables setup
    ## NOTE: Variables named manually -- to change edit the next 3 lines
    var_filename = 'variants.vcf'
    smorf_filename = 'main_smorf.bed'
    overlap_filename = 'overlap.bed'
    flanking_size = 1000 ## NOTE: Flanking size = 1000 means -500 from the first variant position and up to +500 after the last variant position
    if outputname == None: ## default name
        now = datetime.datetime.now().strftime('%Y-%m-%d')
        outputfilename = 'report_'+now+'.html'
    else: 
        outputfilename = outputname + '.html'
    

    ## 1 - generate the input files used by igv-ºreports library
    generate_igv_files(smorf_vars_filename, all_smorfs_coordinates_filename, var_filename, smorf_filename, overlap_filename)


    os.system('create_report '+ var_filename + ' --flanking '+str(flanking_size)+' --info-columns PRED_EFFECT SMORF CLNSIG CLNREVSTAT --tracks '+smorf_filename+' '+var_filename+' '+overlap_filename+' --output ' +outputfilename)

    print('report generated')


def main():
    """
    Main entry point
    """

    parser = argparse.ArgumentParser(description='Script to generate an interactive report for specific smORFs')

    parser.add_argument('-v','--varsFile', metavar='\b', required=True, type=str, help='variants filename')
    parser.add_argument('-a','--allsmORFsFile', metavar='\b', type=str, help='all smorfs to compare filename')

    parser.add_argument('-o', '--outputfile', metavar='\b', type=str, help='output filename prefix. Default format: .html')

    args = parser.parse_args()

    run_igv_reports(args.varsFile, args.allsmORFsFile)




if __name__ == '__main__':
    main()