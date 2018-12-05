#!/bin/python
import pandas as pd
import os
import sys
import argparse

def find_interest(means, op):
    '''
    Function to determine genes of interes
    Inputs:
	means: Vector of transcripts and the mean for that transcipt
               for a given clade
        op:    Unary Predicate to determine which transcirpts to return
    Outpus:
        index of genes where op(mean) is true
    '''
    return means[op(means)].index

def find_match(df, sample_list):
    '''
    Function to return the samples in df are also in sample_list  
    Inputs:
         df - dataframe to find matches in
         sample_list - list of clade to find matches for
    Outpus:
         a list of all the column names where a sample name provided by sample_list is in the column name
    '''
    return [col_name for col_name in df.columns  if col_name in sample_list]

def parse_meta_data(meta_data):
    '''
    Determine which samples belong to which clade.
    Inputs:
        meta_data: String, path to metadata direcotry, assumes structre of:
            meta_data
                |
                |-clade1.txt
                |
                |-clade2.txt
            where the files contains the sample names for that clade seperated
            by newlines('\n')
    Outputs:
         Dictionary that maps a clade to the samples listed within that clades file.
         i.e
           {'clade1': [all, samples, in, clade1.txt], ...}
    '''
    pwd = os.getcwd()
    try:
        meta_data_list = os.listdir(pwd+meta_data)
        #Make a dict of clade name and samples in clade
        #Read in clade_list.txt, make list of sample names
        #strip out white spaces
        #ignore last empty space
        all_samples = {str(sample_list): list(map(lambda s: s.strip(),open(pwd+meta_data+sample_list, 'r').read().split('\n')[:-1])) for sample_list in meta_data_list}
        return all_samples

    except Exception as e:
        print("Error {}, make sure -o input is a direcoty (and ends in a '\\')".format(str(e)))
        exit(-1)

def write_output(transcripts, out_file):
    '''
    Given an iterable of transcrit ids write the results out to outilfe delim
    by newlines
    Inputs:
        transcipts: iterable of transcript ids
        out_file: string for output file path/name
    '''
    with open(out_file, 'w') as f:
        for transcript_id in transcripts:
            f.write("{}\n".format(transcript_id))


def main(pan, meta_data, output):
    meta_data = '/' + meta_data
    if(meta_data[len(meta_data) -1] != '/'):
        meta_data += '/'
    #read in dataframe
    pan_gene_df = pd.read_csv(pan,index_col=0)

    #predicate of interest
    bp = lambda clade_mean: clade_mean == 0

    #Format FPKM out of sample names
    columns = pan_gene_df.columns
    columns = map(lambda x: x.replace('FPKM.', ''), columns)
    pan_gene_df.columns = list(columns)
    
    #parse meta data from meta_data dir
    all_samples = parse_meta_data(meta_data)
    for k,v in all_samples.items():
        print("Searching within clade: {}".format(str(k)))
        match_samples = find_match(pan_gene_df, v)
        print("Within clade  {} found {} matching samples: {}".format(str(k),len(match_samples), str(match_samples) ))
        absent_ids = find_interest(pan_gene_df[match_samples].mean(axis=1), bp)
        print("Found {} interesting transcripts".format(str(len(absent_ids))))
        write_output(absent_ids, str(k) + output)
        print("Wrote output to {}".format(str(k) + output))
    absent_ids = find_interest(pan_gene_df.mean(axis=1), bp)
    print("Searching across all clades")
    print("Found {} interesting transcripts within all clades".format(str(len(absent_ids))))
    write_output(absent_ids, "all_clades" + output)
    #print(pan_gene_df.loc[absent_ids].head())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine and output transcripts of interest(pavs) accross clades", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--counts', required=True, help="CSV file of the transcript counts for all samples")
    parser.add_argument('-o', '--out', help="String to append to all the output file names")
    parser.add_argument('-m', '--meta_data', required=True,  help='Path to Metadata directory \n assumes structre of: \n \
           meta_data \n \
                |\n \
                |-clade1.txt\n \
                |\n \
                |-clade2.txt\n \
where the files contains the sample names for that clade seperated\n \
by newlines\n \
')
    args = parser.parse_args()
    main(args.counts, args.meta_data, args.out)

