import pandas as pd
import numpy as np
import itertools
import argparse
parser = argparse.ArgumentParser(description="Reads GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct and returns avergae transcript length of a gen per tissue (avg_transcript_length_per_tissue.tsv), and median variance of transcript expression per gen (transcript_lengths.tsv) \nFor the transcript lengths you need a file containing ensembl_gene_id ensembl_transcript_id and the transcript lenghts as columns (e.g. from biomaRt, all grch38 ensembl 106 genes)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)         
args  = parser.parse_args()

def presettings(gtex_file = 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct',sample_attributes_file = 'gtex_sample_attributes.txt'):
    cols = pd.read_csv(gtex_file,nrows=0,skiprows=2, sep='\t' ) #just read in the column name basically
    columns = list(cols.columns)
    dtypes = (('float ')*(cols.shape[1]-2)).split() #just make a list of multiple times string 'float' (for the GTEX-SAMPLE_ID columns containing all the TPMs)
    dtypes.insert(0,'object')
    dtypes.insert(0,'object') #for the gene/transcript_id columns

    col_dtypes = dict(itertools.zip_longest(columns,dtypes)) #a dictionary of col types that pd.read_csv can take
    df = pd.read_csv(sample_attributes_file, sep='\t')
    df = df[df['SMAFRZE']== 'RNASEQ'].reset_index(drop=True)
    df2 = df.loc[:,['SAMPID','SMTSD']]
    df2 =  df2[df2['SAMPID'].isin(columns[2:])]

    sampid_dict = df2.groupby(by='SAMPID').groups
    transl_sampid_dict = {}
    for k,v in sampid_dict.items():
        transl_sampid_dict[k] = df2.iloc[v].SMTSD.values[0]

    return col_dtypes, transl_sampid_dict


def get_tpm_medians_per_tissue(column_dtypes,tr_sampid_dict , gtex_file = "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct" ):
    ##realized that pandas might not be so perfect for this dataset, but in order to make it work faster, I want to specify the column dtypes explicitly for each column (I need a dictionary according to read_csv documnetation)
    ## two things to make pd.read_csv faster: specify column dtypes, use chunksize and then concat
     
    dfs = pd.read_csv(gtex_file, dtype=column_dtypes,skiprows=2, sep='\t', chunksize=100 )
    medians_list = []
    for df in dfs:
        
        df.loc[:,'gene_id'] = df.loc[:,'gene_id'].apply(lambda x: x.split('.')[0])#
        df.loc[:,'transcript_id'] = df.loc[:,'transcript_id'].apply(lambda x: x.split('.')[0])#
        
        df.set_index(['gene_id','transcript_id'], inplace=True)
        df.columns = pd.MultiIndex.from_tuples([(i.split('-')[0], tr_sampid_dict[i]) for i in df.columns],names=['donor_id','tissue_name'])
        
        medians = df.groupby(level=1, axis=1).median() #median TPM expression per tissue for every transcript
        medians_list.append(medians)
    medians = pd.concat(medians_list)

  
    return medians
    


def get_avg_transcript_length(median_df):
   
    df2 = pd.read_csv('transcript_lengths.tsv', sep='\t', index_col=[0,1]) #previously generated using biomaRt in R ==> all grch38 ensembl 106 genes, their gene_ids, transcript_ids, and transcript_length
    df2 = df2.sort_index()
    
    g = median_df.groupby(level=0, axis=0) #groupby gene_id (first level of ROW multiindex)
    g = g.transform(lambda x: x/sum(x)) #https://sparkbyexamples.com/python/pandas-percentage-total-with-groupby/ divide every TPM of the transcripts of a gene by the sum of TPM ==>fraction_tpm_per_transcript
  
    concat =g.rename_axis(['ensembl_gene_id','ensembl_transcript_id']).join(df2) #merge the two dataframes (on multiindex 'ensembl_gene_id','ensembl_transcript_id')
    
    concat = concat.iloc[:,:-1].mul(concat.transcript_length, axis=0) #multiply every column with the last column, which is the transcript length
    
    g = concat.groupby(level=0, axis=0).sum() #==> group by gene id and calculate sum of fractions
    # g.reset_index(inplace=True) #reset index to write it as a normal column to csv later
    
    return g


def main():
    col_dtypes,tr_sampid_dict = presettings()

    med_df = get_tpm_medians_per_tissue(col_dtypes,tr_sampid_dict)
    
    avg_transcript_length = get_avg_transcript_length(med_df)
    avg_transcript_length.to_csv('avg_transcript_length_per_tissue.tsv',sep='\t')

if __name__ == "__main__":
    main()