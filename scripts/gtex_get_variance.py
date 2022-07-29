## RUNS WAY TOO LONG, DIDN'T FINISH WITHIN 24H ON 4 NODES ON BIANCA. NOT REALLY FEASIBLE. NEEDS TO BE IMPROVED
import pandas as pd
import numpy as np
import itertools
import argparse
parser = argparse.ArgumentParser(description="Reads GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct and returns median variance of median gene length per tissue (var_length_per_tissue.tsv)\nFor the transcript lengths you need a file containing ensembl_gene_id ensembl_transcript_id and the transcript lenghts as columns (e.g. from biomaRt, all grch38 ensembl 106 genes) \nRequired files in the same directory: \nGTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct \ttranscript_lengths.tsv", formatter_class=argparse.ArgumentDefaultsHelpFormatter)         
args  = parser.parse_args()

def presettings(gtex_file = r'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct',sample_attributes_file = 'gtex_sample_attributes.txt'):
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



def get_avg_transcript_length(col_dtypes, gtex_file=r'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct'):
   
    df2 = pd.read_csv('transcript_lengths.tsv', sep='\t', index_col=[0,1]) #previously generated using biomaRt in R ==> all grch38 ensembl 106 genes, their gene_ids, transcript_ids, and transcript_length
    df2 = df2.sort_index()
    dfs = pd.read_csv(gtex_file ,dtype=col_dtypes,skiprows=2, sep='\t', nrows=100, chunksize=20 )
    df_list = []
    for df in dfs:
        df.loc[:,'gene_id'] = df.loc[:,'gene_id'].apply(lambda x: x.split('.')[0])#
        df.loc[:,'transcript_id'] = df.loc[:,'transcript_id'].apply(lambda x: x.split('.')[0])
        df.set_index(['gene_id','transcript_id'], inplace=True)
        
        
        g = df.groupby(level=0, axis=0) #groupby gene_id (first level of ROW multiindex)
        g = g.transform(lambda x: x/sum(x)) #https://sparkbyexamples.com/python/pandas-percentage-total-with-groupby/ divide every TPM of the transcripts of a gene by the sum of TPM ==>fraction_tpm_per_transcript
       
        concat =g.rename_axis(['ensembl_gene_id','ensembl_transcript_id']).join(df2) #merge the two dataframes (on multiindex 'ensembl_gene_id','ensembl_transcript_id')
        
        concat = concat.iloc[:,:-1].mul(concat.transcript_length, axis=0) #multiply every column with the last column, which is the transcript length
        
        g = concat.groupby(level=0, axis=0).sum() #==> group by gene id and calculate sum of fractions
        df_list.append(g)
    median_length_per_subject = pd.concat(df_list)
    return median_length_per_subject


def get_length_var_per_tissue(tr_sampid_dict, df):      
    
    df.columns = pd.MultiIndex.from_tuples([(i.split('-')[0], tr_sampid_dict[i]) for i in df.columns],names=['donor_id','tissue_name'])
    
    var = df.groupby(level=1, axis=1).var(ddof=0) 
  
    return var

def main():
    col_dtypes,tr_sampid_dict = presettings()

    med_df = get_avg_transcript_length(col_dtypes)
    
    var_length_per_tissue = get_length_var_per_tissue(tr_sampid_dict,med_df)
    var_length_per_tissue.to_csv('var_length_per_tissue.tsv',sep='\t')

if __name__ == "__main__":
    main()