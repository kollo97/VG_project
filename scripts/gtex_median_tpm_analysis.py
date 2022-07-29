import pandas as pd
df = pd.read_csv(r"..\raw_data\GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", sep='\t', skiprows=2, index_col=0) #contains gene TPM medians per tissue

df.index = [i.split('.')[0] for i in df.index]
g = df.drop('Description', axis=1).median(axis=1).reset_index() #calculate median of tissue medians
g.columns = ['ensembl_gene_id','median_tpm']

g.to_csv(r'..\datacollection\output_files\median_gene_tpm.tsv',sep='\t',header=True,index=False) #gene TPM medians over all tissues
attributes = pd.read_csv(r'..\raw_data\gtex_sample_attributes.txt', sep='\t')
attributes = attributes[attributes['SMAFRZE']== 'RNASEQ'].reset_index(drop=True)
df2 = attributes.loc[:,['SAMPID','SMTSD']]

## check, which tissues are in the VG_AE dataset from the VG paper, and translate those to the original tissue name found in the GTEX sample attributes file
abbr =pd.read_excel(r'..\raw_data\VG_AE.xlsx', sheet_name='GTEx Tissue IDs', skiprows=1)
abbr.rename(columns={'TISSUE_NAME_Original':'SMTSD'}, inplace=True)
c = abbr.merge(df2['SMTSD'],on='SMTSD', how='inner')
c.drop_duplicates(inplace=True)
c.reset_index(inplace=True, drop=True)

#make a dictionary of those translated tissues
abbrv_dict = {}
for i, r in c.iterrows():
    abbrv_dict[str(r.SMTSD)] = r.TISSUE_ABBRV

def shared(ls1,ls2):
    set1 = set(ls1)
    set2 = set(ls2)
    shared = sorted(list(set1.intersection(set2))) #('Kidney - Medulla' and 'Cells - Cultured fibroblasts' are in the GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct file (thus in the average_transcript_length_per_tissue_new.tsv) but not in the VG_AE dataset ==> exclude them from this)

    return shared
shared_tissues = shared(abbrv_dict.keys(), df.columns)

df = df.loc[:,shared_tissues] # subset the df to use only the shared tissues
df.columns = [abbrv_dict[i] for i in df.columns] # translate the tissues to abbreviation
df.index.name = 'ensembl_gene_id'
df.to_csv(r'..\datacollection\output_files\median_gene_tpm_per_tissue.tsv', sep='\t',header=True) # gene TPM medians for every single tissue
