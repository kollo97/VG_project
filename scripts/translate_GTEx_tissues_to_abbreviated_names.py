
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description="Translate the GTEx tissues present in GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct to the abbreviated tissue names given in the VG_AE dataset. \nFiles that are read:  ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)         
args  = parser.parse_args()


df = pd.read_csv(r'..\datacollection\output_files\avg_transcript_length_per_tissue_new.tsv', sep='\t', index_col=0)


df2 = pd.read_csv(r'..\raw_data\gtex_sample_attributes.txt', sep='\t')
df2 = df2[df2['SMAFRZE']== 'RNASEQ'].reset_index(drop=True)
df2 = df2.loc[:,['SAMPID','SMTSD']]

abbr =pd.read_excel(r'..\raw_data\VG_AE.xlsx', sheet_name='GTEx Tissue IDs', skiprows=1)
abbr.rename(columns={'TISSUE_NAME_Original':'SMTSD'}, inplace=True)

#VG_AE.xlsx sheet 'GTEX Tissue IDs' contains both the abbreviations and the original tissue names (SMTSD), the gtex_sample_attributes only the original name. Merge dataframes on SMTSD and make a dictionary with the "translation" of SMTSD to abbreviations
c = abbr.merge(df2['SMTSD'],on='SMTSD', how='inner')
c.drop_duplicates(inplace=True)
c.reset_index(inplace=True, drop=True)
c.head()
abbrv_dict = {}
for i, r in c.iterrows():
    abbrv_dict[str(r.SMTSD)] = r.TISSUE_ABBRV
len(abbrv_dict.keys()) #52 tissues are represented in the VG_AE dataset

def diff_list(ls1,ls2):

    set1 = set(ls1)
    set2 = set(ls2)
    diff = sorted(list(set1.symmetric_difference(set2)))
    return diff
diff = diff_list(abbrv_dict.keys(), df.columns) # ['Kidney - Medulla', 'Cells - Cultured fibroblasts'] are in the GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct file (thus in the average_transcript_length_per_tissue_new.tsv) but not in the VG_AE dataset ==> exclude them from this

def shared(ls1,ls2):
    set1 = set(ls1)
    set2 = set(ls2)
    shared = sorted(list(set1.intersection(set2)))

    return shared
shared_tissues = shared(abbrv_dict.keys(), df.columns) 


df = df.loc[:,shared_tissues] # subset the df to use only the shared tissues
df.columns = [abbrv_dict[i] for i in df.columns] # translate the tissues to abbreviation
df.to_csv(r'..\datacollection\output_files\ATLPT_new_abbrv.tsv', sep='\t')


