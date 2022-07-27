import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Reads file GeneHancer_v5.10.gff and returns file enhancer_lengths.csv containing a table with genes, their corresponding genehancer elements, link_score to genehancer element, score of the genehancer element itself, and the length of these elements\n", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
args = parser.parse_args()


df = pd.read_csv(r'GeneHancer_v5.10.gff', sep='\t') #downloaded after request from GeneCards https://www.genecards.org/Download/File?file=GeneHancer_v5.10.gff
df['score'].astype('float')
df = df[df['score']>=0.7]
df = df.sort_values(by='attributes').reset_index(drop=True) #'attributes' column carries information on the genehancer element ID and multiple associated genes
# print(df.shape)

filter_list = list()
for index, row in df.iterrows():
    filter_list.append(row.attributes.split(';'))

# I want to have every gene separately and then its associated genehancer IDs
filter_dict = {}
enhancer_ids = [] 
for i,v in enumerate(filter_list):
    enhancer_ids.append(filter_list[i][0].strip('genehancer_id='))
    for j in range(0,len(v)-1,2):
        # print(i,j)
        gene = str(filter_list[i][j+1]).strip('connected_gene=') #use j+1 because j=0 is the genehancer-ID
        genehancer_id = filter_list[i][0].strip('genehancer_id=')
        score = str(filter_list[i][j+2]).strip('score=')
        # print(gene, genehancer_id, score)
        filter_dict[gene] = {'enhancer_id':'','score':float()}
        filter_dict[gene]['enhancer_id'] = genehancer_id
        filter_dict[gene]['score'] = float(score)

filter_df = pd.DataFrame.from_dict(filter_dict, orient='index').reset_index()
filter_df.columns = ['gene','enhancer_id','link_score']

filter_df = filter_df[filter_df['link_score']>=7]
#print(filter_df.head(), filter_df.shape, len(pd.unique(filter_df['gene']))) #soo every gene is associated with one genehancer element ==> shape[0] of filter_df is equal to number of unique genes

#now i want to merge the two dataframes again, but with the genes as indices,so just append the enhancer attributes, i.e. start-end ==> so that we can calc length
enhancer_series = pd.Series(enhancer_ids)    
df['enhancer_id'] = enhancer_series #assign the enahcner IDs to a separate column

## MERGE THE TWO DATAFRAMES ==> to calculate the enhancer lengths

concatenated = pd.merge(left=filter_df, right=df, how='inner', on='enhancer_id') # used inner (probably left join also works) to only keep the rows that are present in the filter_df
lengths = concatenated['end']-concatenated['start']
concatenated['length'] = pd.Series(lengths)
concatenated = concatenated.drop(columns=['attributes','strand','frame','#chrom','source','feature name','start','end'])
##WRITE TO FILE
concatenated.to_csv('enhancer_lengths.csv',sep='\t',header=1, index=0)  # table containing the genes, their corresponding genehancer elements, link_score to genehancer element, score of the genehancer element itself, and the length of these elements
