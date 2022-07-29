import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Reads file GeneHancer_v5.10.gff and returns file enhancer_lengths.csv containing a table with genes, their corresponding genehancer elements, link_score to genehancer element, score of the genehancer element itself, and the length of these elements\n", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
args = parser.parse_args()


df = pd.read_csv(r'..\raw_data\GeneHancer_v5.10.gff', sep='\t') #downloaded after request from GeneCards https://www.genecards.org/Download/File?file=GeneHancer_v5.10.gff
df['score'].astype('float')
df = df[df['score']>=0.7]
df = df.sort_values(by='attributes').reset_index(drop=True) #'attributes' column carries information on the genehancer element ID and multiple associated genes
# print(df.shape)

filter_list = list()
for index, row in df.iterrows():
    filter_list.append(row.attributes.split(';')) #attribute field are separated by semicolon

#  I want to have every gene separately and then its associated genehancer IDs
filter_dict = {'gene':[],'enhancer_id':[],'link_score':[]}
enhancer_ids = [] 
for i,v in enumerate(filter_list): #for every list in filter list, extract the genehancer_ID, the connected genes, and their link_score 
    enhancer_ids.append(filter_list[i][0].strip('genehancer_id=')) #need this later, want to merge the original dataframe with the dataframe generated in the next steps on 'enhancer_id'
    for j in range(0,len(v)-1,2):
        
        gene = str(filter_list[i][j+1]).strip('connected_gene=') #use j+1 because j=0 is the genehancer-ID (could also specify range(1, len(v)))
        genehancer_id = filter_list[i][0].strip('genehancer_id=')
        score = str(filter_list[i][j+2]).strip('score=')
      
        filter_dict['gene'].append(gene)
        filter_dict['enhancer_id'].append(genehancer_id)
        filter_dict['link_score'].append(float(score))

df2 = pd.DataFrame.from_dict(filter_dict)
df2 = df2[df2['link_score']>=7]
num_enhancers = df2.groupby('gene').size() # get the number of enhancers per gene
df2 =  df2.sort_values(by='link_score', ascending=False).groupby('gene', as_index=False).first()#get the best scoring enhancer per gene ==> this enhancer will be used as the representative enhancer length
df2['num_enhancers'] = num_enhancers.values

#now i want to merge the two dataframes again, but with the genes as indices,so just append the enhancer attributes, i.e. start-end ==> so that we can calc length
enhancer_series = pd.Series(enhancer_ids)    
df['enhancer_id'] = enhancer_series #assign the enahcner IDs to a separate column

## MERGE THE TWO DATAFRAMES ==> to calculate the enhancer lengths
concatenated = pd.merge(left=df2, right=df, how='inner', on='enhancer_id') # used inner (probably left join also works) to only keep the rows that are present in the df2
lengths = concatenated['end']-concatenated['start']
concatenated['length'] = pd.Series(lengths)
concatenated = concatenated.loc[:,['gene','num_enhancers','length']]

##WRITE TO FILE
concatenated.to_csv(r'..\raw_data\enhancer_lengths_new.csv',sep='\t',header=1, index=0)
