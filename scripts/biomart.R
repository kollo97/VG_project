#biomart and merging the data

library('biomaRt')
library('dplyr')
library('data.table')

remove(list =ls())
#load(file='Rdata/biomart.Rdata') #all the things generated in this script are in this Rdata




load(file='Rdata/datacollection.Rdata')
mart = useMart('ensembl',dataset='hsapiens_gene_ensembl')
#listAttributes(mart,what = c("name","description"))[1:80,]
grch38 <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','external_gene_name','external_transcript_name','ensembl_peptide_id','transcript_gencode_basic','chromosome_name','gene_biotype','percentage_gene_gc_content','transcript_length'), mart =mart)

##cleaning the column names of the tables a bit, to make joining of tables easier later
colnames(string_data)
colnames(ncGERP_data)[which(names(ncGERP_data)=='CCDS release 15')] <- 'hgnc_symbol'
colnames(ncGERP_data)[which(names(ncGERP_data) %like% 'ncGERP')] <- 'ncGERP'
colnames(ncGERP_data)[which(names(ncGERP_data)%like%'ncCADD')] <- 'ncCADD'

colnames(loeuf_data)[which(names(loeuf_data)=='gene')] <- 'hgnc_symbol'

colnames(avg_vg_ae_data)[which(names(avg_vg_ae_data)=='GeneID')] <- 'ensembl_gene_id'
colnames(avg_vg_ae_data)[which(names(avg_vg_ae_data)=='Avg_VG')] <- 'Avg_VG_AE'

colnames(avg_vg_eqtl_data)[which(names(avg_vg_eqtl_data)=='GeneID')] <- 'ensembl_gene_id'
colnames(avg_vg_eqtl_data)[which(names(avg_vg_eqtl_data)=='Avg_VG')] <- 'Avg_VG_eQTL'

colnames(RVIS_data)[which(names(RVIS_data)=='HGNC gene')] <- 'hgnc_symbol'
colnames(RVIS_data)[which(names(RVIS_data)=='Residual Variation Intolerance Score')] <- 'RVIS_score'


## genehancer had mixed ensembl_gene_id, hgnc_symbol and other gene names all in one column, so I'm trying to join these 
genehancer_biomart <- grch38
genehancer_biomart <- genehancer_biomart %>%
  left_join(genehancer_data, by=c('ensembl_gene_id'='gene'))%>%
  left_join(genehancer_data, by=c('hgnc_symbol' = 'gene'))%>% 
  left_join(genehancer_data, by=c('external_gene_name' = 'gene'))
#the identifiers overlap, so I'm prioritizing length of genehancer data with 'ensembl_gene_id' (length.x) 
genehancer_biomart[!is.na(genehancer_biomart$length.x),c('length.y','length','num_enhancers.y','num_enhancers')] <- NA  # if there is a value for length.x, discard the other length values
genehancer_biomart[!is.na(genehancer_biomart$length.y), c('length','num_enhancers')] <- NA #for those values where there was no length.x ==> if theres a value for length.y, discard length.y

#now there is only one length value per row available ==> make a rowsum and then keep the rowsums values in one single vector/column 
genehancer_biomart$length <- rowSums(genehancer_biomart[, c('length.x','length.y','length')], na.rm = TRUE)
genehancer_biomart$num_enhancers <- rowSums(genehancer_biomart[, c('num_enhancers.x','num_enhancers.y','num_enhancers')], na.rm = TRUE)

#drop the *.x / *.y columns
drop <- which(colnames(genehancer_biomart)%in% c('length.x','length.y','num_enhancers.x','num_enhancers.y'))
genehancer_biomart <- genehancer_biomart[, -drop] #discard length.x, length.y, num_enhancers.x and num_enhancers.y
colnames(genehancer_biomart)[which(names(genehancer_biomart)=='length')] <- 'enhancer_length'

genehancer_biomart[genehancer_biomart$enhancer_length == 0, 'enhancer_length'] <- NA  #turn lengths of 0 into NA

#median_gene_tpm is fine as it is

raw_final_df <- grch38 %>%
  left_join(genehancer_biomart) %>%
  left_join(RVIS_data) %>%
  left_join(ncGERP_data) %>%
  left_join(string_data) %>%
  left_join(loeuf_data)%>%
  left_join(median_gene_tpm_data)%>%
  left_join(avg_vg_eqtl_data) %>%
  left_join(avg_vg_ae_data)
write.table(raw_final_df, file= 'output_files/raw_final_df.txt', sep='\t', row.names = FALSE)

save(file= 'Rdata/biomart.Rdata',list=ls())

