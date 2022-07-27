#biomart and merging the data

library('biomaRt')
remove(list =ls())
load(file='datacollection.Rdata')
mart = useMart('ensembl',dataset='hsapiens_gene_ensembl')
#listAttributes(mart,what = c("name","description"))[1:80,]
grch38 <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','external_gene_name','external_transcript_name','ensembl_peptide_id','transcript_gencode_basic','chromosome_name','gene_biotype','percentage_gene_gc_content','transcript_length'), mart =mart)

library('dplyr')
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
genehancer_data <- genehancer_data[,c(1,7)]
genehancer_biomart <- grch38
genehancer_biomart <- genehancer_biomart %>%
  left_join(genehancer_data, by=c('ensembl_gene_id'='gene'))%>%
  left_join(genehancer_data, by=c('hgnc_symbol' = 'gene'))%>% 
  left_join(genehancer_data, by=c('external_gene_name' = 'gene'))
#the identifiers overlap, so I'm prioritizing length of genehancer data with 'ensembl_gene_id' (length.x) 
genehancer_biomart[!is.na(genehancer_biomart$length.x),c('length.y','length')] <- NA  # if there is a value for length.x, discard the other length values
genehancer_biomart[!is.na(genehancer_biomart$length.y), 'length'] <- NA #for those values where there was no length.x ==> if theres a value for length.y, discard length.y

library('data.table')
#now there is only one length value per row available ==> make a rowsum and then keep the rowsums values in one single vector/column 
genehancer_biomart$length <- rowSums(genehancer_biomart[, c('length.x','length.y','length')], na.rm = TRUE)
colnames(genehancer_biomart)
genehancer_biomart <- genehancer_biomart[, -c(12,13)] #discard length.x and length.y
colnames(genehancer_biomart)[which(names(genehancer_biomart)=='length')] <- 'enhancer_length'
genehancer_biomart[genehancer_biomart$enhancer_length == 0, 11] <- NA  #turn lengths of 0 into NA 

#median_gene_tpm is fine as it is

raw_final_df <- grch38 %>%
  left_join(avg_vg_ae_data) %>%
  left_join(avg_vg_eqtl_data) %>%
  left_join(genehancer_biomart) %>%
  left_join(RVIS_data) %>%
  left_join(ncGERP_data) %>%
  left_join(string_data) %>%
  left_join(loeuf_data)%>%
  left_join(median_gene_tpm_data)

write.table(raw_final_df, file= 'raw_final_df.txt', sep='\t', row.names = FALSE)

save(file= 'biomart.Rdata',list=ls())

