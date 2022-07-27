library('readxl')
library('data.table')
library('dplyr')
#remove(list=ls())
avg_vg_ae_data <- read_excel('../raw_data/VG_AE.xlsx', sheet = 3, skip=1)
avg_vg_eqtl_data <- read_excel('../raw_data/VG_eqtl.xlsx', sheet = 2, skip=1)

genehancer_data <- read.table('../raw_data/enhancer_lengths.csv', sep='\t', header=1)

RVIS_data <- read_excel('../raw_data/RVIS_per_gene.xlsx') %>% .[,-3]

ncGERP_data <- read_excel('../raw_data/ncGERP.xlsx') %>% .[,c(2,5,7,13)]

library('igraph')
string_data <- read.table('../raw_data/string_db_filtered.txt', sep='\t', header=1)
string_data[,1] <- gsub('9606.','', string_data[,1])
string_data[,2] <- gsub('9606.','', string_data[,2])
string_graph <- graph_from_data_frame(string_data, directed = FALSE, vertices = NULL)
degrees <- degree(string_graph, v =V(string_graph))
local_corr_coeff <-   transitivity(string_graph,type = "local", vids = V(string_graph))
string_data <- data.frame(ensembl_peptide_id= names(degrees), degree = degrees, transitivity = local_corr_coeff)

loeuf_data <- read.table('../raw_data/loeuf_per_gene.txt', sep='\t', header=T) %>% .[,c('gene','oe_lof_upper','pLI','cds_length')]

median_tpm_per_tissue_data <- read.table('../raw_data/median_gene_tpm_per_tissue.tsv', sep='\t', header=T, check.names=F) %>% 
  select(!Description)
median_gene_tpm_data <- read.table('../raw_data/median_gene_tpm.tsv', sep='\t', header=T, check.names=F)


save(file='datacollection.Rdata', list=ls())
