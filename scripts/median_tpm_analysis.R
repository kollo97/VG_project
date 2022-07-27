library('dplyr')
library('data.table')
library('readxl')
library('tidyverse')
library('ggcorrplot')
library('gridExtra')
remove(list=ls())
avg_vg_ae_data <- read_excel('../raw_data/VG_AE.xlsx', sheet = 1, skip=1)
colnames(avg_vg_ae_data)
avg_vg_ae_data <- avg_vg_ae_data %>%
  rename(ensembl_gene_id=IDs)
avg_vg_ae_data[,2:ncol(avg_vg_ae_data)] <- lapply(avg_vg_ae_data[,2:ncol(avg_vg_ae_data)],as.numeric)

avg_vg_eqtl_data <- read_excel('../raw_data/VG_eqtl.xlsx', sheet = 1, skip=1)
avg_vg_eqtl_data <- avg_vg_eqtl_data %>%
  rename(ensembl_gene_id=IDs)
avg_vg_eqtl_data[,2:ncol(avg_vg_eqtl_data)] <- lapply(avg_vg_eqtl_data[,2:ncol(avg_vg_eqtl_data)],as.numeric)

median_tpm_per_tissue_data <- read.table('../raw_data/median_gene_tpm_per_tissue.tsv', sep='\t', header=T, check.names=F)  
##CORRELATION OF VG_AE MIT TRANSCRIPT_LENGTH
corr_list <- list()
i <- 1
shared_cols <-intersect(colnames(avg_vg_ae_data)[2:ncol(avg_vg_ae_data)],colnames(median_tpm_per_tissue_data)[2:ncol(median_tpm_per_tissue_data)])
for (col in shared_cols ) { #FOR EVERY TISSUE, CALCULATE A CORRELATION VG<==>TRANSCRIPT_LENGTH, AND STORE THE CORRELATION COEFF. IN A LIST
  merged <- avg_vg_ae_data %>%
    select(ensembl_gene_id, col) %>%
    left_join(median_tpm_per_tissue_data%>%select(ensembl_gene_id,col), by = 'ensembl_gene_id')
  print(colnames(merged))
  corr2 <- cor(merged[,2],merged[,3], use = 'pairwise.complete.obs', method= 'spearman') #SPEARMAN
  corr_list[[i]]<-corr2
  i <- i+1
  
}

names(corr_list) <- shared_cols
corr_list<- corr_list[order(unlist(corr_list), decreasing=TRUE)]
new_df <- as.data.frame(corr_list, row.names = 'median_tpm') #TABLE WITH ALL THE CORRELATIONS

plist <- list()
j <- 1
for (i in seq(from = 10, to = ncol(new_df), length.out=6)){
  plist[[j]]<-ggcorrplot(t(new_df[,c((i-9):i)]),outline.col = "white", lab=TRUE)
  j <- j+1
}

do.call('grid.arrange', c(plist, nrow=6))

##SAME FOR eQTL DATA 
##CORRELATION OF VG_AE MIT TRANSCRIPT_LENGTH
corr_list <- list()
i <- 1
shared_cols <-intersect(colnames(avg_vg_eqtl_data)[2:ncol(avg_vg_eqtl_data)],colnames(median_tpm_per_tissue_data)[2:ncol(median_tpm_per_tissue_data)])
for (col in shared_cols ) { #FOR EVERY TISSUE, CALCULATE A CORRELATION VG<==>TRANSCRIPT_LENGTH, AND STORE THE CORRELATION COEFF. IN A LIST
  merged <- avg_vg_eqtl_data %>%
    select(ensembl_gene_id, col) %>%
    left_join(median_tpm_per_tissue_data%>%select(ensembl_gene_id,col), by = 'ensembl_gene_id')
  print(colnames(merged))
  corr2 <- cor(merged[,2],merged[,3], use = 'pairwise.complete.obs', method= 'spearman') #SPEARMAN
  corr_list[[i]]<-corr2
  i <- i+1
  
}

names(corr_list) <- shared_cols
corr_list<- corr_list[order(unlist(corr_list), decreasing=TRUE)]
new_df <- as.data.frame(corr_list, row.names = 'median_tpm') #TABLE WITH ALL THE CORRELATIONS

plist <- list()
j <- 1
for (i in seq(from = 10, to = ncol(new_df), length.out=6)){
  plist[[j]]<-ggcorrplot(t(new_df[,c((i-9):i)]),outline.col = "white", lab=TRUE)
  j <- j+1
}

do.call('grid.arrange', c(plist, nrow=6))


