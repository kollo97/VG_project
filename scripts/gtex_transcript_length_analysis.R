library('dplyr')
library('data.table')
library('readxl')
library('tidyverse')
library('ggcorrplot')
library('gridExtra')
load(file = 'gtex_TL_analysis.Rdata') #this file contains the most important stuff generated in this script
avg_vg_ae_data <- read_excel('../raw_data/VG_AE.xlsx', sheet = 1, skip=1)
colnames(avg_vg_ae_data)
avg_vg_ae_data <- avg_vg_ae_data %>%
  rename(ensembl_gene_id=IDs)
avg_vg_ae_data[,2:ncol(avg_vg_ae_data)] <- lapply(avg_vg_ae_data[,2:ncol(avg_vg_ae_data)],as.numeric)
  
avg_vg_eqtl_data <- read_excel('../raw_data/VG_eqtl.xlsx', sheet = 1, skip=1)
avg_vg_eqtl_data <- avg_vg_eqtl_data %>%
  rename(ensembl_gene_id=IDs)
avg_vg_eqtl_data[,2:ncol(avg_vg_eqtl_data)] <- lapply(avg_vg_eqtl_data[,2:ncol(avg_vg_eqtl_data)],as.numeric)


transcript_length <- read.table('../ATLPT_new_abbrv.tsv',sep='\t', header=TRUE)

##CORRELATION OF VG_AE MIT TRANSCRIPT_LENGTH
corr_list <- list()
i <- 1
shared_cols <-intersect(colnames(avg_vg_ae_data)[2:ncol(avg_vg_ae_data)],colnames(transcript_length)[2:ncol(transcript_length)])
for (col in shared_cols ) { #FOR EVERY TISSUE, CALCULATE A CORRELATION VG<==>TRANSCRIPT_LENGTH, AND STORE THE CORRELATION COEFF. IN A LIST
  merged <- avg_vg_ae_data %>%
    select(ensembl_gene_id, col) %>%
    left_join(transcript_length%>%select(ensembl_gene_id,col), by = 'ensembl_gene_id')
  print(colnames(merged))
  corr2 <- cor(merged[,2],merged[,3], use = 'pairwise.complete.obs', method= 'spearman') #SPEARMAN
  corr_list[[i]]<-corr2
  i <- i+1
  
}

names(corr_list) <- shared_cols
corr_list<- corr_list[order(unlist(corr_list), decreasing=TRUE)]
new_df <- as.data.frame(corr_list, row.names = 'transcript_length') #TABLE WITH ALL THE CORRELATIONS

#MAKE A GGCORRPLOT, BUT TAKE ONLY ~10 VALUES (SO THAT IT'S EASIER TO VISUALIZE). LATER SHOW THE MULTIPLE PLOTS IN A GRID 
plist <- list()
j <- 1
for (i in seq(from = 10, to = ncol(new_df), length.out=6)){
  plist[[j]]<-ggcorrplot(t(new_df[,c((i-9):i)]),outline.col = "white", lab=TRUE)
  j <- j+1
}

do.call('grid.arrange', c(plist, nrow=6))


##SAME FOR EQTL ==> CORRELATION OF VG_eQTL MIT TRANSCRIPT_LENGTH
corr_list <- list()
i <- 1
shared_cols <-intersect(colnames(avg_vg_eqtl_data)[2:ncol(avg_vg_eqtl_data)],colnames(transcript_length)[2:ncol(transcript_length)])
for (col in shared_cols ) {
  merged <- avg_vg_eqtl_data %>%
    select(ensembl_gene_id, col) %>%
    left_join(transcript_length%>%select(ensembl_gene_id,col), by = 'ensembl_gene_id')
  print(colnames(merged))
  corr2 <- cor(merged[,2],merged[,3], use = 'pairwise.complete.obs', method= 'spearman')
  corr_list[[i]]<-corr2
  i <- i+1
  
}

names(corr_list) <- shared_cols
corr_list<- corr_list[order(unlist(corr_list), decreasing=TRUE)]
new_df <- as.data.frame(corr_list, row.names = 'transcript_length')

#MAKE A GGCORRPLOT, BUT TAKE ONLY ~10 VALUES (SO THAT IT'S EASIER TO VISUALIZE). LATER SHOW THE MULTIPLE PLOTS IN A GRID 
plist <- list()
j <- 1
for (i in seq(from = 10, to = ncol(new_df), length.out=6)){
  plist[[j]]<-ggcorrplot(t(new_df[,c((i-9):i)]),outline.col = "white", lab=TRUE)
  j <- j+1
}
library('gridExtra')
do.call('grid.arrange', c(plist, nrow=6))


















map(avg_vg_ae_data, ~sum(is.na(.))/nrow(avg_vg_ae_data)) #MANY NAs IN EVERY TISSUE 
map(avg_vg_eqtl_data, ~sum(is.na(.))/nrow(avg_vg_eqtl_data)) #MANY NAs IN EVERY TISSUE 

##NOW LET'S ANALYSE THE CORRELATION OF TRANSCRIPT LENGTH AND VG WITHIN ONE GENE OVER MULTIPLE TISSUES
#PROBLEM ==> FOR MANY GENES ONLY A FEW TISSUES AVAILABLE ==> might filter those with only one (or two/three?==> correlation of two tissues pretty weak) tissue available out
rowsums <- rowSums(!is.na(avg_vg_ae_data))
hist(rowsums, breaks = 50)
avg_vg_ae_data$rowsums <- rowsums

ae_gene_corr <- avg_vg_ae_data %>%
  filter(rowsums>10)
test <- as.data.frame(t(ae_gene_corr))
names(test) <- test[1,]
test <- test[-1,]
test <- test %>%
  rownames_to_column(var = 'tissue')
test[,2:ncol(test)] <- lapply(test[,2:ncol(test)],as.numeric)

test2 <- as.data.frame(t(transcript_length))
names(test2) <- test2[1,]
test2 <- test2[-1,]
test2 <- test2 %>%
  rownames_to_column(var = 'tissue')
test2[,2:ncol(test2)] <- lapply(test2[,2:ncol(test2)],as.numeric)

shared_genes <- intersect(colnames(test)[2:ncol(test)],colnames(test2)[2:ncol(test2)])

corr_list <- list()
i <- 1
#rand_col <- sample(1:length(shared_genes), length(shared_genes),replace=T )

##CORRELATE VG OF ALL GENES OVER ALL TISSUES WITH TRANSCRIPT LENGTH OF GENE OVER ALL TISSUES
for (col in shared_genes ) {
  df <- test2[,c('tissue',col)]
  
  merged <- test %>%
    select(tissue, col) %>%
    left_join(df, by = 'tissue')
  #print(colnames(merged))
  corr2 <- cor(merged[,2],merged[,3], use = 'pairwise.complete.obs', method= 'spearman')
  corr_list[[i]]<-corr2
  i <- i+1
  if ((i/length(shared_genes)*100)%%5==0){
    print(c(i/length(shared_genes),' %'))
  }
  
    
}
names(corr_list) <- shared_genes
corrs <- corr_list[order(unlist(corr_list), decreasing=TRUE)]
new_df <- as.data.frame(corr_list, row.names = 'transcript_length')
new_df <- as.data.frame(t(new_df))
#new_df[sort(abs(new_df$transcript_length))]
nums <- c(1:nrow(new_df))
new_df$genes <- nums
group <- rep.int(1, nrow(new_df))
new_df$group <- group
new_df$group <- as.factor(new_df$group)

plot(new_df$genes, new_df$transcript_length)
library('ggplot2')
p <- ggplot(new_df, aes(x=group,y=transcript_length))+
  geom_violin(trim=FALSE)+
  geom_boxplot()
p
#save(file = 'gtex_TL_analysis.Rdata', list= ls())


##PLAYING AROUND, NOT SO IMPORTANT
#TEST FOR RANDOM CORRELATION ==> SHUFFLE CORRELATION
corr_list <- list()
i <- 1
for (col in shared_genes ) {
  rand <-sample(1:length(shared_genes),1)
  df <- test2[,c('tissue',shared_genes[rand])]
  #print(c(rand,df))
  merged <- test %>%
    select(tissue, col) %>%
    left_join(df, by = 'tissue')
  #print(colnames(merged),df)
  corr2 <- cor(merged[,2],merged[,3], use = 'pairwise.complete.obs', method= 'spearman')
  corr_list[[i]]<-corr2
  i <- i+1
  if ((i/length(shared_genes)*100)%%5==0){
    print(i/length(shared_genes))
  }
  
  
}
names(corr_list) <- shared_genes
corrs <- corr_list[order(unlist(corr_list), decreasing=TRUE)]
new_df <- as.data.frame(corr_list, row.names = 'transcript_length')
new_df <- as.data.frame(t(new_df))
#new_df[sort(abs(new_df$transcript_length))]
nums <- c(1:nrow(new_df))
new_df$genes <- nums
group <- rep.int(1, nrow(new_df))
new_df$group <- group
new_df$group <- as.factor(new_df$group)

plot(new_df$genes, new_df$transcript_length)
library('ggplot2')
library('ggsignif')
p <- ggplot(new_df, aes(x=group,y=transcript_length))+
  geom_violin(trim=FALSE)+
  geom_boxplot()
p


##GENERATION OF BINS FOR TRANSCRIPT LENGTH. FOR EACH TISSUE, DEFINE THE BIN 
#FOR AE DATA
shared_cols <-intersect(colnames(avg_vg_ae_data)[2:ncol(avg_vg_ae_data)],colnames(transcript_length)[2:ncol(transcript_length)])
test3 <- transcript_length %>%
  filter(ensembl_gene_id %in% avg_vg_ae_data$ensembl_gene_id)%>% 
  select(shared_cols) 
test3 <- data.frame(transcript_length = unlist(test3))
bins <- c(quantile(test3$transcript_length, probs = seq(0.05, 1, by = 0.05)))
test3$TL_bin <- cut(test3$transcript_length, breaks=bins, na.rm=TRUE, labels = c(1:19))
                
test3 <- test3%>%
  rownames_to_column(var='tissues')
test4 <- avg_vg_ae_data %>%
  select(shared_cols)
test4 <- data.frame(VG = unlist(test4))

test4 <- test4%>%
  rownames_to_column(var='tissues')

merged <- test4 %>%
  left_join(test3, by = 'tissues')%>%
  drop_na()
group <- rep.int(1, nrow(merged))
merged$group <- group
merged$group <- as.factor(merged$group)
merged$sdg <- sqrt(merged$VG)
library(ggbeeswarm)
quantil_95 <- quantile(merged$sdg, probs = 0.98)
p <- ggplot(merged, aes(x=TL_bin,y=sdg))+
  geom_boxplot()+
  ylim(0,quantil_95)
p

#FOR eQTL DATA 
shared_cols <-intersect(colnames(avg_vg_eqtl_data)[2:ncol(avg_vg_eqtl_data)],colnames(transcript_length)[2:ncol(transcript_length)])
test3 <- transcript_length %>%
  filter(ensembl_gene_id %in% avg_vg_eqtl_data$ensembl_gene_id)%>% 
  select(shared_cols) 
test3 <- data.frame(transcript_length = unlist(test3))
bins <- c(quantile(test3$transcript_length, probs = seq(0.20, 1, by = 0.05)))
test3$TL_bin <- cut(test3$transcript_length, breaks=bins, na.rm=TRUE, labels = c(1:16))
test3 <- test3%>%
  rownames_to_column(var='tissues')
test4 <- avg_vg_eqtl_data %>%
  select(shared_cols)
test4 <- data.frame(VG = unlist(test4))

test4 <- test4%>%
  rownames_to_column(var='tissues')

merged <- test4 %>%
  left_join(test3, by = 'tissues')%>%
  drop_na()
group <- rep.int(1, nrow(merged))
merged$group <- group
merged$group <- as.factor(merged$group)
merged$sdg <- sqrt(merged$VG)
quantil_95 <- quantile(merged$sdg, probs = 0.98)
p <- ggplot(merged, aes(x=TL_bin,y=sdg))+
  geom_boxplot()+
  ylim(0,quantil_95)
p
