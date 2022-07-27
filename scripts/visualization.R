#visualization
library('dplyr')
library('tidyverse')
remove(list=ls())

df1 <- read.table('raw_final_df.txt',sep='\t',header=1)
df <- df1
colnames(df)
colnames(df)[which(names(df)=='oe_lof_upper')] <- 'loeuf_score'
load(file='visualization.Rdata')

##I WANT TO FILTER THE TABLE TO ULTIMATELY ONLY HAVE ON ROW PER ENSEMBL_GENE_ID
#filtering is done based on how much external information we have in each row ==> represented in the "rowsums" column
#if there are multiple rows with same external information==> keep longest transcript
#if there are still duplicate ensembl_gene_ids ==> distinct() ==> removes duplicates 
df$rowsums <- rowSums(!is.na(df[,c(12:ncol(df))]))
df <- df[df$rowsums>=1,] #select only the rows where you have data from the external metrics ==> 0 means, there are 0 values available
df <- df[!duplicated(df),]
df <- df %>%
  group_by(ensembl_gene_id)%>%
  filter(rowsums == max(rowsums)) %>%
  filter(!(rowsums == 1 & !is.na(median_tpm))) %>% #there were so many genes that had only the median_tpm values that I decided to get rid of those
  filter(transcript_length == max(transcript_length))
df <- distinct(df, ensembl_gene_id, .keep_all=TRUE)
df <- df %>% ungroup()

#ncGERP data has some '-' strings in some fields instead of NA ==> clean up
test <- df %>% select_if(negate(is.numeric))
test <-test[is.na(as.numeric(as.character(test$ncGERP))),]
test <- test[!is.na(test$ncGERP),]
df[df$ensembl_gene_id %in%test$ensembl_gene_id, c('ncGERP','ncCADD')] <- NA

## KEEP THE RAW DF AND MAKE A DF WITH ONLY THE NUMERIC METRICS, NO IDENTIFIERS OR GENE_BIOTYPE OR CHROMOSOMAL LOCATION OF THE GENE
for_plot <- df[, c(10:24)]
for_plot[names(for_plot)] <- lapply(for_plot[names(for_plot)], as.numeric)
for_plot$sdg_AE <- sqrt(for_plot$Avg_VG_AE)     #sdg has been often used in the VG paper so I included it too to see if I get similar results. Not relevant for spearman correlation but anyways..
for_plot$sdg_eQTL <- sqrt(for_plot$Avg_VG_eQTL)
lapply(for_plot, class)
colnames(for_plot)
for_plot <- for_plot[,c(1,2,5:15,4,3,17,16)] #just reordering the columns, this way the ggcorrplot looks tidier, easier to focus on the sdg_AE in the last column

##CORRPLOT
library('ggplot2')
library('ggcorrplot')
corr <- cor(for_plot, use = 'pairwise.complete.obs', method= 'spearman') #full corrplot matrix
corr2 <- rbind(data.frame(sdg_AE = 0),data.frame(sdg_AE = corr[,15])) #only sdg_AE correlations to every other variable 

ggc <- ggcorrplot(corr, type = "lower",
           outline.col = "white",
           lab=TRUE)
ggc

#write.table(for_plot,file='for_plot.tsv',sep='\t',row.names = FALSE)

##VENNDIAGRAM, to check for which genes we have data for each metric/variable 
#make vectors of row indices/ row names and then use ggvenndiagram on these row indices
vg_aes <- df %>%
  rownames_to_column()%>%
  mutate(rowname = as.numeric(rowname))%>%
  filter(!is.na(Avg_VG_AE))%>%
  pull(rowname)
vg_eqtls <- df %>%
  rownames_to_column()%>%
  mutate(rowname = as.numeric(rowname))%>%
  filter(!is.na(Avg_VG_eQTL))%>%
  pull(rowname)
loeufs<- df %>%
  rownames_to_column()%>%
  mutate(rowname = as.numeric(rowname))%>%
  filter(!is.na(loeuf_score))%>%
  pull(rowname)
rviss<- df %>%
  rownames_to_column()%>%
  mutate(rowname = as.numeric(rowname))%>%
  filter(!is.na(RVIS_score))%>%
  pull(rowname)
ncgerps<- df %>%
  rownames_to_column()%>%
  mutate(rowname = as.numeric(rowname))%>%
  filter(!is.na(ncGERP))%>%
  pull(rowname)
ncrviss<- df %>%
  rownames_to_column()%>%
  mutate(rowname = as.numeric(rowname))%>%
  filter(!is.na(ncRVIS))%>%
  pull(rowname)
strings<- df %>%
  rownames_to_column()%>%
  mutate(rowname = as.numeric(rowname))%>%
  filter(!is.na(degree))%>%
  pull(rowname)
genehancers <- df %>%
  rownames_to_column()%>%
  mutate(rowname = as.numeric(rowname))%>%
  filter(!is.na(enhancer_length))%>%
  pull(rowname)
tpms <- df %>%
  rownames_to_column()%>%
  mutate(rowname = as.numeric(rowname))%>%
  filter(!is.na(median_tpm))%>%
  pull(rowname)

library('ggVennDiagram')
x <- list( vg_AE = vg_aes,  loeuf_score = loeufs, ncGERP = ncgerps,ncRVIS=ncrviss, RVIS = rviss, TPM = tpms)
#loeuf_score = loeufs, ncGERP = ncgerps, RVIS = rviss, vg_AE = vg_aes, vg_eQTL = vg_eqtls, STRING = strings,vg_eQTL = vg_eqtls, TPM = tpms
ggVennDiagram(x, label='count')



##PLAYING AROUND, NOT SOO IMPORTANT
#histogram with loeuf_score und is.na(VG_AE) als count bzw fill ==> for which genes, based on their LOEUF score, do we have VG data
breaks = seq(0, round(max(for_plot$loeuf_score, na.rm = TRUE)), length.out = 11)
p <- for_plot[loeufs,] %>% 
  ggplot(aes(loeuf_score, fill=!is.na(Avg_VG_AE), label = ..count..))+
  geom_histogram(breaks = breaks)+
  theme_classic()+
  geom_text(stat='bin', size=5,vjust=0, breaks = breaks, position='stack')

p

save(file='visualization.Rdata',list=ls())
