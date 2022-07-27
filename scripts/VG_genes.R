#VG_genes
remove(list=ls())
load(file='visualization.Rdata')
load(file='VG_genes.Rdata')
df$sdg_AE <- for_plot$sdg_AE
df$sdg_eQTL <- for_plot$sdg_eQTL
df <- rownames_to_column(df)
ae_eqtls_comp <- df %>%
  select(ensembl_gene_id,rowname, sdg_AE, sdg_eQTL)%>%
  filter(!is.na(sdg_AE)&!is.na(sdg_eQTL))

ae_eqtl.lm = lm(sdg_AE ~ sdg_eQTL, data=ae_eqtls_comp) 
ae_eqtl.res = resid(ae_eqtl.lm)


ae_eqtls_comp$residuals <- abs(ae_eqtl.res)
ae_eqtls_comp$resid_bin <- cut(ae_eqtls_comp$residuals, breaks=c(quantile(ae_eqtls_comp$residuals, probs = seq(0, 1, by = 0.10), na.rm=TRUE)),
                               labels=c(1:10))
ae_eqtls_comp$resid_bin <- as.factor(ae_eqtls_comp$resid_bin)


p <- ae_eqtls_comp%>% 
  ggplot(aes(x=sdg_AE, y=sdg_eQTL, col=resid_bin))+
  geom_point()+
  xlim(0,0.5)+
  ylim(0,0.5)+
  scale_color_brewer(palette="Set3")

p



ae_eqtls_comp$AE_bin <- cut(ae_eqtls_comp$sdg_AE, breaks=c(quantile(ae_eqtls_comp$sdg_AE, probs = seq(0, 1, by = 0.10), na.rm=TRUE)),
                        labels=c(1:10)
                          )
ae_eqtls_comp$eQTL_bin <- cut(ae_eqtls_comp$sdg_eQTL, breaks=c(quantile(ae_eqtls_comp$sdg_eQTL, probs = seq(0, 1, by = 0.10), na.rm=TRUE)),
                            labels=c(1:10)
                            )
colnames(ae_eqtls_comp)
ae_eqtls_comp[names(ae_eqtls_comp)][7:8] <-lapply(ae_eqtls_comp[names(ae_eqtls_comp)][7:8], as.numeric) 
ae_eqtls_comp$difference <- ae_eqtls_comp$AE_bin-ae_eqtls_comp$eQTL_bin
ae_eqtls_comp$difference <-as.factor(ae_eqtls_comp$difference) 


breaks = seq(0, 10,1)
breks <- as.character(breaks)
p <- ae_eqtls_comp%>% 
  ggplot(aes(difference, label =..count..))+
  geom_bar()+
  theme_classic()+
  geom_text( stat = 'count', size=5,vjust=0,  position='stack')

p







ae_eqtls_comp[is.na(ae_eqtls_comp$difference),]

ae_eqtls_comp %>%
  ggplot(aes(x=sdg_AE, y=sdg_eQTL, col=difference))+
  geom_point()+
  xlim(0,0.5)+
  ylim(0,0.5)+
  scale_color_brewer(palette="Set3")
lmTemp <- lm(sdg_AE~sdg_eQTL, ae_eqtls_comp)
plot(sqrt(lmTemp$residuals))  

ae_eqtls_comp$difference <- as.numeric(ae_eqtls_comp$difference)
big_diff <- ae_eqtls_comp %>%
  filter(as.numeric(difference) >7)%>%
  pull(rowname)
big_diff<-as.numeric(big_diff)

big_diff <- for_plot[big_diff,]

small_diff <- ae_eqtls_comp %>%
  filter(as.numeric(difference) <(-7))%>%
  pull(rowname)
small_diff<-as.numeric(small_diff)

small_diff <- for_plot[small_diff,]

medium_diff <- ae_eqtls_comp %>%
  filter(4< as.numeric(difference) & as.numeric(difference) < 6) %>%
  pull(rowname)
medium_diff<-as.numeric(medium_diff)

medium_diff <- for_plot[medium_diff,]

corr <- cor(small_diff, use = 'pairwise.complete.obs', method= 'spearman')
p.mat <- cor_pmat(small_diff)
#head(p.mat)
ggcorrplot(corr, type = "lower",
           outline.col = "white",
           lab=TRUE)+
  ggtitle('different in quantile between 4 and 6')

save(file='VG_genes.Rdata', list=ls())
