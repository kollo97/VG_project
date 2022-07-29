library('randomForest')
library('caret')
library('tidyverse')
library('dplyr')
remove(list=ls())
load(file='Rdata/supervised.Rdata') #all the things generated in this script are in this Rdata 


#load(file='Rdata/visualization.Rdata')
superv_df <- df %>%
  select(ensembl_gene_id,Avg_VG_AE, loeuf_score, ncGERP, ncRVIS,  RVIS_score,  Avg_VG_eQTL,  median_tpm, num_enhancers)%>% #with or without num_enhancers, does not make a big differences
  drop_na()%>%
  column_to_rownames(var='ensembl_gene_id')
#testing correlation on this subset of genes ==> cor(superv_df$median_tpm, superv_df$transcript_length, method = 'spearman')  #==> only -0.034, so the "correlation" seen when you compare all genes in the raw df containing all genes is not 

#possible variables: transcript_length, percentage_gene_gc_content, loeuf_score, ncGERP, ncRVIS, transitivity, degree, RVIS_score  
#,percentage_gene_gc_content, transcript_length, transitivity
##ADD LOG AE_VG and eQTL_VG, "exchange" it for VG
superv_df$log_AE <- log(superv_df$Avg_VG_AE)
superv_df$log_eQTL <- log(superv_df$Avg_VG_eQTL)
subset <- superv_df %>%
  select(!c(contains('VG')))%>% 
  drop_na()#%>%

#map(superv_df, ~var(.))
set.seed(123)
inTrain <- createDataPartition(y = subset$log_AE, p=0.65, list =F)
training <- subset[inTrain,]%>%
  drop_na()
inValidation <- createDataPartition(y = training$log_AE, p=0.65, list =F)
training <- training[inValidation,]%>%
  drop_na()
validation <- training[-inValidation,]%>%
  drop_na()

test <- subset[-inTrain,]%>%
  drop_na()
rf <- randomForest(log_AE ~ ., data = training, importance=TRUE, keep.inbag = TRUE)
rf
varImpPlot(rf)
#%IncMSE is the increase in mean square error when that variable is being randomly shuffled. A larger value means that the variable is more important for the model. %IncNodePurity is slightly more complicated but relates to the distribution of values in the split, suffice to say that more useful variables tend to have higher purity. Reassuringly enough the values seem to tell pretty much the same story.
#checked varImpPlot, and then tried without percentage_gene_gc_content, transcript_length, transitivity ==> correlation increased 

train.pred <- predict(rf, newdata = training)
cor(training$log_AE, train.pred)
plot(training$log_AE, train.pred)

valid.pred <- predict(rf, newdata = validation)
cor(validation$log_AE, valid.pred)
plot(validation$log_AE, valid.pred)


test.pred <- predict(rf, newdata = test)
cor(test$log_AE, test.pred, method='pearson')
corr <- round(cor(test$log_AE, test.pred, method='pearson'),3)
plot(test$log_AE, test.pred)
df2 <- data.frame(test$log_AE, test.pred)

library('forestError')
mspes <- quantForestError(rf, training, test, what = 'mspe')
rownames(mspes) <-  rownames(test)
mspes$real_log_AE <- test$log_AE
mean(mspes$mspe)
rmse <- sqrt(mean((df2$test.pred - df2$test.log_AE)^2))
rmse/mean(df2$test.log_AE)
ggp<- ggplot(df2, aes(test.log_AE, test.pred))+
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_text( x = (-7), y= (-3.5),label = paste('Pearson correlation: ', corr))

ggp

#save(file='Rdata/supervised.Rdata', list=ls())
