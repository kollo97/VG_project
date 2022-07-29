library('e1071')
library('caret')
library('dplyr')
library('tidyverse')
##TESTING IF SUPPORT VECTOR MACHINE WORKS BETTER THAN RF (I.E. WITH THE KNOWLEDGE OF SVM THAT I HAVE... (==> NOT SO MUCH))

remove(list=ls())
load(file='Rdata/visualization.Rdata')
superv_df <- df %>%
  select(ensembl_gene_id,Avg_VG_AE, loeuf_score, ncGERP, ncRVIS,  RVIS_score,  Avg_VG_eQTL,  median_tpm)%>%
  drop_na()%>%
  column_to_rownames(var='ensembl_gene_id')

#possible variables: transcript_length, percentage_gene_gc_content, loeuf_score, ncGERP, ncRVIS, transitivity, degree, RVIS_score  
#,percentage_gene_gc_content, transcript_length, transitivity
##ADD LOG AE_VG and eQTL_VG
superv_df$log_AE <- log(superv_df$Avg_VG_AE)
superv_df$log_eQTL <- log(superv_df$Avg_VG_eQTL)
subset <- superv_df %>%
  select(!c(contains('VG')))%>% 
  drop_na()#%>%
lapply(subset, class)
subset[,names(subset)] <- lapply(subset[,names(subset)], as.numeric)

#map(superv_df, ~var(.))
set.seed(123)
inTrain <- createDataPartition(y = subset$log_AE, p=0.80, list =F)
training <- subset[inTrain,]%>%
  drop_na()
inValidation <- createDataPartition(y = training$log_AE, p=0.80, list =F)
training <- training[inValidation,]%>%
  drop_na()
validation <- training[-inValidation,]%>%
  drop_na()

test <- subset[-inTrain,]%>%
  drop_na()



svm <- best.svm(x = training[,-6] , y=training$log_AE, gamma = c(0.001, 0.01, 0.1), cost = c(4,8,16), probability = T, tunecontrol = tune.control(cross=3),kernel='linear')
#default_svm <- svm(x = training[,-6] , y=training$log_AE, gamma = 0.1, cost = 8, probability = T, cross=3,kernel='linear') 
train.pred <- predict(svm, newdata = training[,-6])
plot(training$log_AE, train.pred)
cor(training$log_AE, train.pred,  method= 'spearman')

valid.pred <- predict(svm, newdata = validation[,-6])
plot(validation$log_AE, valid.pred)
cor(validation$log_AE, valid.pred,  method= 'spearman')


test.pred <- predict(svm, newdata = test[,-6])
plot(test$log_AE, test.pred)
cor(test$log_AE, test.pred)
