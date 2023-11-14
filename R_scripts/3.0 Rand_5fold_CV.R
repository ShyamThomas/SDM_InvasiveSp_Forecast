library(randomForest)
library(mgcv)
library(tidyverse)
library(caret)
library(pROC)
library(here)

here()

#### Step1. Start by creating a list of all the saved files in the new folder
Train.fileNames = list.files(path=here("Data","TrainData"),pattern=".csv")
Train.fileNames

##### Step 2. Execute 5-fold cross-validation and capture AUCs for each of the 5 RF models

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL


for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  folds = rep_len(1:5,nrow(full.df))
  sample.folds=sample(folds,nrow(full.df))
  full.df$folds=sample.folds
  
  set.seed(007)
  
  for(i in 1:5){test.data=full.df[full.df$folds==i,]
  train.data= full.df[full.df$folds !=i,]
  train.rf = randomForest(train.data[,c(2:4)], as.factor(train.data$EWMSTATUS),
                          importance=TRUE, ntree=5000, type="regression")
  preds.test=predict(train.rf, newdata=test.data, type="prob")
  precrec_obj <- evalmod(scores =  preds.test[,2], labels =test.data$EWMSTATUS,
                         mode="aucroc")
  AUC=precrec_obj$uaucs$aucs
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  
  preds.val=ifelse(preds.test[,2] >0.5,1,0)
  acc=confusionMatrix(as.factor(preds.val), as.factor(test.data$EWMSTATUS))$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, k, acc))
  
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, k, sens))
  
  spec=confusionMatrix(as.factor(preds.val), as.factor(test.data$EWMSTATUS))$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, k, spec))
  
  }
}

### Summarize by GCM, and get rid off all the unwanted letters

MeanAUC_GCMs_RF=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))

#write.table(MeanAUC_GCMs_RF,"Results/AllGCMs_5foldCV_RF_AUCs.txt", sep="\t")

MeanAcc_GCMs_RF=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))

MeanSens_GCMs_RF=Sens_all%>%group_by(Train.fileName)%>%summarise(
meanSens=mean(sens))

#write.table(MeanSens_GCMs_RF,"Results/AllGCMs_5foldCV_RF_Sens.txt", sep="\t")

MeanSpec_GCMs_RF=Spec_all%>%group_by(Train.fileName)%>%summarise(
meanSpec=mean(spec))

MeanAUC_GCMs_RF
MeanAcc_GCMs_RF
MeanSens_GCMs_RF
MeanSpec_GCMs_RF
##############################################################################
### New GAMs; min K and tuned best K: Begins Here!
### GAM s3,r,t3
AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  folds = rep_len(1:5,nrow(sub.df))
  sample.folds=sample(folds,nrow(full.df))
  sub.df$folds=sample.folds
  
  set.seed(007)
  
  for(i in 1:5){test.data=sub.df[sub.df$folds==i,]
  train.data= sub.df[sub.df$folds !=i,]
  fm <- paste('s(', names(sub.df[ 2 ]),' ,k=3)', ' + ',  names(sub.df[ 3 ]), 
              ' + ', 's(', names(sub.df[ 4 ]), ',k= 3)')  ### FOR s=3, t=3
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_s3.r.t3 = gam(fm,data=train.data, method="REML", family = "binomial")
  preds.test=predict(gam_s3.r.t3, newdata=test.data, type="response")
  
  AUC=auc(roc(test.data$EWMSTATUS,preds.test))
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  
  preds.val=ifelse(preds.test >0.5,1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, i, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, i, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, i, spec))
  
  
  }
}

MeanAUC_GCMs_GAM.s3.r.t3=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAcc_GCMs_GAM.s3.r.t3=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))
MeanSens_GCMs_GAM.s3.r.t3=Sens_all%>%group_by(Train.fileName)%>%summarise(
  meanSens=mean(sens))
MeanSpec_GCMs_GAM.s3.r.t3=Spec_all%>%group_by(Train.fileName)%>%summarise(
  meanSpec=mean(spec))

MeanAUC_GCMs_GAM.s3.r.t3
MeanAcc_GCMs_GAM.s3.r.t3
MeanSens_GCMs_GAM.s3.r.t3
MeanSpec_GCMs_GAM.s3.r.t3

###################################### GAM s3.r.t7
AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  folds = rep_len(1:5,nrow(sub.df))
  sample.folds=sample(folds,nrow(full.df))
  sub.df$folds=sample.folds
  
  set.seed(007)
  
  for(i in 1:5){test.data=sub.df[sub.df$folds==i,]
  train.data= sub.df[sub.df$folds !=i,]
  fm <- paste('s(', names(sub.df[ 2 ]),' ,k=3)', ' + ',  names(sub.df[ 3 ]), 
              ' + ', 's(', names(sub.df[ 4 ]), ',k= 7)')  ### FOR s=3, t=7
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_s3.r.t7 = gam(fm,data=train.data, method="REML", family = "binomial")
  preds.test=predict(gam_s3.r.t7, newdata=test.data, type="response")
  
  AUC=auc(roc(test.data$EWMSTATUS,preds.test))
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  
  preds.val=ifelse(preds.test >0.5,1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, i, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, i, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, i, spec))
  }
}

MeanAUC_GCMs_GAM.s3.r.t7=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAcc_GCMs_GAM.s3.r.t7=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))
MeanSens_GCMs_GAM.s3.r.t7=Sens_all%>%group_by(Train.fileName)%>%summarise(
  meanSens=mean(sens))
MeanSpec_GCMs_GAM.s3.r.t7=Spec_all%>%group_by(Train.fileName)%>%summarise(
  meanSpec=mean(spec))

MeanAUC_GCMs_GAM.s3.r.t7
MeanAcc_GCMs_GAM.s3.r.t7
MeanSens_GCMs_GAM.s3.r.t7
MeanSpec_GCMs_GAM.s3.r.t7

MeanAUC_GCMs_GAM.s3.r.t3
MeanAcc_GCMs_GAM.s3.r.t3
MeanSens_GCMs_GAM.s3.r.t3
MeanSpec_GCMs_GAM.s3.r.t3

GAM.minK_perf=left_join(MeanAUC_GCMs_GAM.s3.r.t3, MeanSens_GCMs_GAM.s3.r.t3,by="Train.fileName")%>%
  left_join(MeanSpec_GCMs_GAM.s3.r.t3, by="Train.fileName")%>%
  left_join(MeanAcc_GCMs_GAM.s3.r.t3, by="Train.fileName")
GAM.bestK_perf=left_join(MeanAUC_GCMs_GAM.s3.r.t7, MeanSens_GCMs_GAM.s3.r.t7,by="Train.fileName")%>%
  left_join(MeanSpec_GCMs_GAM.s3.r.t7, by="Train.fileName")%>%
  left_join(MeanAcc_GCMs_GAM.s3.r.t7, by="Train.fileName")
RF_perf=left_join(MeanAUC_GCMs_RF,MeanSens_GCMs_RF,by="Train.fileName")%>%
  left_join(MeanSpec_GCMs_RF, by="Train.fileName")%>%left_join(MeanAcc_GCMs_RF, by="Train.fileName")
AllSDMs_5foldCV=bind_rows(GAM.minK_perf, GAM.bestK_perf,RF_perf )%>%
  mutate(SDM=c(rep("GAM.minK",5),rep("GAM.bestK",5),rep("RF",5)))
AllSDMs_5foldCV

write.table(AllSDMs_5foldCV,"Results/AllSDMs_5foldCV.txt", sep="\t")
###############################################################################################################
