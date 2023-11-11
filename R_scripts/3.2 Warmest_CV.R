library(randomForest)
library(mgcv)
library(tidyverse)
library(caret)
library(pROC)
library(here)

#################################################################################################################
######################## INDEPENDENT VALIDATION, 90% AND HIGHER GDD VALUES BLOCKED ############################
#################################################################################################################

Train.fileNames = list.files(path=here("Data","TrainData"),pattern=".csv")
Train.fileNames

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  q=quantile(sub.df[,4], probs=0.9) ### set the threshold at 90th percentile
  test.df=sub.df%>%filter(.[[4]]>q) ### test data containing lakes with upper 10th percent temperatures
  train.df=sub.df%>%filter(.[[4]]<q) ### train data all but
  
  rf = randomForest(as.factor(EWMSTATUS)~., train.df[,-1 ], ntree = 500, data=train.df, keep.forest=TRUE)
  test.df$preds=NULL
  test.df$preds=predict(rf, test.df[,-1], type = "prob")[,2]
  
  auc_obj=evalmod(scores=test.df$preds,labels=test.df$EWMSTATUS, mode="aucroc")
  AUC=auc_obj$uaucs$aucs
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, AUC)) 
  
  preds.val=ifelse(test.df$preds > mean(test.df$preds),1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS))$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS))$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS))$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, spec))
  
}

#AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
#AUC_all
#AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)

AUC_RF_all=AUC_all
Acc_RF_all=Acc_all
Sens_RF_all=Sens_all
Spec_RF_all=Spec_all

AUC_RF_all
Acc_RF_all
Sens_RF_all
Spec_RF_all

#####################################################################################################
#####################################################################################################
######### Now for GAMs with min K and optimized K
######### Let's begin with s3.r.t3

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  q=quantile(sub.df[,4], probs=0.9)
  test.df=sub.df%>%filter(.[[4]]>q)
  train.df=sub.df%>%filter(.[[4]]<q)
  
  fm <- paste('s(', names(train.df[ 2 ]),' ,k=3)', ' + ', names(train.df[ 3 ]), 
              ' + ', 's(', names(train.df[ 4 ]), ',k= 3)')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  GAM_s3.r.t3 = gam(fm,data=train.df, method="REML", family = "binomial")
  test.df$GAMpreds=NULL
  test.df$GAMpreds <- predict(GAM_s3.r.t3, test.df[, -1], type = "response") # predict the test set
  
  auc_obj=evalmod(scores=test.df$GAMpreds,labels=test.df$EWMSTATUS, mode = "aucroc")
  AUC=auc_obj$uaucs$aucs
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, AUC)) 
  
  preds.val=ifelse(test.df$GAMpreds >mean(test.df$GAMpreds),1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, spec))
  
}



AUC_GAM.minK_all=AUC_all
Acc_GAM.minK_all=Acc_all
Sens_GAM.minK_all=Sens_all
Spec_GAM.minK_all=Spec_all

AUC_GAM.minK_all
Acc_GAM.minK_all
Sens_GAM.minK_all
Spec_GAM.minK_all

######### Next with s3.r.t7; the best K
AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  q=quantile(sub.df[,4], probs=0.9)
  test.df=sub.df%>%filter(.[[4]]>q)
  train.df=sub.df%>%filter(.[[4]]<q)
  
  fm <- paste('s(', names(train.df[ 2 ]),' ,k=3)', ' + ', names(train.df[ 3 ]), 
              ' + ', 's(', names(train.df[ 4 ]), ',k= 7)')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  GAM_s3.r.t7 = gam(fm,data=train.df, method="REML", family = "binomial")
  test.df$GAMpreds=NULL
  test.df$GAMpreds <- predict(GAM_s3.r.t7, test.df[, -1], type = "response") # predict the test set
  
  auc_obj=evalmod(scores=test.df$GAMpreds,labels=test.df$EWMSTATUS, mode = "aucroc")
  AUC=auc_obj$uaucs$aucs
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, AUC)) 
  
  preds.val=ifelse(test.df$GAMpreds >mean(test.df$GAMpreds),1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, spec))
  
}



AUC_GAM.bestK_all=AUC_all
Acc_GAM.bestK_all=Acc_all
Sens_GAM.bestK_all=Sens_all
Spec_GAM.bestK_all=Spec_all

AUC_GAM.bestK_all
Acc_GAM.bestK_all
Sens_GAM.bestK_all
Spec_GAM.bestK_all



GAM.minK_Warmest_perf=left_join(AUC_GAM.minK_all,Sens_GAM.minK_all,by="Train.fileName")%>%
  left_join(Spec_GAM.minK_all, by="Train.fileName")%>%left_join(Acc_GAM.minK_all, by="Train.fileName")
GAM.bestK_Warmest_perf=left_join(AUC_GAM.bestK_all,Sens_GAM.bestK_all,by="Train.fileName")%>%
  left_join(Spec_GAM.bestK_all, by="Train.fileName")%>%left_join(Acc_GAM.bestK_all, by="Train.fileName")
RF_Warmest_perf=left_join(AUC_RF_all,Sens_RF_all,by="Train.fileName")%>%
  left_join(Spec_RF_all, by="Train.fileName")%>%left_join(Acc_RF_all, by="Train.fileName")

AllSDMs_Warmest_CV=bind_rows(GAM.minK_Warmest_perf, GAM.bestK_Warmest_perf,RF_Warmest_perf )%>%
  mutate(SDM=c(rep("GAM.minK",5),rep("GAM.bestK",5),rep("RF",5)))
AllSDMs_Warmest_CV

write.table(AllSDMs_Warmest_CV,"Results/AllSDMs_Warmest_CV.txt", sep="\t")

######################################################################################################
######################################################################################################
##### Putting it all together into a table
here()
AUC.fileNames = list.files(pattern="CV.txt", recursive = TRUE)
AUC.fileNames
result = lapply(AUC.fileNames, function(x) read.table(x, header = TRUE))
result
names(result[[3]])[2:5]=c("meanAUC", "meanSens", "meanSpec","meanAcc")

AllAUCs_combined=do.call(rbind, result)
AllAUCs_combined

Val.Results=AllAUCs_combined%>%mutate(SDM=recode(SDM, GAM.bestK = "GAM (k=7)", GAM.minK = "GAM (k=3)"))%>%
  mutate(Validation = c(rep("Random 5-fold", 15),rep("Spatial-block", 15),rep("Warmest 10%", 15)))%>%
  group_by(Validation,SDM)%>%summarise(
    avgAUC=round(mean(meanAUC),2),
    avgAcc=round(mean(meanAcc),2),
    avgSens=round(mean(meanSens),2),
    avgSpec=round(mean(meanSpec),2))

library(flextable)
flextable(Val.Results) %>%
  theme_vanilla()%>%save_as_pptx(path="Results/Val.Results_wAcc.table.pptx")

######################################################################################################
######################################################################################################

























######################################################################################################
############################# Previous GAM models that are NO more used ##########################
######################################################################################################
################## NOW WITH GAM, K=10
AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  q=quantile(sub.df[,4], probs=0.9)
  test.df=sub.df%>%filter(.[[4]]>q)
  train.df=sub.df%>%filter(.[[4]]<q)
  
  fm <- paste('s(', names(train.df[ -1 ]), ',k=10)', sep = "", collapse = ' + ')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  GAM_k10 = gam(fm,data=train.df, method="REML", family = "binomial")
  test.df$GAMpreds=NULL
  test.df$GAMpreds <- predict(GAM_k10, test.df[, -1], type = "response") # predict the test set
  
  auc_obj=evalmod(scores=test.df$GAMpreds,labels=test.df$EWMSTATUS, mode = "aucroc")
  AUC=auc_obj$uaucs$aucs
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, AUC)) 
  
  preds.val=ifelse(test.df$GAMpreds >mean(test.df$GAMpreds),1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, spec))
  
}

AUC_GAM.k10_all=AUC_all
Acc_GAM.k10_all=Acc_all
Sens_GAM.k10_all=Sens_all
Spec_GAM.k10_all=Spec_all

AUC_GAM.k10_all
Acc_GAM.k10_all
Sens_GAM.k10_all
Spec_GAM.k10_all

################## NOW WITH GAM, K=3
AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  q=quantile(sub.df[,4], probs=0.9)
  test.df=sub.df%>%filter(.[[4]]>q)
  train.df=sub.df%>%filter(.[[4]]<q)
  
  fm <- paste('s(', names(train.df[ -1 ]), ',k=3)', sep = "", collapse = ' + ')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  GAM_k3 = gam(fm,data=train.df, method="REML", family = "binomial")
  test.df$GAMpreds=NULL
  test.df$GAMpreds <- predict(GAM_k3, test.df[, -1], type = "response") # predict the test set
  
  auc_obj=evalmod(scores=test.df$GAMpreds,labels=test.df$EWMSTATUS, mode = "aucroc")
  AUC=auc_obj$uaucs$aucs
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, AUC)) 
  
  preds.val=ifelse(test.df$GAMpreds >mean(test.df$GAMpreds),1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, spec))
  
}



AUC_GAM.k3_all=AUC_all
Acc_GAM.k3_all=Acc_all
Sens_GAM.k3_all=Sens_all
Spec_GAM.k3_all=Spec_all

AUC_GAM.k3_all
Acc_GAM.k3_all
Sens_GAM.k3_all
Spec_GAM.k3_all
