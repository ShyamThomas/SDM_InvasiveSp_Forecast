library(sf)
library(randomForest)
library(mgcv)
library(tidyverse)
library(corrplot)
library(pdp)
library(miceadds)

####################################################################################################################
######### R script: Training EWM data with Random forest (RF) & GAM SDMs, and evaluate model performance   #########
########  with random 5-fold cross-validation.                                                             #########
####################################################################################################################

#### Step 1. Load the merged water temperature and EWM data developed in the previous modules and subset all
#### 5 GCM based temperature estimates.
EWM.GCMs.data=read_csv("Data/EWM.prsabs95to15_AllGCMs_v2.csv")
EWM.GCMs.data

-#### Step 1a. GCM-ACCESS estimates
EWM.train.data_ACCESS.WtrTemp=EWM.GCMs.data[,-c(1:6,9,12:15)]
EWM.train.data_ACCESS.WtrTemp
write_csv(EWM.train.data_ACCESS.WtrTemp, "Data/TrainData/EWM.train.data_ACCESS.WtrTemp.csv")

#### Step 1b. GCM-MIROC5 estimates
EWM.train.data_MIROC5.WtrTemp=EWM.GCMs.data[,-c(1:6,9,11,13:15)]
EWM.train.data_MIROC5.WtrTemp
write_csv(EWM.train.data_MIROC5.WtrTemp, "Data/TrainData/EWM.train.data_MIROC5.WtrTemp.csv")

#### Step 1c. GCM-IPSL estimates
EWM.train.data_IPSL.WtrTemp=EWM.GCMs.data[,-c(1:6,9,11:12,14:15)]
EWM.train.data_IPSL.WtrTemp
write_csv(EWM.train.data_IPSL.WtrTemp, "Data/TrainData/EWM.train.data_IPSL.WtrTemp.csv")

#### Step 1d. GCM-GFDL estimates
EWM.train.data_GFDL.WtrTemp=EWM.GCMs.data[,-c(1:6,9,11:13,15)]
EWM.train.data_GFDL.WtrTemp
write_csv(EWM.train.data_GFDL.WtrTemp, "Data/TrainData/EWM.train.data_GFDL.WtrTemp.csv")

#### Step 1e. GCM-MRI estimates
EWM.train.data_MRI.WtrTemp=EWM.GCMs.data[,-c(1:6,9,11:14)]
EWM.train.data_MRI.WtrTemp
write_csv(EWM.train.data_MRI.WtrTemp, "Data/TrainData/EWM.train.data_MRI.WtrTemp.csv")

#### Step 2. Iterating RF models across all the above created data sets for 5 GCMs:
#### start by creating a list of all the saved files in the new folder
Train.fileNames = list.files(path="Data/TrainData/",pattern=".csv")
Train.fileNames

########################################################################
### A loop that reads all the files and executes RANDOM FOREST algorithm; saves RF obj and OOB errors as text file,
### and finally, builds and saves partial plots for each variable within each model

results_top3=NULL

for(Train.fileName in Train.fileNames) {
  sample = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  rf = randomForest(sample[,(2:4)], sample$EWMSTATUS,importance=TRUE, ntree=5000, type="regression")
  save(rf, file=paste("Data/TrainData/", sub('....$','',Train.fileName), "RF.Rdata", sep=""))
  
  meanMSE = mean(rf$mse)
  results_top3 = rbind(results_top3, data.frame(Train.fileName, meanMSE))
  write.table(results_top3,"Results/RF_MSE_Top3.txt",sep = "\t")
  
  Top3Preds.Names=colnames(sample)[c(2:4)]
  
  for (Pred.Name in Top3Preds.Names){
    partial_plot=autoplot(partial(rf, pred.var = Pred.Name, ice=TRUE, rug=TRUE, train = sample, prob = TRUE),xlab=Pred.Name, ylab="Invasion risk", alpha=0.1)
    ggsave(filename=paste(sub('....$','',Train.fileName),"_",Pred.Name,"_IcePlot.png", sep=""), partial_plot, path="Figures/",units="in", width=9, height=6, dpi=900)
    }
}

##### Step 3. Execute 5-fold cross-validation and capture AUCs for each of the 5 RF models
library(pROC)

AUC_all=NULL


for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  folds = rep_len(1:5,nrow(full.df))
  sample.folds=sample(folds,nrow(sample))
  full.df$folds=sample.folds
  
  set.seed(007)

  for(i in 1:5){test.data=full.df[full.df$folds==i,]
                train.data= full.df[full.df$folds !=i,]
                train.rf = randomForest(train.data[,c(2:4)], train.data$EWMSTATUS,importance=TRUE, ntree=5000, 
                                        type="regression")
                preds.test=predict(train.rf, newdata=test.data)
                AUC=auc(roc(test.data$EWMSTATUS,preds.test))
                AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  }
}

### Get rid off all the unwanted letters
AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
AUC_all

MeanAUC_GCMs_RF=AUC_all%>%group_by(Train.fileName)%>%summarise(
 meanAUC=mean(AUC))
MeanAUC_GCMs_RF
#write.table(MeanAUC_GCMs_RF,"Results/AllGCMs_5foldCV_RF_AUCs.txt", sep="\t")


###############################################################################################
#### Optimizing k value for GAMs
### The minimum k value that gives a non-linearish response curves; k= 3 

for(Train.fileName in Train.fileNames) {
  sample = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  fm <- paste('s(', names(sample[ 2 ]),' ,k=3)', ' + ',  names(sample[ 3 ]), 
              ' + ', 's(', names(sample[ 4 ]), ',k= 3)')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_s3.r.t3 = gam(fm,data=sample, method="REML", family = "binomial")
  
  save(gam_s3.r.t3, file=paste("Data/TrainData/", sub('....$','',Train.fileName), "GAM_s3.r.t3.Rdata", sep=""))
}

load.Rdata("Data/TrainData/EWM.train.data_MIROC5.WtrTempGAM_s3.r.t3.Rdata", "MIROC5.GAM_minK.model")
gam.check(MIROC5.GAM_minK.model)

load.Rdata("Data/TrainData/EWM.train.data_MRI.WtrTempGAM_s3.r.t3.Rdata", "MRI.GAM_minK.model")
gam.check(MRI.GAM_minK.model)

load.Rdata("Data/TrainData/EWM.train.data_ACCESS.WtrTempGAM_s3.r.t3.Rdata", "ACCESS.GAM_minK.model")
gam.check(ACCESS.GAM_minK.model)

load.Rdata("Data/TrainData/EWM.train.data_GFDL.WtrTempGAM_s3.r.t3.Rdata", "GFDL.GAM_minK.model")
gam.check(GFDL.GAM_minK.model)

load.Rdata("Data/TrainData/EWM.train.data_IPSL.WtrTempGAM_s3.r.t3.Rdata", "IPSL.GAM_minK.model")
gam.check(IPSL.GAM_minK.model)

#### The tuned k values that where p-values of gam check is less than 0.05

for(Train.fileName in Train.fileNames) {
  sample = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  fm <- paste('s(', names(sample[ 2 ]),' ,k=3)', ' + ',  names(sample[ 3 ]), 
              ' + ', 's(', names(sample[ 4 ]), ',k= 7)')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_s3.r.t7 = gam(fm,data=sample, method="REML", family = "binomial")
  
  save(gam_s3.r.t7, file=paste("Data/TrainData/", sub('....$','',Train.fileName), "GAM_s3.r.t20.Rdata", sep=""))
}

load.Rdata("Data/TrainData/EWM.train.data_MIROC5.WtrTempGAM_s3.r.t7.Rdata", "MIROC5.GAM_bestK.model")
gam.check(MIROC5.GAM_bestK.model)
plot.gam(MIROC5.GAM_bestK.model)

load.Rdata("Data/TrainData/EWM.train.data_MRI.WtrTempGAM_s3.r.t7.Rdata", "MRI.GAM_bestK.model")
gam.check(MRI.GAM_bestK.model)
plot.gam(MRI.GAM_bestK.model)

load.Rdata("Data/TrainData/EWM.train.data_ACCESS.WtrTempGAM_s3.r.t7.Rdata", "ACCESS.GAM_bestK.model")
gam.check(ACCESS.GAM_bestK.model)
plot.gam(ACCESS.GAM_bestK.model)

load.Rdata("Data/TrainData/EWM.train.data_GFDL.WtrTempGAM_s3.r.t7.Rdata", "GFDL.GAM_bestK.model")
gam.check(GFDL.GAM_bestK.model)
plot.gam(GFDL.GAM_bestK.model)

load.Rdata("Data/TrainData/EWM.train.data_IPSL.WtrTempGAM_s3.r.t7.Rdata", "IPSL.GAM_bestK.model")
gam.check(IPSL.GAM_bestK.model)
plot.gam(IPSL.GAM_bestK.model)

##########################################################################################################
####### 5-fold cross validation based on the above min K and K-optimized GAMs
library(pROC)

#### First for minimum k values that gives non-linearish response
AUC_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  folds = rep_len(1:5,nrow(sub.df))
  sample.folds=sample(folds,nrow(sample))
  sub.df$folds=sample.folds
  
  set.seed(007)
  
  for(i in 1:5){test.data=sub.df[sub.df$folds==i,]
  train.data= sub.df[sub.df$folds !=i,]
  fm <- paste('s(', names(sub.df[ 2 ]),' ,k=3)', ' + ',  names(sub.df[ 3 ]), 
              ' + ', 's(', names(sub.df[ 4 ]), ',k= 3)')  ### FOR s=3, t=3
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_s3.r.t3 = gam(fm, data=train.data, method="REML", family = "binomial")
  preds.test=predict(gam_s3.r.t3, newdata=test.data, type="response")
  
  AUC=auc(roc(test.data$EWMSTATUS,preds.test))
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  }
}

### Get rid off all the unwanted letters
AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
AUC_all

MeanAUC_GCMs_s3.r.t3=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC)
)
MeanAUC_GCMs_s3.r.t3
#write.table(MeanAUC_GCMs_s3.r.t3,"Results/AllGCMs_5foldCV_GAM.minK_AUCs.txt", sep="\t")

#################################
#### First for best k values that gives best fit non-linear response
AUC_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  folds = rep_len(1:5,nrow(sub.df))
  sample.folds=sample(folds,nrow(sample))
  sub.df$folds=sample.folds
  
  set.seed(007)
  
  for(i in 1:5){test.data=sub.df[sub.df$folds==i,]
  train.data= sub.df[sub.df$folds !=i,]
  fm <- paste('s(', names(sub.df[ 2 ]),' ,k=3)', ' + ',  names(sub.df[ 3 ]), 
              ' + ', 's(', names(sub.df[ 4 ]), ',k= 7)')  ### FOR s=3, t=3
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_s3.r.t7 = gam(fm, data=train.data, method="REML", family = "binomial")
  preds.test=predict(gam_s3.r.t7, newdata=test.data, type="response")
  
  AUC=auc(roc(test.data$EWMSTATUS,preds.test))
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  }
}

### Get rid off all the unwanted letters
AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
AUC_all

MeanAUC_GCMs_s3.r.t7=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC)
)
MeanAUC_GCMs_s3.r.t7
#write.table(MeanAUC_GCMs_s3.r.t7,"Results/AllGCMs_5foldCV_GAM.BestK_AUCs.txt", sep="\t")

######################################################################################################
######################################################################################################


