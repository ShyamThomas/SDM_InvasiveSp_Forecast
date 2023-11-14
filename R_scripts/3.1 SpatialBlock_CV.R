
library(precrec)
library(maptools)
library(spatstat)
library(blockCV)
library(sf)
library(tidyverse)
library(here)

here()
#### Plot distribution of all EWM data across Minnesota
Minn.sf=read_sf(dsn=here("Data","GIS_Data"), layer="Minn.map")
plot(Minn.sf$geometry)
Minn.sf
st_crs(Minn.sf)

EWM.GCMs.data=read_csv("Data/EWM.prsabs95to15_AllGCMs_v2.csv")
EWM.GCMs.data

EWM.GCMs.sf=st_as_sf(EWM.GCMs.data, coords=c("LON","LAT"), crs=32615)
EWM.GCMs.sf
st_crs(EWM.GCMs.sf)

Minn.sf=st_transform(Minn.sf, crs=32615)
ggplot(Minn.sf)+geom_sf()+
  geom_sf(data=EWM.GCMs.sf)


################ SPATIALLY BLOCKED CV ALONG LATITUDINAL GRADIENT
CheqBrd.trial.blocks = cv_spatial(x = EWM.GCMs.sf, # sf or SpatialPoints
                                 column = "EWMSTATUS", # the response column (binomial or multi-class)
                                 #theRange = 100000, # size of the blocks in meters
                                 hexagon = FALSE,
                                 rows_cols = c(4,1), # number of folds
                                 selection = "checkerboard",
                                 iteration = 100, # find evenly dispersed folds
                                 biomod2 = FALSE)

CB=ggplot(CheqBrd.trial.blocks$blocks)+
            geom_sf(data = Minn.sf, fill=NA, colour="gray", lwd=1)+
            geom_sf(data = EWM.GCMs.sf, alpha = 0.25, aes(colour=as.factor(EWMSTATUS)))+theme(legend.title=element_blank())
CB+theme_minimal()+theme(legend.title=element_blank())+ggtitle("Spatial blocks", subtitle = "EWM distribution")

p=ggplot()+geom_sf(data = Minn.sf, fill=NA, colour="gray", lwd=1)+geom_sf(data = EWM.GCMs.sf, alpha = 0.5, aes(colour=as.factor(EWMSTATUS)))+theme(legend.title=element_blank())
p+theme_minimal()+theme(legend.title=element_blank())+ggtitle("EWM distribution")


############## START WITH RANDOM FOREST MODELS
Train.fileNames = list.files(here("Data","TrainData"),pattern=".csv")
Train.fileNames


AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL


for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  full.df$EWMSTATUS=as.factor(full.df$EWMSTATUS)

  folds=CheqBrd.trial.blocks$folds_list

  testTable <- full.df

  testTable$RFpreds <- NA

  for(k in seq_len(length(folds))){
    # extracting the training and testing indices
    # this way works with folds list (but not foldID)
    trainSet <- unlist(folds[[k]][1]) # training set indices
    testSet <- unlist(folds[[k]][2]) # testing set indices
    rf <- randomForest(EWMSTATUS~.,full.df[trainSet,], ntree = 500,importance=TRUE, type="regression") # model fitting on training set
    testTable$RFpreds[testSet] <- predict(rf, newdata = full.df[testSet, -1], type="prob")[,2] # predict the test set
   
    precrec_obj <- evalmod(scores = testTable$RFpreds[testSet], labels = testTable$EWMSTATUS[testSet],
                           mode="aucroc")
    AUC=precrec_obj$uaucs$aucs
    AUC_all = rbind(AUC_all, data.frame(Train.fileName, k, AUC))
    
    preds.val=ifelse(testTable$RFpreds[testSet] >0.5,1,0)
    acc=confusionMatrix(as.factor(preds.val),testTable$EWMSTATUS[testSet])$overall[1]
    Acc_all = rbind(Acc_all, data.frame(Train.fileName, k, acc))
    sens=confusionMatrix(as.factor(preds.val),testTable$EWMSTATUS[testSet])$byClass[1]
    Sens_all = rbind(Sens_all, data.frame(Train.fileName, k, sens))
    spec=confusionMatrix(as.factor(preds.val),testTable$EWMSTATUS[testSet])$byClass[2]
    Spec_all = rbind(Spec_all, data.frame(Train.fileName, k, spec))

  }
}


MeanAUC_SpBlk_RF=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAcc_SpBlk_RF=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))
MeanSens_SpBlk_RF=Sens_all%>%group_by(Train.fileName)%>%summarise(
  meanSens=mean(sens))
MeanSpec_SpBlk_RF=Spec_all%>%group_by(Train.fileName)%>%summarise(
  meanSpec=mean(spec))

MeanAUC_SpBlk_RF
MeanAcc_SpBlk_RF
MeanSens_SpBlk_RF
MeanSpec_SpBlk_RF

#write.table(MeanAUC_SpBlk_RF,"Results/AllGCMs_SpatialBlockCV_RF_AUCs.txt", sep="\t")

#### Optional codes to trim file names
#AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
#AUC_all
#AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
#AUC_all

##################################################################################################
######################NEW GAMs with min K and fine-tuned K
##################### min K BEGINS: s=3, r=1, t=3

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  #sub.df=full.df
  
  folds=CheqBrd.trial.blocks$folds_list
  testTable <- full.df
  testTable$GAMpreds <- NA
  
  for(k in seq_len(length(folds))){
    trainSet <- unlist(folds[[k]][1]) # training set indices
    testSet <- unlist(folds[[k]][2]) # testing set indices
    sample_sub=full.df[trainSet, ]
    fm <- paste('s(', names(sample_sub[ 2 ]),' ,k=3)', ' + ', names(sample_sub[ 3 ]), 
                ' + ', 's(', names(sample_sub[ 4 ]), ',k= 3)')
    fm <- as.formula(paste('EWMSTATUS ~', fm))
    GAM_s3.r.t3 = gam(fm,data=sample_sub, method="REML", family = "binomial")
    testTable$GAMpreds[testSet] <- predict(GAM_s3.r.t3, full.df[testSet, ], type = "response") # predict the test set
    
    precrec_obj <- evalmod(scores = testTable$GAMpreds[testSet], labels = testTable$EWMSTATUS[testSet],mode="aucroc")
    AUC=precrec_obj$uaucs$aucs
    AUC_all = rbind(AUC_all, data.frame(Train.fileName, k, AUC))
    
    preds.val=ifelse(testTable$GAMpreds[testSet] >0.5,1,0)
    acc=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$overall[1]
    Acc_all = rbind(Acc_all, data.frame(Train.fileName, k, acc))
    sens=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$byClass[1]
    Sens_all = rbind(Sens_all, data.frame(Train.fileName, k, sens))
    spec=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$byClass[2]
    Spec_all = rbind(Spec_all, data.frame(Train.fileName, k, spec))
  }
}


MeanAUC_SpBlk_GAM.minK=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAcc_SpBlk_GAM.minK=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))
MeanSens_SpBlk_GAM.minK=Sens_all%>%group_by(Train.fileName)%>%summarise(
  meanSens=mean(sens))
MeanSpec_SpBlk_GAM.minK=Spec_all%>%group_by(Train.fileName)%>%summarise(
  meanSpec=mean(spec))

MeanAUC_SpBlk_GAM.minK
MeanAcc_SpBlk_GAM.minK
MeanSens_SpBlk_GAM.minK
MeanSpec_SpBlk_GAM.minK

##################### Best K BEGINS: s=3, r=1, t=7
AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("Data/TrainData/",Train.fileName, sep=""))
  #sub.df=full.df
  
  folds=CheqBrd.trial.blocks$folds_list
  testTable <- full.df
  testTable$GAMpreds <- NA
  
  for(k in seq_len(length(folds))){
    trainSet <- unlist(folds[[k]][1]) # training set indices
    testSet <- unlist(folds[[k]][2]) # testing set indices
    sample_sub=full.df[trainSet, ]
    fm <- paste('s(', names(sample_sub[ 2 ]),' ,k=3)', ' + ', names(sample_sub[ 3 ]), 
                ' + ', 's(', names(sample_sub[ 4 ]), ',k= 7)')
    fm <- as.formula(paste('EWMSTATUS ~', fm))
    GAM_s3.r.t7 = gam(fm,data=sample_sub, method="REML", family = "binomial")
    testTable$GAMpreds[testSet] <- predict(GAM_s3.r.t7, full.df[testSet, ], type = "response") # predict the test set
    
    precrec_obj <- evalmod(scores = testTable$GAMpreds[testSet], labels = testTable$EWMSTATUS[testSet],mode="aucroc")
    AUC=precrec_obj$uaucs$aucs
    AUC_all = rbind(AUC_all, data.frame(Train.fileName, k, AUC))
    
    preds.val=ifelse(testTable$GAMpreds[testSet] >0.5,1,0)
    acc=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$overall[1]
    Acc_all = rbind(Acc_all, data.frame(Train.fileName, k, acc))
    sens=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$byClass[1]
    Sens_all = rbind(Sens_all, data.frame(Train.fileName, k, sens))
    spec=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$byClass[2]
    Spec_all = rbind(Spec_all, data.frame(Train.fileName, k, spec))
  }
}


MeanAUC_SpBlk_GAM.bestK=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAcc_SpBlk_GAM.bestK=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))
MeanSens_SpBlk_GAM.bestK=Sens_all%>%group_by(Train.fileName)%>%summarise(
  meanSens=mean(sens))
MeanSpec_SpBlk_GAM.bestK=Spec_all%>%group_by(Train.fileName)%>%summarise(
  meanSpec=mean(spec))

MeanAUC_SpBlk_GAM.bestK
MeanAcc_SpBlk_GAM.bestK
MeanSens_SpBlk_GAM.bestK
MeanSpec_SpBlk_GAM.bestK

MeanAUC_SpBlk_GAM.minK
MeanAcc_SpBlk_GAM.minK
MeanSens_SpBlk_GAM.minK
MeanSpec_SpBlk_GAM.minK

GAM.minK_SpBlk_perf=left_join(MeanAUC_SpBlk_GAM.minK,MeanSens_SpBlk_GAM.minK,by="Train.fileName")%>%
  left_join(MeanSpec_SpBlk_GAM.minK, by="Train.fileName")%>%
  left_join(MeanAcc_SpBlk_GAM.minK, by="Train.fileName")

GAM.bestK_SpBlk_perf=left_join(MeanAUC_SpBlk_GAM.bestK,MeanSens_SpBlk_GAM.bestK,by="Train.fileName")%>%
  left_join(MeanSpec_SpBlk_GAM.bestK, by="Train.fileName")%>%
  left_join(MeanAcc_SpBlk_GAM.bestK, by="Train.fileName")

RF_SpBlk_perf=left_join(MeanAUC_SpBlk_RF, MeanSens_SpBlk_RF,by="Train.fileName")%>%
  left_join(MeanSpec_SpBlk_RF, by="Train.fileName")%>%
  left_join(MeanAcc_SpBlk_RF, by="Train.fileName")

AllSDMs_SpBlk_CV=bind_rows(GAM.minK_SpBlk_perf, GAM.bestK_SpBlk_perf,RF_SpBlk_perf )%>%
  mutate(SDM=c(rep("GAM.minK",5),rep("GAM.bestK",5),rep("RF",5)))

write.table(AllSDMs_SpBlk_CV,"Results/AllSDMs_SpBlk_CV.txt", sep="\t")

######################################################################################################################
######################################################################################################################
