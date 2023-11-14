####################################################################################################################
####### R script: Putting together all the SDM results with lake-specific covariates, and prediction domain #######
####### These scripts were used to put together the SINGLE LARGE DATASET needed by USpatial team @ UMN      #######
####################################################################################################################

library(tidyverse)
library(sf)
here()

### 1. Get response data and predictors: two sets of predictors for water GDD - present and future for 5 GCMs
####### 2 unchanging covariates, 10 GDD estimates (5GCMs*2 time-periods), and 5 domains (1*5GCMs) = 17 columns

EWM.CurrTrain.Data=read_csv(here("Data","EWM.prsabs95to15_AllGCMs_v2.csv"))
EWM.CurrTrain.Data
colnames(EWM.CurrTrain.Data)[11:15]=c("ACCESS.Curr", "MIROC5.Curr", "IPSL.Curr",   "GFDL.Curr" ,  "MRI.Curr")

EWM.FutrTest.Data=read_csv(here("Data", "EWM.prsabs40to60_AllGCMs_v2.csv"))
colnames(EWM.FutrTest.Data)[11:15]=c("ACCESS.Futr", "MIROC5.Futr", "IPSL.Futr",   "GFDL.Futr" ,  "MRI.Futr")
EWM.CurrFutr.Preds=left_join(EWM.CurrTrain.Data, EWM.FutrTest.Data[,c(1,11:15)], by="DOWLKNUM")

EWM.CurrFutr.Preds.Doms=EWM.CurrFutr.Preds%>%mutate(ACCESS.domain=case_when(ACCESS.Futr < max(ACCESS.Curr) ~ 'Analog', ACCESS.Futr > max(ACCESS.Curr)~ 'NonAnalog'),
                            MIROC5.domain=case_when(MIROC5.Futr < max(MIROC5.Curr) ~ 'Analog', MIROC5.Futr > max(MIROC5.Curr) ~ 'NonAnalog'),
                            IPSL.domain=case_when(IPSL.Futr < max(IPSL.Curr) ~ 'Analog', IPSL.Futr > max(IPSL.Curr) ~ 'NonAnalog'),
                            GFDL.domain=case_when(GFDL.Futr < max(GFDL.Curr) ~ 'Analog', GFDL.Futr > max(GFDL.Curr) ~ 'NonAnalog'),
                            MRI.domain=case_when(MRI.Futr < max(MRI.Curr) ~ 'Analog', MRI.Futr > max(MRI.Curr) ~ 'NonAnalog')
)
        
EWM.CurrFutr.Preds.Doms                    
write_csv(EWM.CurrFutr.Preds.Doms, "Data/EWM.CurrFutr.TempPreds.Doms.csv")

### 2. Predictions of EWM occurrence: two sets for each 5 GCM temp estimates and 3 SDM modeling algorithms
####### 30 columns: SDM_GCM_PERIOD
########### Start with GAM minK s3.r.t3
fut.preds.df=read_csv("Results/GAM.minK_Fut.Predictions.csv") ## from 2.0 GAM_Predictions output
curr.preds.df=read_csv("Results/GAM.minK_Curr.Predictions.csv")
currANDfut_preds=bind_rows(curr.preds.df,fut.preds.df)
head(currANDfut_preds) 
dim(currANDfut_preds)

## Split into 2 seperate datasets of 578 rows each
curr.preds.GAM_minK=currANDfut_preds[1:578,]
futr.preds.GAM_minK=currANDfut_preds[579:1156,]
### Give each column a unique name
colnames(futr.preds.GAM_minK)[1:6]=c("ACCESS.GAMminK.FutrPred" ,  "GFDL.GAMminK.FutrPred"   ,  
                                     "IPSL.GAMminK.FutrPred"  ,   "MIROC5.GAMminK.FutrPred"  ,
                                     "MRI.GAMminK.FutrPred"   ,   "DOWLKNUM")
colnames(curr.preds.GAM_minK)[1:6]=c("ACCESS.GAMminK.CurrPred" ,  "GFDL.GAMminK.CurrPred"   ,  
                                   "IPSL.GAMminK.CurrPred"  ,   "MIROC5.GAMminK.CurrPred"  , 
                                   "MRI.GAMminK.CurrPred"   ,   "DOWLKNUM")

########### Next with GAM best K, s3.r.t7
fut.preds.df=read_csv("Results/GAM.bestK_Fut.Predictions.csv") #from 2.1 GAM_Predictions output
curr.preds.df=read_csv("Results/GAM.bestK_Curr.Predictions.csv")
currANDfut_preds_bestK=bind_rows(curr.preds.df,fut.preds.df)
currANDfut_preds_bestK
dim(currANDfut_preds_bestK)

curr.preds.GAM_bestK=currANDfut_preds_bestK[1:578,]
futr.preds.GAM_bestK=currANDfut_preds_bestK[579:1156,]

colnames(futr.preds.GAM_bestK)[1:6]=c("ACCESS.GAMbestK.FutrPred" ,  "GFDL.GAMbestK.FutrPred"   , 
                                      "IPSL.GAMbestK.FutrPred"  ,   "MIROC5.GAMbestK.FutrPred"  , 
                                      "MRI.GAMbestK.FutrPred"   ,   "DOWLKNUM")
colnames(curr.preds.GAM_bestK)[1:6]=c("ACCESS.GAMbestK.CurrPred" ,  "GFDL.GAMbestK.CurrPred"   ,  
                                      "IPSL.GAMbestK.CurrPred"  ,   "MIROC5.GAMbestK.CurrPred"  , 
                                      "MRI.GAMbestK.CurrPred"   ,   "DOWLKNUM")

########### Finally with Random Forest model predictions
RF_curr.preds=read_csv("Results/RF_Curr.Predictions.csv")
RF_futr.preds=read_csv("Results/RF_Futr.Predictions.csv")
RF_curr.preds
RF_futr.preds
colnames(RF_curr.preds)[1:5]=c("ACCESS.RF.CurrPred" ,  "GFDL.RF.CurrPred"   ,  "IPSL.RF.CurrPred"  , 
                               "MIROC5.RF.CurrPred"  , "MRI.RF.CurrPred")
colnames(RF_futr.preds)[1:5]=c("ACCESS.RF.FutrPred" ,  "GFDL.RF.FutrPred"   ,  "IPSL.RF.FutrPred"  ,  
                               "MIROC5.RF.FutrPred"  , "MRI.RF.FutrPred")

RF_futr.preds$DOWLKNUM=curr.preds.GAM_minK$DOWLKNUM
RF_futr.preds
RF_curr.preds$DOWLKNUM=curr.preds.GAM_minK$DOWLKNUM
RF_curr.preds


Final_MergedData=left_join(EWM.CurrFutr.Preds.Doms,curr.preds.GAM_minK[,1:6], by="DOWLKNUM")%>%
                  left_join(.,futr.preds.GAM_minK[,1:6], by="DOWLKNUM")%>%
                        left_join(.,curr.preds.GAM_bestK[,1:6], by="DOWLKNUM")%>%
                        left_join(.,futr.preds.GAM_bestK[,1:6], by="DOWLKNUM")%>%
                             left_join(., RF_curr.preds[,-6],by="DOWLKNUM")%>%
                              left_join(., RF_futr.preds[,-6],by="DOWLKNUM")

Final_MergedData%>%View()

##### Now merge with Bayesian predictions of best GAM
Final_MergedData_V2=left_join(Final_MergedData,Curr.Brm.Preds_EWMGeoIndex[,c(5:7)], by="DOWLKNUM")%>%
                        left_join(.,Fut.Brm.Preds_EWMGeoIndex[,c(5,9,10)], by="DOWLKNUM")
Final_MergedData_V2
write_csv(Final_MergedData_V2, "Results/Final_MergedData.csv") ### THE FINAL DATASET

#######################################################################################################################
