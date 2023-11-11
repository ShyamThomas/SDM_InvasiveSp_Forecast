library(tidyverse)
library(miceadds)
library(ggplot2)
library(randomForest)
library(reshape2)


####################################################################################################################
######### R script: Predicting EWM occurrence with the trained Random forest (RF) & GAM SDMs for current   #########
########  and future water temperature conditions.                                                         #########
####################################################################################################################

### List all the models
RFmodels = list.files(path ="Data/TrainData/", pattern = "*RF.Rdata$")
### Revised GAMs
GAM_minK_models = list.files(path ="Data/TrainData/", pattern = "*s3.r.t3.Rdata$")
GAM_bestK_models = list.files(path ="Data/TrainData/", pattern = "*s3.r.t7.Rdata$")


### Read the saved future temperature data for years 2040-2060
EWM.futr.data=read_csv("Data/EWM.prsabs40to60_AllGCMs_v2.csv")
colnames(EWM.futr.data)
### Filter and subset the data by 5 GCMs
EWM.test.data_ACCESS.WtrTemp=EWM.futr.data%>%select('avg_secchi','roaddensity_density_mperha','ACCESS.avg.ann.gdd')
EWM.test.data_ACCESS.WtrTemp
write_csv(EWM.test.data_ACCESS.WtrTemp, "Data/TestData/ForecastData/EWM.forecast.data_ACCESS.WtrTemp.csv")

EWM.test.data_MIROC5.WtrTemp=EWM.futr.data%>%select('avg_secchi','roaddensity_density_mperha','MIROC5.avg.ann.gdd')
EWM.test.data_MIROC5.WtrTemp
write_csv(EWM.test.data_MIROC5.WtrTemp, "Data/TestData/ForecastData/EWM.forecast.data_MIROC5.WtrTemp.csv")

EWM.test.data_IPSL.WtrTemp=EWM.futr.data%>%select('avg_secchi','roaddensity_density_mperha','IPSL.avg.ann.gdd')
EWM.test.data_IPSL.WtrTemp
write_csv(EWM.test.data_IPSL.WtrTemp, "Data/TestData/ForecastData/EWM.forecast.data_IPSL.WtrTemp.csv")


EWM.test.data_GFDL.WtrTemp=EWM.futr.data%>%select('avg_secchi','roaddensity_density_mperha','GFDL.avg.ann.gdd')
EWM.test.data_GFDL.WtrTemp
write_csv(EWM.test.data_GFDL.WtrTemp, "Data/TestData/ForecastData/EWM.forecast.data_GFDL.WtrTemp.csv")

EWM.test.data_MRI.WtrTemp=EWM.futr.data%>%select('avg_secchi','roaddensity_density_mperha','MRI.avg.ann.gdd')
EWM.test.data_MRI.WtrTemp
write_csv(EWM.test.data_MRI.WtrTemp, "Data/TestData/ForecastData/EWM.forecast.data_MRI.WtrTemp.csv")

################################################################################################ 
### List all test data files (i.e future water temperature)
TestDataFiles = list.files(path ="Data/TestData/ForecastData", pattern = "*.csv$")
TestDataFiles

#### Start with GAMs with  min K value and make predictions for future temps

new.data=list()

for (TestData in TestDataFiles){
  Test = read.csv(paste("Data/TestData/ForecastData/",TestData, sep=""))
  #Test.preds=Test[,c(5,11,12)]
  colnames(Test)[3]=paste(str_sub(TestData, 19,-13),".avg.ann.gdd", sep="")
  new.data[[TestData]]=Test
}

fut.predictions=list()

for (i in 1:5){
  load.Rdata(paste("Data/TrainData/",GAM_minK_models[i], sep=""), "train.model")
  fut.predictions[[i]]=predict.gam(train.model, newdata = new.data[[i]], type = "response")
}
head(fut.predictions)

fut.preds.df=as.data.frame(do.call(cbind,fut.predictions))
head(fut.preds.df)
tail(fut.preds.df)

#### Repeat the same above codes with Current temperature(i.e. the training data)
TrainDataFiles=list.files(path ="Data/TrainData", pattern = "*.csv$")
TrainDataFiles

new.data=list()

for (TrainData in TrainDataFiles){
  Train = read.csv(paste("Data/TrainData/",TrainData, sep=""))
  #Train.preds=Train[,c(5,11,12)]
  new.data[[TrainData]]=Train
}

curr.predictions=list()

for (i in 1:5){
  load.Rdata(paste("Data/TrainData/",GAM_minK_models[i], sep=""), "train.model")
  curr.predictions[[i]]=predict(train.model, newdata = new.data[[i]],type = "response")
}
head(curr.predictions)
curr.preds.df=as.data.frame(do.call(cbind,curr.predictions))
head(curr.preds.df)
tail(curr.preds.df)

ModelNames=c("ACCESS","GFDL","IPSL","MIROC5","MRI")
ModelNames
colnames(fut.preds.df)=ModelNames
colnames(fut.preds.df)
colnames(curr.preds.df)=ModelNames
colnames(curr.preds.df)

fut.preds.df$DOWLKNUM=EWM.futr.data$DOWLKNUM
curr.preds.df$DOWLKNUM=EWM.futr.data$DOWLKNUM

write_csv(fut.preds.df, "Results/GAM.minK_Fut.Predictions.csv")
write_csv(curr.preds.df, "Results/GAM.minK_Curr.Predictions.csv")

### Box plot comparing risk current and future
curr.preds.df$Period=rep("Current",578)
fut.preds.df$Period=rep("Future",578)
currANDfut_preds=bind_rows(curr.preds.df,fut.preds.df)
currANDfut_preds_noDOW=currANDfut_preds[,-6]
head(currANDfut_preds_noDOW)
currANDfut_preds.melt=melt(currANDfut_preds_noDOW)
head(currANDfut_preds.melt)

### Finally, plot the predicted EWM risk estimates for current and future temperatures
InvasionRisk_Plot_GAM_minK=ggplot(currANDfut_preds.melt, aes(x=variable, y=value, fill=Period))+
  geom_boxplot()+ylab("Invasion Risk")+xlab("GCMs")
InvasionRisk_Plot_GAM_minK

###############################################################################
######### Repeat again for best K value based GAMs
new.data=list()

for (TestData in TestDataFiles){
  Test = read.csv(paste("Data/TestData/ForecastData/",TestData, sep=""))
  #Test.preds=Test[,c(5,11,12)]
  colnames(Test)[3]=paste(str_sub(TestData, 19,-13),".avg.ann.gdd", sep="")
  new.data[[TestData]]=Test
}

fut.predictions=list()

for (i in 1:5){
  load.Rdata(paste("Data/TrainData/",GAM_bestK_models[i], sep=""), "train.model")
  fut.predictions[[i]]=predict.gam(train.model, newdata = new.data[[i]], type = "response")
}
head(fut.predictions)

fut.preds.df=as.data.frame(do.call(cbind,fut.predictions))
head(fut.preds.df)
tail(fut.preds.df)

######## Repeat the same above codes with Current temperature(i.e. the training data)
TrainDataFiles=list.files(path ="Data/TrainData", pattern = "*.csv$")
TrainDataFiles

new.data=list()

for (TrainData in TrainDataFiles){
  Train = read.csv(paste("Data/TrainData/",TrainData, sep=""))
  #Train.preds=Train[,c(5,11,12)]
  new.data[[TrainData]]=Train
}

curr.predictions=list()

for (i in 1:5){
  load.Rdata(paste("Data/TrainData/",GAM_bestK_models[i], sep=""), "train.model")
  curr.predictions[[i]]=predict(train.model, newdata = new.data[[i]],type = "response")
}
head(curr.predictions)
curr.preds.df=as.data.frame(do.call(cbind,curr.predictions))
head(curr.preds.df)
tail(curr.preds.df)

ModelNames=c("ACCESS","GFDL","IPSL","MIROC5","MRI")
ModelNames
colnames(fut.preds.df)=ModelNames
colnames(fut.preds.df)
colnames(curr.preds.df)=ModelNames
colnames(curr.preds.df)

fut.preds.df$DOWLKNUM=EWM.futr.data$DOWLKNUM
curr.preds.df$DOWLKNUM=EWM.futr.data$DOWLKNUM

write_csv(fut.preds.df, "Results/GAM.bestK_Fut.Predictions.csv")
write_csv(curr.preds.df, "Results/GAM.bestK_Curr.Predictions.csv")

### Box plot comparing risk current and future
curr.preds.df$Period=rep("Current",578)
fut.preds.df$Period=rep("Future",578)
currANDfut_preds=bind_rows(curr.preds.df,fut.preds.df)
currANDfut_preds_noDOW=currANDfut_preds[,-6]
head(currANDfut_preds_noDOW)
currANDfut_preds.melt=melt(currANDfut_preds_noDOW)
head(currANDfut_preds.melt)

### Finally, plot the predicted EWM risk estimates for current and future temperatures
InvasionRisk_Plot_GAM_bestK=ggplot(currANDfut_preds.melt, aes(x=variable, y=value, fill=Period))+
  geom_boxplot()+ylab("Invasion Risk")+xlab("GCMs")
InvasionRisk_Plot_GAM_bestK

############ Clean up and give columns names

ModelNames=c("ACCESS","GFDL","IPSL","MIROC5","MRI")
ModelNames
colnames(fut.preds.df)=ModelNames
colnames(fut.preds.df)
colnames(curr.preds.df)=ModelNames
colnames(curr.preds.df)

fut.preds.df$DOWLKNUM=EWM.futr.data$DOWLKNUM
curr.preds.df$DOWLKNUM=EWM.futr.data$DOWLKNUM
write_csv(fut.preds.df, "Results/GAM.k3_Fut.Predictions.csv")
write_csv(curr.preds.df, "Results/GAM.k3_Curr.Predictions.csv")

curr.preds.df$Period=rep("Current",578)
fut.preds.df$Period=rep("Future",578)
currANDfut_preds=bind_rows(curr.preds.df,fut.preds.df)
currANDfut_preds_noDOW=currANDfut_preds[,-6]
head(currANDfut_preds_noDOW)
currANDfut_preds.melt=melt(currANDfut_preds_noDOW)
head(currANDfut_preds.melt)

### Finally, plot the predicted EWM risk estimates for current and future temperatures
InvasionRisk_Plot_GAM_k3=ggplot(currANDfut_preds.melt, aes(x=variable, y=value, fill=Period))+
  geom_boxplot()+ylab("Invasion Risk")+xlab("GCMs")
InvasionRisk_Plot_GAM_k3

################################################################################################ 
#### RF model predictions
### Loop all test data files to make a new data file with just the three predictors
### Then load each RF model and iterate model predictions for each new future dataset
new.data=list()

for (TestData in TestDataFiles){
  Test = read.csv(paste("processed_data/TestData/ForecastData/",TestData, sep=""))
  Test.preds=Test[,c(1:3)]
  colnames(Test.preds)[3]=paste(str_sub(TestData, 19,-13),".avg.ann.gdd", sep="")
  new.data[[TestData]]=Test.preds
}

fut.predictions=list()

for (i in 1:5){
  load.Rdata(paste("processed_data/TrainData/",RFmodels[i], sep=""), "train.model")
  fut.predictions[[i]]=predict(train.model, newdata = new.data[[i]])
}

head(fut.predictions)
fut.preds.df=as.data.frame(do.call(cbind,fut.predictions))
head(fut.preds.df)
write_csv(fut.preds.df, "Results/RF_Futr.Predictions.csv")

###########
#### Repeat the same above codes with Current temperature(i.e. the training data)
TrainDataFiles=list.files(path ="processed_data/TrainData", pattern = "*.csv$")
TrainDataFiles

new.data=list()

for (TrainData in TrainDataFiles){
  Train = read.csv(paste("processed_data/TrainData/",TrainData, sep=""))
  Train.preds=Train[,c(2:4)]
  new.data[[TrainData]]=Train.preds
}

curr.predictions=list()

for (i in 1:5){
  load.Rdata(paste("processed_data/TrainData/",RFmodels[i], sep=""), "train.model")
  curr.predictions[[i]]=predict(train.model, newdata = new.data[[i]])
}
head(curr.predictions)
curr.preds.df=as.data.frame(do.call(cbind,curr.predictions))
head(curr.preds.df)
write_csv(curr.preds.df, "Results/RF_Curr.Predictions.csv")

### Read the saved model predictions, index by DOW ids, and summarize predictions by GCMs
curr.preds.df=read_csv("Results/Curr.Predictions.csv")
curr.preds.df
fut.preds.df=read_csv("Results/Futr.Predictions.csv")
fut.preds.df
fut.preds.df$DOWLKNUM=EWM.futr.data$DOWLKNUM
curr.preds.df$DOWLKNUM=EWM.futr.data$DOWLKNUM
ModelNames=c("ACCESS","GFDL","IPSL","MIROC5","MRI")
ModelNames
colnames(fut.preds.df)=ModelNames
colnames(fut.preds.df)
colnames(curr.preds.df)=ModelNames
colnames(curr.preds.df)

curr.preds.df$Period=rep("Current",578)
fut.preds.df$Period=rep("Future",578)
currANDfut_preds=bind_rows(curr.preds.df,fut.preds.df)
currANDfut_preds.melt=melt(currANDfut_preds)
head(currANDfut_preds.melt)

### Finally, plot the predicted EWM risk estimates for current and future temperatures
InvasionRisk_Plot=ggplot(currANDfut_preds.melt, aes(x=variable, y=value, fill=Period))+
  geom_boxplot()+ylab("Invasion Risk")+xlab("GCMs")
InvasionRisk_Plot

###################################################################################################################### 
###################################################################################################################### 

