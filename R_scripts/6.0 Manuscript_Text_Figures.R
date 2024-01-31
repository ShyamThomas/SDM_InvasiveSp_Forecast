library(tidyverse) 
library(sf)
library(patchwork)
library(gratia)
library(mgcv)
library(cowplot)
library(here)
library(miceadds)
library(randomForest)
library(pdp)


########################################################################################################
##### R scripts of all the figures that went into the manuscript main text and appendix ###############
#######################################################################################################

### The final data with all SDM results, lake-level covariates, and Minnesota map shapefile
EWM.clmchng_data=read_csv("Results/Final_MergedData.csv") 
Minn.sf=read_sf(here("Data", "GIS_Data"), layer="Minn.map")

############################  
### Figure 1: Boxplots comparing predicted EWM invasion risk under current and future water temperature/GDD
EWM.clmchng_data ### from Line 6 above
EWM.climchng.preds_wide=EWM.clmchng_data%>%select(c(1,3,4),ends_with("Pred"))%>%
    pivot_longer(!c(DOWLKNUM,LON,LAT), names_to = "models", values_to = "preds")
EWM.clim.chng.preds=EWM.climchng.preds_wide%>%separate(models, c("GCM","SDM","PERIOD"))
EWM.clim.chng.preds

EWM.clim.chng.preds2=EWM.clim.chng.preds%>%mutate(SDM=recode(SDM, GAMminK="GAM (k=3)",GAMbestK="GAM (k=7)", 
                                                    RF="Random Forest"),
                                            PERIOD=recode(PERIOD,CurrPred="Current", FutrPred="Future"))
EWM.clim.chng.preds2


PeriodWise=EWM.clim.chng.preds%>%ggplot(aes(x = SDM, y = preds, fill = GCM)) +
            geom_boxplot(outlier.shape = NA) +
            facet_grid(~PERIOD) +
            theme(legend.position = "bottom")

PeriodWise+scale_fill_viridis_d(option="inferno") ### not using this plot version

SDMWise=EWM.clim.chng.preds2%>%ggplot(aes(x = PERIOD, y = preds, fill = GCM)) +
        geom_boxplot(outlier.shape = NA) +
        facet_grid(~SDM) +
        theme(legend.position = c(0.93,0.22), legend.box.background = element_rect(color = "black"))+
        xlab(" ")+ylab("M. spicatum habitat suitability")+
        theme(text = element_text(size=12))

SDMWise+scale_fill_viridis_d(option="inferno")
ggsave("Mnspt.Text.Fig1.tif", path="Figures/", device="tiff",width=9, height=4.5, dpi=600)

####################################################################################
### Figure 2: Arrow plots showing change in risk over time
EWM.clmchng_data ### from above
### Estimate change in  mean risk predictions
EWM.climchng.MeanPredsBySDM.Year=EWM.clmchng_data%>%select(c(1,3,4),ends_with("Pred"))%>%rowwise%>%mutate(
  MeanGAM_minK_Curr=mean(c_across(ends_with("GAMminK.CurrPred"))),
  MeanGAM_minK_Futr=mean(c_across(ends_with("GAMminK.FutrPred"))),
  
  MeanGAM_bestK_Curr=mean(c_across(ends_with("GAMbestK.CurrPred"))),
  MeanGAM_bestK_Futr=mean(c_across(ends_with("GAMbestK.FutrPred"))),
  
  MeanRF_Curr=mean(c_across(ends_with("RF.CurrPred"))),
  MeanRF_Futr=mean(c_across(ends_with("RF.FutrPred")))
)

EWM.clim.chng.PredsChange=EWM.climchng.MeanPredsBySDM.Year%>%select(c(1:3),starts_with("Mean"))%>%mutate(
                            GAM_minK_Change=MeanGAM_minK_Futr-MeanGAM_minK_Curr,
                            GAM_bestK_Change=MeanGAM_bestK_Futr-MeanGAM_bestK_Curr,
                            RF_Change=MeanRF_Futr-MeanRF_Curr)
EWM.clim.chng.PredsChange%>%View()

### Change in mean temp conditions
EWM.climchng.MeanTempByYear=EWM.clmchng_data%>%select(c(1,3,4),ends_with(c("Curr","Futr")))%>%rowwise%>%mutate(
  MeanTemp_Curr=mean(c_across(ends_with("Curr"))),
  MeanTemp_Futr=mean(c_across(ends_with("Futr")))
)

EWM.clim.chng.TempChange=EWM.climchng.MeanTempByYear%>%select(c(1:3),starts_with("Mean"))%>%rowwise%>%mutate(
                        TempDiff=MeanTemp_Futr-MeanTemp_Curr
)
        
EWM.clim.chng.Preds.Temp.Change=left_join(EWM.clim.chng.PredsChange,EWM.clim.chng.TempChange[,c(1,4:6)], by="DOWLKNUM")
EWM.clim.chng.Preds.Temp.Change ## the final table capturing all the changes

### Plot for GAM mink (s3.r.t3)
GAM_minK=ggplot()+
  geom_segment(data=EWM.clim.chng.Preds.Temp.Change, aes(x=MeanTemp_Curr, y=MeanGAM_minK_Curr, xend=MeanTemp_Curr+TempDiff, 
  yend=MeanGAM_minK_Curr+GAM_minK_Change,color=GAM_minK_Change), arrow=arrow(), size=0.5) +
  geom_point(data=EWM.clim.chng.Preds.Temp.Change, mapping=aes(x=MeanTemp_Curr, y=MeanGAM_minK_Curr), size=1, shape=21, fill="white")+
  scale_color_viridis_c()+xlab("Annual growing degree days (GDD)")+ylab("EWM invasion risk")+
  labs(colour="Change\nin risk")+geom_vline(xintercept = 2200, lty=2)+theme(text=element_text(size=20))+
  theme(legend.position = c(0.06,0.8))
GAM_minK

### GAM minK response curve inset plot

load.Rdata("Data/TrainData/EWM.train.data_MRI.WtrTempGAM_s3.r.t3.Rdata", "MRI.GAM_minK.model")

gdd.resp_mink=draw(MRI.GAM_minK.model, residuals = FALSE, select=2, rug=FALSE, ci_alpha = 0.00, title="")
GAM.minK.GDD_respcurv=gdd.resp_mink+xlab("GDD")+ylab(bquote(italic("M. spicatum")))+theme_classic(36)+ggtitle("")+theme(
axis.text.x = element_blank(),
axis.text.y = element_blank(),
axis.ticks = element_blank())
GAM.minK.GDD_respcurv+geom_line(size=2)
ggsave("GDD.resp_minK.tif", path="Figures/", device="tiff",width = 6, height = 6, dpi=600)


### Plot for GAM bestK s3.r.t7
GAM_bestK=ggplot()+
  geom_segment(data=EWM.clim.chng.Preds.Temp.Change, aes(x=MeanTemp_Curr, y=MeanGAM_bestK_Curr, xend=MeanTemp_Curr+TempDiff, 
  yend=MeanGAM_bestK_Curr+GAM_bestK_Change,color=GAM_bestK_Change), arrow=arrow(), size=0.5) +
  geom_point(data=EWM.clim.chng.Preds.Temp.Change, mapping=aes(x=MeanTemp_Curr, y=MeanGAM_bestK_Curr), size=1, shape=21, fill="white")+
  scale_color_viridis_c()+xlab("Annual growing degree days (GDD)")+ylab("EWM invasion risk")+
  labs(colour="Change\nin risk")+geom_vline(xintercept = 2200, lty=2)+theme(text=element_text(size=20))+
  theme(legend.position = c(0.06,0.8))
GAM_bestK

### GAM bestK (s3.r.t7) response curve inset plot
load.Rdata("Data/TrainData/EWM.train.data_MRI.WtrTempGAM_s3.r.t7.Rdata", "ACCESS.GAM_bestK.model")
gam.bestk=draw(ACCESS.GAM_bestK.model, residuals = FALSE, select=2, rug=FALSE, ci_alpha = 0.00, title="")
GAM.bestK.GDD_respcurv=gam.bestk+xlab("GDD")+ylab(bquote(italic("M. spicatum")))+theme_classic(24)+ggtitle("")+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())
GAM.bestK.GDD_respcurv+geom_line(size=2)
ggsave("GDD.resp_bestK.tif", path="Figures/", device="tiff",width = 6, height = 6, dpi=600)

### Plot for Random forest plot
RF=ggplot()+
geom_segment(data=EWM.clim.chng.Preds.Temp.Change, aes(x=MeanTemp_Curr, y=MeanRF_Curr, xend=MeanTemp_Curr+TempDiff,
yend=MeanRF_Curr+RF_Change,color=RF_Change),arrow=arrow(), size=0.5) +
  geom_point(data=EWM.clim.chng.Preds.Temp.Change, mapping=aes(x=MeanTemp_Curr, y=MeanRF_Curr), size=1, shape=21, fill="white")+
  scale_color_viridis_c()+xlab("Annual growing degree days (GDD)")+ylab("EWM invasion risk")+
  labs(colour="Change \nin risk")+geom_vline(xintercept = 2200, lty=2)+theme(text=element_text(size=20))+
  theme(legend.position = c(0.06,0.8))

### Random forest response curve inset plot
EWM.train.data_ACCESS.WtrTemp=read_csv("Data/TrainData/EWM.train.data_ACCESS.WtrTemp.csv")
EWM.train.data_ACCESS.WtrTemp
rf=randomForest(EWM.train.data_ACCESS.WtrTemp[,c(2:4)], EWM.train.data_ACCESS.WtrTemp$EWMSTATUS,importance=TRUE, ntree=5000, type="regression")
par.RF= pdp::partial(rf, pred.var = c("ACCESS.avg.ann.gdd"))

RF.GDD_respcurv=autoplot(par.RF)+xlab("GDD")+ylab(bquote(italic("M. spicatum")))+theme_classic(24)+theme(
axis.text.x = element_blank(),
axis.text.y = element_blank(),
axis.ticks = element_blank())
RF.GDD_respcurv+geom_line(size=2)
ggsave("GDD.resp_RF.tif", path="Figures/", device="tiff",width = 6, height = 6, dpi=600)

### Making one single multi-panel plot with a common legend
### Pivot the table from wide to long
EWM.clim.chng.Preds.Temp.Change2=EWM.clim.chng.Preds.Temp.Change%>%
  pivot_longer(cols= c(4,6,8), names_to = 'Curr_Models' , values_to = 'Curr.Preds')
t1=EWM.clim.chng.Preds.Temp.Change2[,c(1:3,10:14)]
EWM.clim.chng.Preds.Temp.Change3=EWM.clim.chng.Preds.Temp.Change%>%
  pivot_longer(cols= c(5,7,9), names_to = 'Futr_Models' , values_to = 'Futr.Preds')
t2=EWM.clim.chng.Preds.Temp.Change3[,c(1,13,14)]
EWM.clim.chng.Preds.Temp.Change4=EWM.clim.chng.Preds.Temp.Change%>%
  pivot_longer(cols= c(10:12), names_to = 'Chng_Models' , values_to = 'Chng.Preds')
t3=EWM.clim.chng.Preds.Temp.Change4[,c(1,13,14)]

t1$Futr.preds=t2$Futr.Preds
t1$Chng.preds=t3$Chng.Preds
final_table=t1%>%mutate(Model= recode(Curr_Models, MeanGAM_minK_Curr = "GAM (k=3)", MeanGAM_bestK_Curr = "GAM (k=7)",
                                      MeanRF_Curr = "Random Forest" ))
final_table%>%View()

### the final multi-panel plot using the final table in long format
AllModel_ChangeinRisk=ggplot()+
            geom_segment(data=final_table, aes(x=MeanTemp_Curr, y=Curr.Preds, 
            xend=MeanTemp_Curr+TempDiff,yend=Curr.Preds+Chng.preds,color=Chng.preds), 
            arrow=arrow(), size=0.5) +
              geom_point(data=final_table, mapping=aes(x=MeanTemp_Curr, y=Curr.Preds), size=1, shape=21, fill="white")+
                scale_color_viridis_c(limits = c(-0.4, 0.9),
                  breaks = c(-0.4,0.00, 0.50, 0.9),
                  labels = c(-0.4,0.00, 0.50, 0.9))+
                    xlab("Annual growing degree days (GDD)")+ylab(bquote(italic("M. spicatum ")*"habitat suitability"))+
                    labs(colour="Change\nin risk")+geom_vline(xintercept = 2200, lty=2)+
                    theme(text=element_text(size=12))+
                      facet_wrap(~Model)+theme_bw()+
                      theme(legend.position="right", legend.box.background = element_rect(color = "black"))

AllModel_ChangeinRisk
ggsave("AllChangeInRisk.tif", path="Figures/", device="tiff",width = 18, height = 12, dpi=600, units="cm")

##################################################################################################################################
### Figure 3: Final EWM invasion risk predictions and uncertainty maps, analog/non-analog domains
EWM.clim.chng.PredsChange ### from line 100 above
EWM.clim.chng.PredsChange=EWM.climchng.MeanPredsBySDM.Year%>%select(c(1:3),starts_with("Mean"))%>%mutate(
GAM_minK_Change=MeanGAM_minK_Futr-MeanGAM_k3_Curr,
GAM_bestK_Change=MeanGAM_bestK_Futr-MeanGAM_bestK_Curr,
RF_Change=MeanRF_Futr-MeanRF_Curr)

### Estimate change in EWM suitability status
EWM.clim.chng.ChangeStatus=EWM.clim.chng.PredsChange%>%select(1:3, ends_with("Change"))%>%
  mutate(GAM.minK_Change_Status=case_when(GAM_minK_Change < -0.1 ~ 'loss',
          GAM_minK_Change < 0.1 ~ 'no change',
          GAM_minK_Change < 0.85 ~ 'gain'), 

        GAM.bestK_Change_Status=case_when(GAM_bestK_Change < -0.1 ~ 'loss',
          GAM_bestK_Change < 0.1 ~ 'no change',
          GAM_bestK_Change < 0.85 ~ 'gain'),

        RF_Change_Status=case_when(RF_Change < -0.1 ~ 'loss',
          RF_Change < 0.1 ~ 'no change',
          RF_Change < 0.9 ~ 'gain'),
)

### Estimate consensus among models in change in EWM suitability status
StatusCount=EWM.clim.chng.ChangeStatus%>%select(ends_with("Status"))%>%rowwise()%>%
  do(data.frame(., StatusCount = n_distinct(unlist(.))))%>%pull(StatusCount)
EWM.clim.chng.ChangeStatus$StatusCount=StatusCount
EWM.clim.chng.ChangeStatus%>%mutate(Status=recode(StatusCount, "1"= "Increase", "2"= "Two", "3"="Uncertain"))

IncreasersUncertain_Lakes=EWM.clim.chng.ChangeStatus%>%filter(StatusCount ==1 | StatusCount>1)
IncreasersUncertain_Lakes
IncreasersUncertain_Lakes_sf=st_as_sf(IncreasersUncertain_Lakes, coords=c("LON","LAT"), crs=32615)
IncreasersUncertain_Lakes_sf_WGS=st_transform(IncreasersUncertain_Lakes_sf, crs=4326)
IncreasersUncertain_Lakes_sf_WGS_recoded=IncreasersUncertain_Lakes_sf_WGS%>%
  mutate(Status=recode(StatusCount, "1"= "Increase","2"="Uncertain", "3"="Uncertain"))

### Plot the change in invasion risk status (trajectory)
Fig.3a=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=IncreasersUncertain_Lakes_sf_WGS_recoded, aes(color=as.factor(Status),shape=as.factor(Status)), cex=4)+
  theme_light()+
    scale_color_viridis_d(option="turbo",alpha = 0.5)+theme(legend.title = element_blank())+
    theme(legend.position = c(0.85,0.4))+theme(text=element_text(size=12))+
      ggtitle("a) Future invasion risk trajectory")
Fig.3a
ggsave("Mnspt.Text.Fig3a.tif", path="Figures/", device="tiff",width=14, height=20,units="cm", dpi=900)

### Subset bayesian predictions from the best fitting GAM k=10 model 
EWM.futrpreds.bayes.GAM_bestK=EWM.clmchng_data%>%select(1,3,4,56,57)%>%
  rename(Estimate = mean.fut.preds, Variance = var.fut.preds)
EWM.futrpreds.bayes.GAM_bestK

EWM.futrpreds.bayes.GAM_bestK_sf=st_as_sf(EWM.futrpreds.bayes.GAM_bestK, coords = c("LON", "LAT"),crs=32615)
EWM.futrpreds.bayes.GAM_bestK_sf_WGS=st_transform(EWM.futrpreds.bayes.GAM_bestK_sf, crs=4326)
EWM.futrpreds.bayes.GAM_bestK_sf_WGS

### Plot the lake-level bayesian predictions and variance associated with it
Fig.3b=ggplot(data=Minn.sf)+geom_sf()+
          geom_sf(data=EWM.futrpreds.bayes.GAM_bestK_sf_WGS, aes(color=Estimate, size=Variance))+theme_light()+
          scale_color_viridis_c(option="turbo",alpha = 0.75)+
          guides(size=guide_legend(override.aes=list(shape=1,size=c(2,4,6,8))))+
          theme(legend.position = c(0.85,0.375))+theme(text=element_text(size=12))+
          ggtitle("b) Future invasion risk predictions")
Fig.3b
ggsave("Mnspt.Text.Fig3b.tif", path="Figures/", device="tiff",width=14, height=20,units="cm", dpi=900)

### Subset data to identify the temperature analog and non-analog domains
EWM.domains=EWM.climchng.MeanTempByYear%>%select(1:3,14,15)%>%
                mutate(Domain = case_when(MeanTemp_Futr < 2225 ~ 'Analog', MeanTemp_Futr > 2225 ~ 'Non-analog'))
EWM.domains.sf=st_as_sf(EWM.domains, coords=c("LON", "LAT"), crs=32615)
EWM.domains.sf.WGS=st_transform(EWM.domains.sf, crs=4326)

### Map the lakes in analog and non-analog domains
Fig.3c=ggplot(data=Minn.sf)+geom_sf()+
          geom_sf(data=EWM.domains.sf.WGS, aes(color=Domain,shape=Domain), cex=4)+theme_light()+
          scale_color_viridis_d(option="turbo",alpha = 0.5)+theme(legend.title = element_blank())+
          theme(legend.position = c(0.8,0.45))+theme(text=element_text(size=12))+ggtitle("c) Future temperature domains")
Fig.3c
ggsave("Mnspt.Text.Fig3c.tif", path="Figures/", device="tiff",width=14, height=20,units="cm", dpi=900)

Fig.3d=Fut.Brm.Preds%>%ggplot(.,)+geom_point(aes(avgGDD,mean.fut.preds, size=var.fut.preds), alpha=0.25)+
geom_smooth(aes(x=avgGDD,y=mean.fut.preds), method="loess", span=0.3,se=FALSE, col="black", linewidth=0.5)+
geom_smooth(aes(x= avgGDD, y=min.fut.preds),method="loess", span=0.3, se=FALSE, lty=2, linewidth=0.5)+
geom_smooth(aes(x= avgGDD, y=max.fut.preds), method="loess", span=0.3,se=FALSE, lty=2, linewidth=0.5)+
geom_vline(xintercept = 2200, lty=3, linewidth=1)+xlab("Growing degree days (GDD)")+
ylab("Future M. spicatum habitat suitability")+labs(size="Variance")+theme(legend.position = c(0.85,0.2))+
theme(text=element_text(size=12))+ggtitle("d) Predicted response curve")
Fig.3d
ggsave("Mnspt.Text.Fig3d.tif", path="Figures/", device="tiff",width=20, height=14,units="cm", dpi=900)

### The numbers of certainty and uncertainty in the change in suitability
IncreasersUncertain_Lakes_sf_WGS_recoded%>%st_drop_geometry()
Domain_IncUnc.Lakes=left_join(EWM.domains, IncreasersUncertain_Lakes, by="DOWLKNUM")

Domain_IncUnc.Lakes%>% group_by(Domain)%>%tally()
Domain_IncUnc.Lakes%>% group_by(Domain, Status)%>%tally()


#####################################################################################################################
#####################################################################################################################
