#######################################################################################################
###################### R scripts of all the figures that went into appendix  ##########################
#######################################################################################################

#### Appendix figures
#### Figure S1: Maps showing EWM distribution and lake temperatures in GDD
library(Matrix)
library(tidyverse)
library(here)
library(reshape2)
library(sf)

Minn.sf
ggplot(data=Minn.sf)+geom_sf()

EWM.clmchng_data=read_csv("Results/Final_MergedData.csv") 
EWM.clmchng_data
EWM.climchng.Meantemps=EWM.clmchng_data%>%select(c(1,3,4,6),ends_with(c("Curr","Futr")))%>%rowwise%>%mutate(
  MeanCurr=mean(c_across(ends_with("Curr"))),
  MeanFutr=mean(c_across(ends_with("Futr"))),
  PerChange=((MeanFutr-MeanCurr)/MeanCurr)*100
)
EWM.climchng.Meantemps

EWM.climchng.Meantemps_sf=st_as_sf(EWM.climchng.Meantemps, coords = c("LON", "LAT"))
EWM.climchng.Meantemps_sf=st_set_crs(EWM.climchng.Meantemps_sf, 32615)
EWM.climchng.Meantemps_sf
EWM.climchng.Meantemps_WGS=st_transform(EWM.climchng.Meantemps_sf, crs=4326)
EWM.climchng.Meantemps_WGS2=EWM.climchng.Meantemps_WGS%>%mutate(EWM=recode(EWMSTATUS, "0" = "Abs", "1" = "Prs"))
EWM.climchng.Meantemps_WGS2


CurrGDDmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=MeanCurr))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Current GDD")+
  theme(text=element_text(size=16))

FutrGDDmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=MeanFutr))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Future GDD")+
  theme(text=element_text(size=16))

ChangeGDDmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=PerChange))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("% Increase in GDD")+
  theme(text=element_text(size=16))

library("patchwork")
CurrGDDmap|FutrGDDmap|ChangeGDDmap
ggsave("MnsptAppndx.S1_Fig2.png", path="Figures/", device="png",width=12, height=4, dpi=900)
#########################################################################################################

### Figure S2: Correlogram of the 3 predictors
library(PerformanceAnalytics)

EWM.clmchng_data$AvgCurrTemp=rowMeans(EWM.clmchng_data[, c(9:13)])
Corr.data=EWM.clmchng_data[,c(7,8,58)]
chart.Correlation(Corr.data, histogram=TRUE, pch=19)
ggsave("MnsptAppndx.S1_Fig4.png", path="Figures/", device="png",width=10, height=12, dpi=900)

######################################################################################################
### Figure S3: Maps showing the distribution of the response variable and predictors
EWMmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=as.factor(EWM)))+theme_light()+ 
  ggtitle("M. spicatum distribution")+
  scale_color_viridis_d(option = "turbo", alpha = 0.5)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+
  theme(text=element_text(size=16))

CurrGDDmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=MeanCurr))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Current GDD")+
  theme(text=element_text(size=16))

EWM.secchi.roads=EWM.clmchng_data%>%select(c(1,3,4,6:8))
EWM.secchi.roads_sf=st_as_sf(EWM.secchi.roads, coords = c("LON", "LAT"))
EWM.secchi.roads_sf=st_set_crs(EWM.secchi.roads_sf, 32615)
EWM.secchi.roads_WGS=st_transform(EWM.secchi.roads_sf, crs=4326)
EWM.secchi.roads_WGS2=EWM.secchi.roads_WGS%>%mutate(EWM=recode(EWMSTATUS, "0" = "Abs", "1" = "Prs"))

SecchiMap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.secchi.roads_WGS2, aes(color=avg_secchi))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Secchi depth")+
  theme(text=element_text(size=16))

RoadsMap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.secchi.roads_WGS2, aes(color=roaddensity_density_mperha))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Road density")+
  theme(text=element_text(size=16))

(EWMmap|CurrGDDmap)/(RoadsMap|SecchiMap)
ggsave("MnsptAppndx.S1_Fig4.png", path="Figures/", device="png",width=10, height=12, dpi=900)

#########################################################################################################
### Figure S5: Results of ANOVA
EWM.CurrFutr.TempPreds.Doms.csv=read_csv(here("Data", "EWM.CurrFutr.TempPreds.Doms.csv"))

MergedCurrData=left_join(EWM.CurrFutr.Preds.Doms[,1:15],curr.preds.GAM_minK[,1:6], by="DOWLKNUM")%>%
  left_join(.,curr.preds.GAM_bestK[,1:6], by="DOWLKNUM")%>%
  left_join(., RF_curr.preds[,-6],by="DOWLKNUM")
MergedCurrData  

AllResults_CurrPredslong=MergedCurrData[,c(1,3,4,16:30)]%>%pivot_longer(
  cols=ends_with("CurrPred"),
  values_to = "CurrPreds")
AllResults_CurrPredslong

SDM_string=c(rep("GAMk3",5),rep("GAMk10",5),rep("RF",5))
SDM_string

AllResults_CurrPredslong_newFacts=AllResults_CurrPredslong%>%mutate(SDM=rep(SDM_string, 578), 
                                                                    GCM=rep(c("ACCESS","GFDL","IPSL","MIROC5","MRI"), 1734))
AllResults_CurrPredslong_newFacts

### Now lets do the  Future predictions data
MergedFutrData=left_join(EWM.CurrFutr.Preds.Doms[,c(1:10,16:25)],futr.preds.GAM_minK[,1:6], by="DOWLKNUM")%>%
  left_join(.,futr.preds.GAM_bestK[,1:6], by="DOWLKNUM")%>%
  left_join(., RF_futr.preds[,-6],by="DOWLKNUM")
MergedFutrData

AllResults_FutrPredslong=MergedFutrData[,c(1,3,4,21:35)]%>%pivot_longer(
  cols=ends_with("FutrPred"),
  values_to = "FutrPreds")

AllResults_FutrPredslong_newFacts=AllResults_FutrPredslong%>%mutate(SDM=rep(SDM_string, 578), 
                                                                    GCM=rep(c("ACCESS","GFDL","IPSL","MIROC5","MRI"), 1734))
AllResults_FutrPredslong_newFacts

final_changeinrisk=bind_cols(AllResults_CurrPredslong_newFacts[,-4],AllResults_FutrPredslong_newFacts[,c(1,5)])%>%
  mutate(ChangeInRisk=FutrPreds-CurrPreds)
final_changeinrisk

Sum.Sq_list=final_changeinrisk%>%group_split(DOWLKNUM...1)%>%
  map(~ summary(aov(ChangeInRisk ~ SDM+GCM, data = .))[[1]]$'Sum Sq')
class(Sum.Sq_list)
Sum.Sq_DF=as.data.frame(do.call(rbind, Sum.Sq_list))
colnames(Sum.Sq_DF)=c("SDM", "GCM", "Residual")
Sum.Sq_prpn=Sum.Sq_DF/rowSums(Sum.Sq_DF)
Sum.Sq_prpn$DOWLKNUM=curr.preds.df$DOWLKNUM
Sum.Sq_prpn$LON=EWM.CurrFutr.Preds.Doms$LON
Sum.Sq_prpn$LAT=EWM.CurrFutr.Preds.Doms$LAT
head(Sum.Sq_prpn)
boxplot(Sum.Sq_prpn[,1:3])

AvgTemp.Domain=EWM.CurrFutr.Preds.Doms[,1:4]
AvgTemp.Domain
AvgTemp.Domain$AvgFutrTemp=rowMeans(EWM.CurrFutr.Preds.Doms[,16:20])
AvgTemp.Domain
range(AvgTemp.Domain$AvgFutrTemp)
AvgTemp.Domain=AvgTemp.Domain%>%mutate(Domain=case_when(AvgFutrTemp < 2220 ~ 'Analog', 
                                                        AvgFutrTemp > 2220 ~ 'Non-analog'))
Sum.Sq_prpn.Dom=left_join(Sum.Sq_prpn,AvgTemp.Domain[,c(1,5,6)], by="DOWLKNUM")
head(Sum.Sq_prpn.Dom)

Sum.Sq_prpn.Dom_melt=melt(Sum.Sq_prpn.Dom, id=c("DOWLKNUM", "LON", "LAT", "Domain", "AvgFutrTemp"))
head(Sum.Sq_prpn.Dom_melt)
ggplot(Sum.Sq_prpn.Dom_melt,aes(x=variable, y=value, fill=Domain))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(0.2), alpha=0.2)+ylab("Propn. of Total Sum of Squares")+
  xlab(" ")+ theme(text=element_text(size=16))+theme(legend.position = c(0.85, 0.85))

ggsave("MnsptAppndx.S1_Fig5.png", path="Figures/", device="png",width=6, height=4.5)

Minn.sf=st_transform(Minn.sf, crs=32615)
Minn.sf
Sum.Sq_prpn_sf=st_as_sf(Sum.Sq_prpn,coords=c("LON", "LAT"),crs=32615)
Sum.Sq_prpn_sf

ggplot(Minn.sf)+geom_sf()+geom_sf(data=Sum.Sq_prpn_sf, aes(col=SDM))+scale_color_viridis_c(name="Variance \n(SDMs)")+theme_minimal()+
  theme(text=element_text(size=16))+theme(legend.position = c(0.8, 0.4))
ggsave("PropnTSS_SDM_map.png", path="Figures/", device="png",width=4.5, height=6)
ggplot(Minn.sf)+geom_sf()+geom_sf(data=Sum.Sq_prpn_sf, aes(col=GCM))+scale_color_viridis_c(name="Propn. of TSS \n(GCMs)")+theme_minimal()+
  theme(text=element_text(size=16))+theme(legend.position = c(0.85, 0.4))


#########################################################################################################

###Figure S5: 

left_join(EWM.domains, Fut.Brm.Preds_EWMGeoIndex)%>%
  ggplot(., aes(mean.fut.preds, color= Domain, ..scaled..))+geom_density()+
  xlab("Median of predicted posterior estimates")+ylab("")+theme(legend.position = c(0.125,0.875))
ggsave("MnsptAppndx.S1_Fig6a.png", path="Figures/", device="png",width=6, height=4.5)

left_join(EWM.domains, Fut.Brm.Preds_EWMGeoIndex)%>%
  ggplot(., aes(sd.fut.preds, color= Domain, ..scaled..))+geom_density()+
  xlab("SD of predicted posterior estimates")+ylab("")+theme(legend.position = "NA")
ggsave("MnsptAppndx.S1_Fig6b.png", path="Figures/", device="png",width=6, height=4.5)


#####################################################################################################################
#####################################################################################################################
