rm(list=ls())
source("Assessment.R")
source("ReadBounds.R")
source("ReadIndicatorParms.R")
source("IndicatorFunctions.R")
source("ReadIndicatorType.R")
source("CalculateIndicatorSupport.R")
source("ReadMonitoringData.R")


library(tidyverse)
library(haven)
library(lme4)
library(lubridate)
library(shiny)
library(dplyr)
library(prodlim)


#df<-read.table("data/data.txt", fileEncoding = "UTF-8", sep=";", stringsAsFactors=F, header=T)
df <- ReadMonitoringDataSMHI("data/danafjord_2013_2016.sas7bdat")

df.wb<-read.table("data/waterbodies.txt", fileEncoding = "UTF-8", sep="\t", stringsAsFactors=F, header=T)
df<-df %>% left_join(select(df.wb,WaterbodyID,DistrictID), by=c("WB_ID"="WaterbodyID")) 
df$WB<-paste0(df$WB_ID," ",df$WB_name)  
df$obspoint<-df$station

df.indicators<-read.table("data/IndicatorList.txt", fileEncoding = "UTF-8", sep="\t", stringsAsFactors=F, header=T)

nSimMC=100



df.bounds<-ReadBounds()
df.indicators<-ReadIndicatorType()

# Read general parameters for the indicator
parmlist <- ReadParms_chla()
variance_list <- list(V_station=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "stati(vandom*period)"],V_obspoint=0,
                      V_year=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "year(vandomr*period)"],V_yearmonth=0,
                      V_tempres=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "Residual"],
                      V_stationyear=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "stat*year(vand*peri)"],V_stationmonth=0,
                      V_institution=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "proevetager"])

#df<-filter(df,!is.na(chla))
RefCond_sali<-SalinityReferenceValues(df,df.bounds,"ChlaEQR",missing=50)
res.chlaEQR<-CalculateIndicator("ChlaEQR",df,RefCond_sali,variance_list,n_iter=100)
RefCond_sali<-SalinityReferenceValues(df,df.bounds,"TNwinterEQR",missing=50)
res<-CalculateIndicator("TNwinterEQR",df,RefCond_sali,variance_list,n_iter=100)
RefCond_sali<-SalinityReferenceValues(df,df.bounds,"TNsummerEQR",missing=50)
res<-CalculateIndicator("TNsummerEQR",df,RefCond_sali,variance_list,n_iter=100)

df.temp<-data.frame(Mean=res.chlaEQR$period$mean,StdErr=res.chlaEQR$period$stderr)
df.temp$Indicator<-"chlaEQR"
res.ind.temp<-df.temp
