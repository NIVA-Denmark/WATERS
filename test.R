

rm(list=ls())
source("IndicatorFunctions.R")
source("ReadIndicatorType.R")
source("CalculateIndicatorSupport.R")
source("Assessment.R")
source("ReadBounds.R")
source("ReadIndicatorParms.R")
source("ReadMonitoringData.R")


library(tidyverse)
library(haven)
library(lme4)
library(lubridate)
library(shiny)
library(dplyr)
library(prodlim)

df<-read.table("data/data.txt", fileEncoding = "UTF-8", sep=";", stringsAsFactors=F, header=T)
df<-filter(df,!is.na(sali))
#df <- ReadMonitoringDataSMHI("data/danafjord_2001_2006.sas7bdat")
#df <- ReadMonitoringDataSMHI("data/danafjord_2007_2012.sas7bdat")
#df <- ReadMonitoringDataSMHI("data/danafjord_2013_2016.sas7bdat")
#df <- ReadMonitoringDataSMHI("data/byfjorden_2007_2012.sas7bdat")


df.wb<-read.table("data/waterbodies.txt", fileEncoding = "UTF-8", sep="\t", stringsAsFactors=F, header=T)
df<-df %>% left_join(select(df.wb,WaterbodyID,DistrictID), by=c("WB_ID"="WaterbodyID")) 
df$WB<-paste0(df$WB_ID," ",df$WB_name)  
df$obspoint<-df$station

# wblist<-distinct(df,WB,typology)
# wbcount<-nrow(wblist)
# df <- df %>% filter(WB == wblist$WB[1])

df.indicators<-read.table("data/IndicatorList.txt", fileEncoding = "UTF-8", sep="\t", stringsAsFactors=F, header=T)

nSimMC=100
source("Assessment.R")

# Start the clock!
ptm <- proc.time()

AssessmentResults<-Assessment(df,nsim=nSimMC)

proc.time() - ptm

IndicatorResults<-AssessmentResults
save(IndicatorResults,file="ExampleIndResults.Rda")
rm(IndicatorResults)
rm(AssessmentResults)

# 
# nSimMC
# df.resultsOverall<-AssessmentResults[[1]]
# df.resultsQE<-AssessmentResults[[2]]
# df.resultsInd<-AssessmentResults[[3]]
# 
# df.chl<-AssessmentResults[[4]]

load("C:/Data/GitHub/Waters/code/WATERS/ExampleIndResults.Rda")

IndicatorResults <- IndicatorResults %>% 
  select(WB,Type,Indicator,sim,Unit,Ref,HG,GM,MP,PB,Min,Value,ClassID,Class,EQR,Code) %>%
  left_join(df.indicators)
test<-IndicatorResults%>%filter(sim==1)

SubelementCount<-IndicatorResults %>%
  filter(!is.na(Value)) %>%
  group_by(WB,Type,Quality.element,Quality.subelement,sim) %>%
  summarise(IndCount=n()) %>%
  mutate(IndWeight=1/IndCount) %>%
  select(-IndCount)

IndicatorResults <- IndicatorResults %>% 
  left_join(SubelementCount)

QEresults<-IndicatorResults %>%
  mutate(EQRwt=EQR*IndWeight) %>%
  group_by(WB,Type,Quality.element,sim) %>%
  summarise(sumEQR=sum(EQRwt,na.rm=T),sumWt=sum(IndWeight,na.rm=T),EQR=sum(EQR*IndWeight)/sum(IndWeight,na.rm=T)) %>%
  mutate(EQRcheck=sumEQR/sumWt)
