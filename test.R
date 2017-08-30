#Quality Elements:
#Phytoplankton
#SAV
#Benthic Fauna
#Fish
#Supporting elements
rm(list=ls())
source("Assessment.R")
source("ReadBounds.R")
source("ReadIndicatorParms.R")
source("CalculateIndicator.R")

library(tidyverse)
library(haven)
library(lme4)
library(lubridate)
library(shiny)
library(dplyr)

df<-read.table("data/data.txt", fileEncoding = "UTF-8", sep=";", stringsAsFactors=F, header=T)
df.wb<-read.table("data/waterbodies.txt", fileEncoding = "UTF-8", sep="\t", stringsAsFactors=F, header=T)
df<-df %>% left_join(select(df.wb,WaterbodyID,DistrictID), by=c("WB_ID"="WaterbodyID")) 
df$WB<-paste0(df$WB_ID," ",df$WB_name)  

df.indicators<-read.table("data/IndicatorList.txt", fileEncoding = "UTF-8", sep="\t", stringsAsFactors=F, header=T)


nSimMC=100

AssessmentResults<-Assessment(df,nsim=nSimMC)
df.resultsOverall<-AssessmentResults[[1]]
df.resultsQE<-AssessmentResults[[2]]
df.resultsInd<-AssessmentResults[[3]]

df.chl<-AssessmentResults[[4]]