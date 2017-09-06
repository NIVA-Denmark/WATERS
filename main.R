library(tidyverse)
library(haven)
library(lme4)
library(lubridate)
library(prodlim)
library(matrixStats)

library(extrafont)
loadfonts(quiet = TRUE)
theme_set(theme_bw(base_size = 12, base_family = "Open Sans"))

# source("R/extract_modis.R") # ATTENTION: takes few minutes

# Clean the workspace
rm(list = ls())

source("ReadMonitoringData.R")
source("ReadIndicatorParms.R")
source("CalculateIndicator.R")
source("CalculateIndicatorSupport.R")
source("IndicatorFunctions.R")

# Read data set for specific waterbody and period
df <- ReadMonitoringDataSMHI("data/Gullmarn_2007_2012.sas7bdat")
df <- ReadMonitoringDataSMHI("data/danafjord_2013_2016.sas7bdat")
df <- ReadMonitoringDataSMHI("data/danafjord_2001_2006.sas7bdat")
df <- ReadMonitoringDataSMHI("data/byfjorden_2007_2012.sas7bdat")
# Read covariance parameters for the indicator
parmlist <- ReadParms_chla()
variance_list <- list(V_station=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "stati(vandom*period)"],V_obspoint=0,
                      V_year=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "year(vandomr*period)"],V_yearmonth=0,
                      V_tempres=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "Residual"],
                      V_stationyear=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "stat*year(vand*peri)"],V_stationmonth=0,
                      V_institution=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "proevetager"])
# Calculate the indicator
CalculateIndicator("Chla",df,RefCond_sali,variance_list,2007,2012)

# type with no salinity correction, e.g. type=6
RefCond_sali <- c(rep(0.9,36))
# type with salinity correction, e.g. type=8
RefCond_sali <- c(15.7,12.4,9.5,6.9,4.8,3.0,1.7,rep(1.29,29))
CalculateIndicator("ChlaEQR",df,RefCond_sali,variance_list,2007,2012)

# type with salinity correction, e.g. type=8 in summer
RefCond_sali <- c(56,50,43,37,31,24,18,rep(15,29))
CalculateIndicator("TNsummer",df,RefCond_sali,variance_list,2007,2012)

# type with salinity correction, e.g. type=8 in winter
RefCond_sali <- c(56,50,44,38,32,26,20,rep(17,29))
CalculateIndicator("TNwinter",df,RefCond_sali,variance_list,2007,2012)

# Calculate the indicator for Secchi depth
CalculateIndicator("Secchi",df,RefCond_sali,variance_list,2007,2012)

# type with salinity correction, e.g. type=8
RefCond_sali <- c(1.1,1.4,1.7,2.2,3.1,4.5,7.5,rep(10,29))
CalculateIndicator("SecchiEQR",df,RefCond_sali,variance_list,2007,2012)


CalculateIndicator("Chla",df,RefCond_sali,variance_list)
RefCond_sali <- c(15.7,12.4,9.5,6.9,4.8,3.0,1.7,rep(1.29,29))
CalculateIndicator("ChlaEQR",df,RefCond_sali,variance_list,2007,2012)

RefCond_sali <- c(56,50,44,38,32,26,20,rep(17,29))
CalculateIndicator("TNwinter",df,RefCond_sali,variance_list,2007,2012)
CalculateIndicator("TNwinterEQR",df,RefCond_sali,variance_list,2007,2012)

# BQI indicator
# Read data set for specific waterbody and period
variance_list <- list(V_station=2.5,V_obspoint=0.64,
                      V_year=0.10,V_yearmonth=0,
                      V_tempres=0.5,
                      V_stationyear=0.63,V_stationmonth=0,
                      V_institution=0.5)
df <- ReadMonitoringDataSMHI("data/halse_2007_2012.sas7bdat")
df <- ReadMonitoringDataSMHI("data/marstrand_2007_2012.sas7bdat")
CalculateIndicator("BQI",df,RefCond_sali,variance_list,2007,2012)


