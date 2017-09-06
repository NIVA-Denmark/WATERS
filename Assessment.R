
#' Assessment
#' 
#' 
#' @param df A dataframe with monitoring data from the Swedish Monitoring program. 
#'   The dataframe should contain the
#'   following variables:
#'   
#'   \item{station}{An identifier for the monitoring station.} 
#'   \item{date}{Date of the sampling.} 
#'   \item{depth}{Depth of the observation in meter.} 
#'   \item{chla}{Chlorophyll a concentration in sample.} 
#'   \item{WB}{Waterbody} 
#'   
#' @param nsim Number of iterations for Monte Carlo simulation 
#'   
#' 
#' @examples
Assessment <-
  function(df.all,nsim=1000) {
    
    df.bounds<-ReadBounds()
    df.indicators<-ReadIndicatorType()
    
    df.all$typology<-gsub("SE_", "", df.all$typology)
    
    wblist<-distinct(df.all,WB,typology)
    wbcount<-nrow(wblist)
    
    # Read general parameters for the indicator
    parmlist <- ReadParms_chla()
    variance_list <- list(V_station=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "stati(vandom*period)"],V_obspoint=0,
                          V_year=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "year(vandomr*period)"],V_yearmonth=0,
                          V_tempres=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "Residual"],
                          V_stationyear=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "stat*year(vand*peri)"],V_stationmonth=0,
                          V_institution=parmlist$covparams_CumCover$Estimate[parmlist$covparams_CumCover$CovParm == "proevetager"])
    
    
    #cat(paste0("ALL: ",nrow(df.all),"\n"))
    for(iWB in 1:wbcount){
      df.temp<-df.all %>% filter(WB == wblist$WB[iWB])
      plist<-distinct(df.all,period)
      pcount<-nrow(plist)
      
      for(iPeriod in 1:pcount){
        df <- df.all %>% filter(WB == wblist$WB[iWB],period == plist$period[iPeriod])
        cat(paste0("WB: ",wblist$WB[iWB],"  Period: ",plist$period[iPeriod],"\n"))
        
        startyear<-as.numeric(substr(as.character(plist$period[iPeriod]),1,4))
        endyear<-as.numeric(substr(as.character(plist$period[iPeriod]),6,9))
        # Calculate the indicator
        #function(df,MonthInclude,variance_list,n_iter=1000)
        #res.chl<-CalculateIndicator_Chla(df,c(6,7,8),variance_list,n_iter=nsim)
        
        # Get salinity correction values for Chla EQR indicator and calculate indicator
        RefCond_sali<-SalinityReferenceValues(df,df.bounds,"chlaEQR",missing=0.9)
        res.chlaEQR<-CalculateIndicator("ChlaEQR",df,RefCond_sali,variance_list,startyear,endyear,n_iter=nsim)
        
        # Get salinity correction values for TNsummer indicator
        RefCond_sali<-SalinityReferenceValues(df,df.bounds,"TNsummer",missing=50)
        RefCond_sali <- c(56,50,43,37,31,24,18,rep(15,29))
        res.TNsummer<-CalculateIndicator("TNsummer",df,RefCond_sali,variance_list,startyear,endyear,n_iter=nsim)
        
        # Get salinity correction values for TNwinter indicator
        RefCond_sali<-SalinityReferenceValues(df,df.bounds,"TNwinter",missing=50)
        RefCond_sali <- c(56,50,44,38,32,26,20,rep(17,29))
        res.TNwinter<-CalculateIndicator("TNwinter",df,RefCond_sali,variance_list,startyear,endyear,n_iter=nsim)
        
        # Calculate the indicator for Secchi depth
        res.Secchi<-CalculateIndicator("Secchi",df,c(6,7,8),variance_list,startyear,endyear,n_iter=nsim)
        #CalculateIndicator <-
        #   function(Indicator,df,RefCond_sali,var_list,n_iter=1000,confidence_lvl=0.95)
        
        
        #Combine results of different indicators - Mean and StdErr
        df.temp<-data.frame(Mean=res.chlaEQR$period$mean,StdErr=res.chlaEQR$period$stderr)
        df.temp$Indicator<-"chlaEQR"
        res.ind.temp<-df.temp
        
        df.temp<-data.frame(Mean=res.TNsummer$period$mean,StdErr=res.TNsummer$period$stderr)
        df.temp$Indicator<-"TNsummerEQR"
        res.ind.temp<-bind_rows(res.ind.temp,df.temp)
        
        df.temp<-data.frame(Mean=res.TNwinter$period$mean,StdErr=res.TNwinter$period$stderr)
        df.temp$Indicator<-"TNwinterEQR"
        res.ind.temp<-bind_rows(res.ind.temp,df.temp)
        
        df.temp<-data.frame(Mean=res.Secchi$period$mean,StdErr=res.Secchi$period$stderr)
        df.temp$Indicator<-"Secchi"
        res.ind.temp<-bind_rows(res.ind.temp,df.temp)
        
        res.ind.temp$WB<-wblist$WB[iWB]
        res.ind.temp$Type<-wblist$typology[iWB]
        
        #Combine results of different indicators - Random values
        df.temp<-data.frame(Estimate=res.chlaEQR$indicator_sim,Code=res.chlaEQR$result_code)
        df.temp$Indicator<-"chlaEQR"
        df.temp$sim<-1:nsim
        res.rnd.temp<-df.temp
        
        df.temp<-data.frame(Estimate=res.TNsummer$indicator_sim,Code=res.TNsummer$result_code)
        df.temp$Indicator<-"TNsummerEQR"
        df.temp$sim<-1:nsim
        res.rnd.temp<-bind_rows(res.rnd.temp,df.temp)
        
        df.temp<-data.frame(Estimate=res.TNwinter$indicator_sim,Code=res.TNwinter$result_code)
        df.temp$Indicator<-"TNwinterEQR"
        df.temp$sim<-1:nsim
        res.rnd.temp<-bind_rows(res.rnd.temp,df.temp)
        
        df.temp<-data.frame(Estimate=res.Secchi$indicator_sim,Code=res.Secchi$result_code)
        df.temp$Indicator<-"Secchi"
        df.temp$sim<-1:nsim
        res.rnd.temp<-bind_rows(res.rnd.temp,df.temp)
        
        res.rnd.temp$WB<-wblist$WB[iWB]
        res.rnd.temp$Type<-wblist$typology[iWB]
        
        if(iWB==1){
          res.ind<-res.ind.temp
          res.rnd<-res.rnd.temp
        }else{
          res.ind<-bind_rows(res.ind,res.ind.temp)
          res.rnd<-bind_rows(res.rnd,res.rnd.temp)
        }
      }  
    }    
    
    # Get indicator categories based on mean values
    res.ind<- res.ind %>% select(WB,Type,Indicator,Mean,StdErr)
    res.ind<- res.ind %>% left_join(df.bounds, by=c("Indicator"="Indicator","Type"="Type"))
    
    res.ind$Value<-res.ind$Mean
    # Do we show mean concentrations where the indicator is EQR value?
    #res.ind$Value<-ifelse(res.ind$UseEQR==1,(res.ind$Mean/res.ind$Ref),res.ind$Mean)
    res.ind<-GetClass(res.ind)
    
    # Get indicator categories for MC results
    res.rnd<- res.rnd %>% select(WB,Type,Indicator,sim,Estimate,Code)
    res.rnd<- res.rnd %>% left_join(df.bounds, by=c("Indicator"="Indicator","Type"="Type"))
    names(res.rnd)[names(res.rnd)=="Estimate"]<-"Value"
    
    #res.rnd$Value<-ifelse(res.rnd$UseEQR==1,(res.rnd$Estimate/res.rnd$Ref),res.rnd$Estimate)
    res.rnd<-GetClass(res.rnd)
    cat(paste0("Sim results: ",nrow(res.rnd),"\n"))
    
    #Find counts for each category
    
    res.rnd.count <- res.rnd %>% filter(!is.na(ClassID)) %>%
      group_by(WB,Type,Indicator,ClassID) %>% summarise(n=n())
    res.rnd.count$ClassID<-paste0("C",res.rnd.count$ClassID)
    res.rnd.count$n<-res.rnd.count$n/nsim
    
    res.rnd.count<-spread(res.rnd.count, ClassID, n, fill = NA)
    res.ind<-left_join(res.ind,res.rnd.count)
    
    names(res.ind)[names(res.ind)=="C1"]<-"fBad"
    names(res.ind)[names(res.ind)=="C2"]<-"fPoor"
    names(res.ind)[names(res.ind)=="C3"]<-"fMod"
    names(res.ind)[names(res.ind)=="C4"]<-"fGood"
    names(res.ind)[names(res.ind)=="C5"]<-"fHigh"
    
    res<-list(data.frame)
    return(res.rnd)
    #Overall results
    if(FALSE){
      
      res[[1]]<-df.all %>% group_by(WB) %>%
        summarise(n=n())
      
      #QE results
      res[[2]]<-res.ind %>% group_by(WB,QualityElement) %>%
        summarise(n=n())
      
      #Indicators
      
      res[[3]]<-res.ind %>% 
        select( WB,Type,QualityElement,Indicator,Mean,StdErr,EQR,Class,fBad,fPoor,fMod,fGood,fHigh )
      #res[[3]]<-res.ind %>% select(-c(ClassID,QE)) #(WB,Type,Indicator,Mean,StdErr,UseEQR,Ref,HG,GM,MP,PB)
      
      res[[4]]<-res.rnd
      
      return(res)
    }
    
  }

#' GetCalss
#' 
#' 
#' @param df A dataframe 
#'   The dataframe should contain the
#'   following variables:
#'   
#'   \item{Value} 
#'   \item{Ref}{ } 
#'   \item{HG}{ } 
#'   \item{GM}{ } 
#'   \item{MP}{ } 
#'   \item{PB}{ } 
#'   \item{Resp}{ } 
#'   #'   
#' 
#' @examples
GetClass<-function(df){
  Categories<-c("Bad","Poor","Mod","Good","High","Ref")
  names(df)[names(df)=="RefCond"]<-"Ref"
  names(df)[names(df)=="H.G"]<-"HG"
  names(df)[names(df)=="G.M"]<-"GM"
  names(df)[names(df)=="M.P"]<-"MP"
  names(df)[names(df)=="P.B"]<-"PB"
  #names(df)[names(df)==""]<-""
  
  df$Resp<-ifelse(df$HG > df$GM,-1,1)
  df$class1<-ifelse(df$Resp==1,df$Value<df$Ref,df$Value>df$Ref)
  df$class2<-ifelse(df$Resp==1,df$Value<df$HG,df$Value>df$HG)
  df$class3<-ifelse(df$Resp==1,df$Value<df$GM,df$Value>df$GM)
  df$class4<-ifelse(df$Resp==1,df$Value<df$MP,df$Value>df$MP)
  df$class5<-ifelse(df$Resp==1,df$Value<df$PB,df$Value>df$PB)
  df$ClassID<-df$class1+df$class2+df$class3+df$class4+df$class5+1
  df$Class<-Categories[df$ClassID]
  df$Bnd1<-df$MP-2*(df$MP-df$PB)
  df$Bnd2<-df$PB
  df$Bnd1<-ifelse(df$ClassID==2,df$PB,df$Bnd1)
  df$Bnd2<-ifelse(df$ClassID==2,df$MP,df$Bnd2)
  df$Bnd1<-ifelse(df$ClassID==3,df$MP,df$Bnd1)
  df$Bnd2<-ifelse(df$ClassID==3,df$GM,df$Bnd2)
  df$Bnd1<-ifelse(df$ClassID==4,df$GM,df$Bnd1)
  df$Bnd2<-ifelse(df$ClassID==4,df$HG,df$Bnd2)
  df$Bnd1<-ifelse(df$ClassID==5,df$HG,df$Bnd1)
  df$Bnd2<-ifelse(df$ClassID==5,df$Ref,df$Bnd2)
  df$EQR<-0.2*((df$ClassID-1)+(df$Value-df$Bnd1)/(df$Bnd2-df$Bnd1))
  df$EQR<-ifelse(df$ClassID>5,1,df$EQR)
  df$ClassID<-ifelse(df$ClassID>5,5,df$ClassID)
  df<-select(df,-c(Resp,class1,class2,class3,class4,class5,Bnd1,Bnd2))
  return(df)
}

SalinityReferenceValues <- function(df.data,df.bounds,indicator,missing=1){
  typology <- substr(as.character(distinct(df.data,typology)[1,1]),4,5)
  refcond<-filter(df.bounds,Type==typology,Indicator==indicator)
  refcond<-refcond[,grep("Sali_", names(refcond), value=TRUE)]
  refcond<-as.numeric(refcond[1,])
  refcond[is.na(refcond)]<-missing
  return(refcond)
}

