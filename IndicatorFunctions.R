# Aggregation principle used for Chlorophyll
Aggregate_year_station <- function(df) {
  yearmeans <- df %>%    group_by(year,station) %>%
                         summarise(xvar = mean(xvar)) %>%
                         group_by(year) %>%
                         summarise(xvar = mean(xvar))
                         
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans)
  return(res)
}

# Aggregation principle used for nutrients
Aggregate_year <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvar))
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans)
  return(res)
}

# Aggregation principle used for Secchi depth
Aggregate_period <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvar))
  
  periodmean <- mean(df$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans)
  return(res)
}

# Aggregation principle used for Chlorophyll
AggregateEQR_year_station <- function(df) {
  
  df <- mutate(df,xvarEQR = ifelse(xvar<RefCond,1,RefCond/xvar))

  yearmeans <- df %>%    group_by(year,station) %>%
    summarise(xvarEQR = mean(xvarEQR)) %>%
    group_by(year) %>%
    summarise(xvar = mean(xvarEQR))   # should be returned in xvar
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans)
  return(res)  
}

# Aggregation principle used for nutrients
AggregateEQR_year <- function(df) {
  
  df <- mutate(df,xvarEQR = RefCond/xvar)
  
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvarEQR))   # should be returned in xvar
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans)
  return(res)  
}

# Indicator response negative to degradation, i.e. Secchi depth
AggregateEQR_N_period <- function(df) {
  
  df <- mutate(df,xvarEQR = xvar/RefCond)
  
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvarEQR))   # should be returned in xvar
  
  periodmean <- mean(df$xvarEQR)
  res <- list(periodmean=periodmean,yearmeans=yearmeans)
  return(res)  
}

# Calculation of BQI indicator according to Handbook
BQIbootstrap <- function(df) {
  yearstatmeans <- df %>%    group_by(year,station) %>%
    summarise(xvar = mean(xvar))
  Nyearstat <- yearstatmeans %>% group_by(year) %>% summarise(n_station = length(xvar))
  BQIsimyear <- mat.or.vec(length(Nyearstat$n_station), 1)
  for(i in 1:length(Nyearstat$n_station)) {
     BQIsim <- trunc(runif(9999,1,Nyearstat$n_station[i]+1))
     BQIsim <- yearstatmeans$xvar[yearstatmeans$year == Nyearstat$year[i]][BQIsim]
     BQIsimyear[i] <- quantile(BQIsim,probs=0.2)
  }

  periodmean <- mean(BQIsimyear)
  yearmeans <- data.frame(year=Nyearstat$year,xvar = BQIsimyear)
  res <- list(periodmean=periodmean,yearmeans=yearmeans)
  return(res)  
}

#' Generic routine for calculating indicator statistics
#' 
#' @param df A dataframe with monitoring data from the Swedish Monitoring program. 
#'   The dataframe should contain the
#'   following variables:
#'   
#'   \describe{ 
#'   \item{station}{An identifier for the monitoring station.} 
#'   \item{date}{Date of the observation.} 
#'   \item{institution}{Provider of the observation.} 
#'   \item{chla}{Chlorophyll a concentration in sample.} 
#'   
#' @param MonthInclude A list of months to be included in the indicator
#' @param var_list List of variance components
#' @param n_iter Number of iterations for Monte Carlo simulation
#'   
#' @return
#' @export
#' 
#' @examples 

CalculateIndicator <-
  function(Indicator,df,RefCond_sali,var_list,startyear,endyear,n_iter=1000,confidence_lvl=0.95) {
    # Set flag to zero and change it for error handling below
    flag <- 0
# Select the observation variable for the indicator
    xvar <- switch(Indicator,
                   Chla         = df$chla,
                   ChlaEQR      = df$chla,
                   Secchi       = df$secchi,
                   SecchiEQR    = df$secchi,
                   DINsummer    = df$DIN,
                   DINsummerEQR = df$DIN,
                   DIPsummer    = df$DIP,
                   DIPsummerEQR = df$DIP,
                   TNsummer     = df$TN,
                   TNsummerEQR  = df$TN,
                   TNwinter     = df$TN,
                   TNwinterEQR  = df$TN,
                   TPsummer     = df$TP,
                   TPsummerEQR  = df$TP,
                   TPwinter     = df$TP,
                   TPwinterEQR  = df$TP,
                   BQI          = df$BQI)
    df <- mutate(df,xvar=xvar)
# Associating indicators with transformation from observations
    f_fun <- switch(Indicator,
                    Chla         = Aggregate_year_station,
                    ChlaEQR      = AggregateEQR_year_station,
                    Secchi       = Aggregate_period,
                    SecchiEQR    = AggregateEQR_N_period,
                    DINsummer    = Aggregate_year,
                    DINsummerEQR = AggregateEQR_year,
                    DIPsummer    = Aggregate_year,
                    DIPsummerEQR = AggregateEQR_year,
                    TNsummer     = Aggregate_year,
                    TNsummerEQR  = AggregateEQR_year,
                    TNwinter     = Aggregate_year,
                    TNwinterEQR  = AggregateEQR_year,
                    TPsummer     = Aggregate_year,
                    TPsummerEQR  = AggregateEQR_year,
                    TPwinter     = Aggregate_year,
                    TPwinterEQR  = AggregateEQR_year,
                    BQI          = BQIbootstrap)
# Assigning transformations for measurements to obtain normal distributed variates
    g_fun <- switch(Indicator,
                    Chla         = log,
                    ChlaEQR      = log,
                    Secchi       = log,
                    SecchiEQR    = log,
                    DINsummer    = log,
                    DINsummerEQR = log,
                    DIPsummer    = log,
                    DIPsummerEQR = log,
                    TNsummer     = log,
                    TNsummerEQR  = log,
                    TNwinter     = log,
                    TNwinterEQR  = log,
                    TPsummer     = log,
                    TPsummerEQR  = log,
                    TPwinter     = log,
                    TPwinterEQR  = log,
                    BQI          = identity)    
# Assigning inverse transformations of g_fun
    g_fun_inv <- switch(Indicator,
                    Chla         = exp,
                    ChlaEQR      = exp,
                    Secchi       = exp,
                    SecchiEQR    = exp,
                    DINsummer    = exp,
                    DINsummerEQR = exp,
                    DIPsummer    = exp,
                    DIPsummerEQR = exp,
                    TNsummer     = exp,
                    TNsummerEQR  = exp,
                    TNwinter     = exp,
                    TNwinterEQR  = exp,
                    TPsummer     = exp,
                    TPsummerEQR  = exp,
                    TPwinter     = exp,
                    TPwinterEQR  = exp,
                    BQI          = identity) 
# Select the months for including in indicator calculation
    MonthInclude <- switch(Indicator,
                           Chla         = c(6,7,8),
                           ChlaEQR      = c(6,7,8),
                           Secchi       = c(6,7,8),
                           SecchiEQR    = c(6,7,8),
                           DINsummer    = c(6,7,8),
                           DINsummerEQR = c(6,7,8),
                           DIPsummer    = c(6,7,8),
                           DIPsummerEQR = c(6,7,8),
                           TNsummer     = c(6,7,8),
                           TNsummerEQR  = c(6,7,8),
                           TNwinter     = c(11,12,1,2),
                           TNwinterEQR  = c(11,12,1,2),
                           TPsummer     = c(6,7,8),
                           TPsummerEQR  = c(6,7,8),
                           TPwinter     = c(11,12,1,2),
                           TPwinterEQR  = c(11,12,1,2),
                           BQI          = c(1,2,3,4,5,6,7,8,9,10,11,12))
# Switch year for winter months (Nov+Dec) to in
    if (Indicator %in% c("TNwinterEQR","TPwinterEQR")) {
      df <- mutate(df,year=ifelse(month %in% c(11,12),year+1,year))
    }
# Filter dataframe to include observations used in indicator only
    df <- Filter_df(df,MonthInclude,startyear,endyear)    
# setting RefCond depending on salinity for indicators with salinity correction
    RefCond <- mat.or.vec(nrow(df), 1)
    if (Indicator %in% c("ChlaEQR","SecchiEQR","DINsummerEQR","DIPsummerEQR","TNsummerEQR","TPsummerEQR","TNwinterEQR","TPwinterEQR")) {
       df <- filter(df,!is.na(sali))
       RefCond <- mat.or.vec(nrow(df), 1)
       sali_class <- findInterval(df$sali, c(seq(0, 35)))
       for (i in 1:nrow(df)) {RefCond[i] <- RefCond_sali[sali_class[i]]}
       }
    df <- mutate(df,RefCond = RefCond) 
# Calculate number of years, stations, months, institutions and combinations thereof in df 
    ndf <- DF_Ncalculation(df)
    # Return from function if no observations for calculation
    if (ndf$n_obs == 0) return(list(result_code=-90))
# Estimate mean of the log-transformed chla obs
    alpha <- mean(g_fun(df$xvar))
    alpha <- df %>% group_by(year) %>% summarise(mean = mean(g_fun(xvar)))
# Calculate indicator
    mu_indicator <- f_fun(df)
# Simulate system with random variables for estimating the variance of the indicator
    simres <- vector("numeric",n_iter)
    simresyear <- matrix(nrow=n_iter,ncol=ndf$n_year)
    simresyear <- matrix(nrow=ndf$n_year,ncol=n_iter)
# simulation loop - simres contains the residuals from n_iter simulations
    for (isim in 1:n_iter) {
      # simulate variations in the random factors using the data structure
      simulobs <- SetVector_IndicatorSim(alpha$mean,ndf,var_list,df,length(MonthInclude))
      # backtransform simulations from log domain to original domain
      simulobs <- g_fun_inv(simulobs)
      # add simulated observation to df
      simul_df <- data.frame(year=df$year,station=df$station,xvar=simulobs,RefCond=df$RefCond)
      # Calculate indicator value for each year and period
      simul_indicator <- f_fun(simul_df)
      simresyear[,isim]=simul_indicator$yearmeans$xvar
      simres[isim] <- simul_indicator$periodmean
    } # end simulation loop
    
# Adjust simulations to have zero mean and then add indicator means
    simres <- g_fun_inv(g_fun(simres)-g_fun(mean(simres))+g_fun(mu_indicator$periodmean))
    simresyear <- g_fun_inv(g_fun(simresyear)-g_fun(apply(simresyear,1,mean))+g_fun(mu_indicator$yearmeans$xvar))
    
# Calculate statistics
    period <- data.frame(mean=mean(simres),stderr=sd(simres),lower = quantile(simres,probs=c((1-confidence_lvl)/2)),upper = quantile(simres,probs=c(1-(1-confidence_lvl)/2)),row.names = NULL)
    annual <- data.frame(year = mu_indicator$yearmeans$year,
                         mean = apply(simresyear,1,mean),
                         stderr = apply(simresyear,1,sd),
                         lower = apply(simresyear,1,quantile,probs=c((1-confidence_lvl)/2)),
                         upper = apply(simresyear,1,quantile,probs=c(1-(1-confidence_lvl)/2)))
    flag <- ifelse(length(annual$mean)<3,-1,flag)
    obs_sim <- as.numeric(simres)
    
    res <- list(period=period,annual=annual,indicator_sim=simres,n_list=ndf,result_code=flag)
    return(res)
  }