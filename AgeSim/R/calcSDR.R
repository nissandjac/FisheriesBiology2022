calcSRparams <- function(SSB, R, params){
  
  SSB <- ts$SSB

  if(all(is.na(SSB))){
    SSB <- ts$TBbest
  }

  R <- ts$R

  #  Fit a standard BH function, estimate R0, h, and SDR 
  
  #  Prepare data 
  parms <- list(
      logh = log(0.7),
      logSSBzero = log(max(SSB, na.rm = TRUE)),
      logRzero = log(0.5*max(R, na.rm = TRUE)),
      logSDR = log(0.7)
      
    )
  

  
  df <- data.frame(
    SSB = SSB,
    R = R
  )
  df <- na.omit(df)
  
  if(nrow(df) < 18){
    pred <- NA
    warning('Not enough years to estimate SR')
  }else{
  # Remove NAs 
  df.tmb <- list(SSBn = df$SSB, Robs = df$R, nobs = nrow(df))
  # Phase 1 
  obj <-MakeADFun(df.tmb,parms,DLL="recruitment")
  tt <- obj$report()
  
  lower <- obj$par-Inf
  lower[names(lower) == "logh"]<- log(0.1)
  # lower[names(lower) == "logSSBzero"]<- log(min(SSB, na.rm = TRUE))
  # 
  upper <- obj$par+Inf
  # upper[names(upper) == "logSSBzero"]<- log(max(SSB, na.rm = TRUE))
  upper[names(upper) == "logh"]<- log(.99)
  # 
  # 
  opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper, 
              control = list(iter.max = 4000,
                             eval.max = 4000)) #
  
  
  # Plot the predicted SSR 
  pred <- exp(opt$par)
  
  SSBpred <- seq(0, max(SSB, na.rm = TRUE), length.out = 100)
  
  Rpred <- 4.0*pred[1]*pred[3]*SSBpred/(pred[2]*(1.0-pred[1])+ SSBpred*(5.0*pred[1]-1.0))
  
  plot(SSB,R)
  lines(SSBpred,Rpred)
  }
  
  return(pred)
}
