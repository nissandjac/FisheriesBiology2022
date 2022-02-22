getOM_MSY <- function(df.msy, plotfigs = FALSE, yr = df.msy$tEnd){
  
 
  
  Ftest <- seq(0.01,max(df.msy$M0)*4, length.out = 100)
  Catches <- rep(NA, length(Ftest))
  Biomass <- rep(NA, length(Ftest))
  
  # Make model deterministic 
  df.msy$b <-0*df.msy$b
  df.msy$M0 <- df.msy$M0[yr]
  df.msy$parms$Rin <- df.msy$parms$Rin*0
  df.msy$catch.true <- FALSE
  
  # Stupid hack for a stupid programmer 
  df.msy$F0 <- rep(1,df.msy$nyear)

    # 
  # if(df.msy$nyear>30){
  # df.msy$nyear <- 30 # Just run 30 years
  # 
  
  for(i in 1:length(Ftest)){
    
    df.msy$F0 <- df.msy$F0*0+Ftest[i]
    
    sim.data.tmp <- run.agebased.true.catch(df.msy)
    Catches[i] <- sim.data.tmp$Catch[df.msy$nyear]
    Biomass[i] <- sim.data.tmp$V.save[,1,1][df.msy$nyear]
    
    
  }
  
  MSY <- max(Catches)
  Fmsy <- Ftest[which.max(Catches)]
  Emsy <- max(Catches)/Biomass[which.max(Catches)]
  Bmsy <- Biomass[which.max(Catches)]
  
  if(plotfigs == TRUE){
  plot(Ftest,Catches, type = 'l')
  lines(Ftest, rep(MSY, length(Ftest)))
  }
  
  #plot(Ftest,Catches)  
  
return(list(MSY = MSY, Fmsy = Emsy, Bmsy = Bmsy))
  
}