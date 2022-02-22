getParameters_RAM <- function(params, ts, nstocks = length(unique(params$stockid))){
  
  ##
  #df <- load('data/RAM/processed/SurplusProdStocks.RData')
  
  # Prepare a parameter set for the first species in the data frame 
  spp <- unique(spstocks_params$stockid)
  
  
  compile("src/recruitment.cpp") # Load the TMB model 
  dyn.load(dynlib("src/recruitment"))
  source('R/calcSRparams.R')
  source('R/load_data_seasons.R')
  
  df.rec <- data.frame(nms = unique(spstocks_params$stockid),
                       h = NA,
                       SSB0 = NA,
                       R0 = NA,
                       SDR = NA)
  id <- unique(spstocks_params$stockid)
  
  ls.df <- list()
  
  ls.catch <- data.frame(stockid = NA,
                         year = NA, 
                         Catch = NA)
  
  
  for(i in 1:nrow(df.rec)){
    
    
    
    ts <- filter(spstocks_ts, stockid == id[i])
    params <- filter(spstocks_params, stockid == id[i])
    
    
    
    # Estimate SR params from recruitment 
    SR <- calcSRparams(ts)
    genus <- strsplit(params$sciname[1], split = ' ')[[1]][1]
    species <- strsplit(params$sciname[1], split = ' ')[[1]][2]
    
    if(id[i] == 'SARDWSE'){
      
      params$biovalue[params$bioid ==  'M-1/yr'] <- 0.6
      
    }
    
    
    if(species == 'pallasii'){ # Fishbase lookup error 
      species <- 'pallasii pallasii'
    }
    
    if(species == 'melanostictus'){
      species <- 'sagax'
    }
    
    if(species == 'teres'){
      species <- 'sadina' 
    }
    
    if(ts$stocklong[1] == 'Sardinella West Africa'){
      genus <- 'Sardinops'
      species <- 'sagax'
    }
    
    if(ts$stocklong[1] == "Sand eel Sandeel Area 1" | ts$stocklong[1] == "Sand eel Sandeel Area 2" 
       | ts$stocklong[1] == "Sand eel Sandeel Area 3" | ts$stocklong[1] == "Sand eel Sandeel Area 4"){
      
      genus <- 'Ammodytes'
      species <- 'marinus'
      
    }
    
    if(genus == 'Strangomera'){
      genus <- 'Clupea'
    }
    
    FL <- Plot_taxa( Search_species(Genus=genus,Species=species)$match_taxonomy ,plot_lines = FALSE, verbose = FALSE)
    
    df.rec[i,2:5] <- as.numeric(SR[1:4])
    
    
    # Maturity params 
    isAmat <- is.element('A50', params$bioshort)
    
    if(isAmat){
      
      if(nrow(Amat)>1){
        Amatin <- mean(as.numeric(Amat$biovalue))
      } 
    }else{
      Amatin <- as.numeric(exp(FL[[1]]$Mean_pred['tm']))
    }
    #  
    #     Amat <- params[grep('A50', x = params$bioshort),]
    # 
    # if(length(grep('sel', Amat$stocklong))>0){
    #   stop()
    # }
    # 
    
  
    
    # Selectivity params 
    sel <- filter(params, grepl('sel', biolong))
    tau_F <- NA
    if(nrow(sel) == 0){
      selF <- filter(params, bioid ==  "F-AGE-yr-yr")
      
      if(nrow(selF)>0){
        tau_F <- as.numeric(strsplit(selF$biovalue, split = '+')[[1]][1])*0.5
        
      }  
    }
    
    if(is.na(tau_F)){
      tau_F <- Amatin
    }
    # natural mortality 
    
    isM <- is.element('M-1/yr', params$bioid)
    
    if(isM){
    Mdb <- filter(params, bioid == "M-1/yr")
    
    if(nrow(Mdb)>0){
    Min <- mean(as.numeric(Mdb$biovalue), na.rm =TRUE)
    
    }
    
    if(Mdb$biounitsshort[1] != '1/yr'){
      Mdb <- filter(params, bioid == "M-1/month")
      Min <- Mdb$biovalue*12 
    }
      
    }else{
      Min <- exp(FL[[1]]$Mean_pred["M"])
    }
     
   # Get Linf
    
    isLinf <- is.element('Linf', params$bioshort)
    
    if(isLinf){
    Linf <- filter(params, bioshort == 'Linf')
    Linf <- as.numeric(Linf$biovalue)
    
    }else{
      
    isWinf <- is.element('Winf', params$bioshort)
      if(isWinf){
      Winf <- filter(params, bioshort == 'Linf')
      Linf <- (Winf$biovalue/0.01)^(1/3)
      
      }
    else{
      Linf <- exp(FL[[1]]$Mean_pred['Loo'])
      }
    }
    
    # Get Von B growth rate 
    isK <- is.element('VB-k', params$bioshort)
    
    if(!isK){
      Kin <- exp(FL[[1]]$Mean_pred['K'])
    }else{
      Kin <- as.numeric(params$biovalue[grep('VB-k',params$bioshort)])
    }
    
    
    ist0 <- is.element('VB-t0', params$bioshort)
    
    if(ist0){
      t0in <- as.numeric(params$biovalue[grep('VB-t0',params$bioshort)]) 
    }else{
      t0in <- -1 # Average of all forage fish
    }
    
    
    # max age 
    isMmax <- is.element('Maximum age', params$biolong)
    
    if(isMmax){
      Mmaxin <- as.numeric(params$biovalue[which(params$biolong == 'Maximum age')]) 
    
      if(is.na(Mmaxin)){
        Mmaxin <- as.numeric(strsplit(params$biovalue[which(params$biolong == 'Maximum age')], split = '\\+')[[1]][1])
        
        
      }
      }else{
      Mmaxin <- exp(FL[[1]]$Mean_pred['tmax'])
    }
    
    
    
    # Recalculate exploitation rate to fishing mortality 
    if(is.na(df.rec$SDR[i])){
      Rin <- as.numeric(FL[[1]]$Mean_pred['r'])
    }else{
      Rin <- df.rec$SDR[i]
    }
    
    Fin <- ts$TCbest/ts$TBbest
    
    if(is.nan(mean(Fin, na.rm = TRUE))){
      Fin <- Min/2
      }
    
   
     if(is.nan(mean(Fin, na.rm = TRUE))){
      Fstart <- Min/4
    
    }else{
      Fstart <- Fin[1]
    }
    
    SDF <- sd(Fin)*2
    
    if(is.na(SDF)){
      SDF <- 0.1
    }
    
    if(is.na(df.rec$h[i])){
      hin <- as.numeric(FL[[1]]$Mean_pred['h'])
    }else{
      hin <- df.rec$h[i]
    }
    
    if(is.na(Fstart)){
      Fstart <- 0
    }
    
    df <- load_data_seasons(nseason = 1, # Number of seasons in OM
                            nspace = 1,  # Number of spatial areas (only works with two right now)
                            SDR = Rin,  # Recruitment deviations 
                            SDM = Rin*.5, 
                            myear = 50,
                            h = hin,
                            Linf = Linf,
                            tau_sel = tau_F,
                            tau = Amatin,
                            smul = 0,
                            maxage = Mmaxin,
                            K = Kin,
                            t0 = t0in,
                            catch.true = FALSE,  # input catch (TRUE) or fishing mortality (FALSE)
                            SDsurv = 0.2, # Survey uncertainty 
                            bfuture = 1,
                            b = 0,
                            Fpast = Fstart,
                            M = Min,
                            fishing.type = 'AR',
                            F0 = mean(Fin, na.rm = TRUE),
                            SDF = SDF
                            )# AR 1 natural mortality   
    
    # Estimate R0 

    ls.df[[i]] <- df
    
    # Save the catch data.frame 
    Catchin <- ts$TCbest
    
    if(i == 1){
    Catchdf <- data.frame(
      stock = unique(ts$stockid),
      year = ts$year,
      Catch = ts$TCbest)
    }else{
      Catchdf <- rbind(Catchdf,
                       data.frame(
                         stock = unique(ts$stockid),
                         year = ts$year,
                         Catch = ts$TCbest))
      
    }
    
    
    
  }
  
  # Sort by highest average catch
  Catchtot <- Catchdf[Catchdf$year > 2010,] %>%
    group_by(stock) %>%
    summarise(Catch = mean(Catch, na.rm = TRUE)) %>%
    arrange(-Catch)
    
  # Decide which species to sort
  
  spp.save <- Catchtot[1:nstocks,]
  idxnames <- spp.save$stock[1:nstocks]
  # Add names to list
  names(ls.df) <- spp
  
  ls.df <- ls.df[which(df.rec$nms %in% idxnames)]
  
  
 return(ls.df)
  }
  