# Calculate the operating model reference points for a range of different mortalities # 

calcMSY <- function(df, dirs){
  
  tsim <- 30
  saving <- FALSE
  # Run a model 
  df$M0[length(df$M0)] <- df$M0[1] # For consistency in reference point calcs

  # Run operating model 
  
  sim.data <- run.agebased.true.catch(df)
  
  
  Mtest <- seq(min(df$M0)*0.1, max(df$M0)*5, length.out = 300)
  df.save <- data.frame(M = Mtest, MSY = NA, Fmsy = NA, Bmsy = NA)
  
  for(i in 1:length(Mtest)){
    df$M0[df$tEnd] <- Mtest[i]
    
    reftemp <- getOM_MSY(df,plot = FALSE,yr = df$tEnd)
    
    df.save$MSY[i] <- reftemp$MSY
    df.save$Fmsy[i] <- reftemp$Fmsy
    df.save$Bmsy[i] <- reftemp$Bmsy
    
    
  }
  
  M2 <- ggplot(melt(df.save, id.vars = 'M', measure.vars = 2:4), aes(x = M, y = value))+
    facet_wrap(~variable, scales = 'free')+
    geom_line()+
    theme_classic()
  
  print(M2)
  
  write.csv(df.save,file = paste(dirs,'reference_points.csv', sep = '/'), row.names = FALSE)
  
}

