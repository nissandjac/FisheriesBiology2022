getOM_MSY_deter <- function(df, tEnd = df$tEnd, dirs){
  
  
  p <- read.csv(file.path(dirs,'reference_points.csv'))
  
  
  idx <- which.min((df$M0[tEnd]-p$M)^2)
  
  
  if(df$M0[df$tEnd] > max(p$M)){
    print('recalculating reference points')
    
    print(df$M0[df$tEnd])
    
    ls <- getOM_MSY(df)
  }else{
    ls <- list(MSY = p$MSY[idx], Fmsy = p$Fmsy[idx], Bmsy = p$Bmsy[idx])
  }
  
  return(ls) 
}