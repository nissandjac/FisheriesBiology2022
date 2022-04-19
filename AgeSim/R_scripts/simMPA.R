#### Run an MPA example #### 
source('R/load_files.R')
library(tidyverse)

# Set up parameters with life histories 

df <- load_data_seasons(nseason = 1, # Number of seasons 
                        nspace = 2,
                        myear = 50, # Number of years to run
                        Linf = 30,  # Linfinity
                        K = 1,  # Von Bertalanffy K
                        tau = 3, # age at 50% maturity
                        tau_sel = 3, # Age at 50% selectivity 
                        maxage = 10, # 
                        M = 0.2, # Natural mortality (can be a vector of size myear)
                        b = 1, # bias adjustment 
                        R0 = 1e4, # Unfished recruitment 
                        h = 0.8, # Steepness 
                        nsurvey = 1,  # Number of surveys 
                        SDR = 1,  # Recruitment deviations 
                        catch.true = FALSE, # Use fishing mortality or actual catch
                        F0 = c(.5,0),  # Fishing mortality in the two areas 
                        mortality = 'constant', # Mortality type. AR, constant
                        # Movement parameters 
                        movemaxinit = c(.5, .2), # Fraction of fish moving into the other cell
                        R0_dist = c(0.1, 1) # A number between zero and 1 describing the fraction of R0 in a given area 
                        
                        
)




tmp <- run.agebased.true.catch(df)

names(tmp)

df.plot <- data.frame(time = rep(1:df$tEnd,2), 
                      SSB = c(tmp$SSB[,1], tmp$SSB[,2]),
                      R = c(tmp$R.save[,1], tmp$R.save[,2]),
                      area = rep(c('fished','mpa'), each = df$tEnd)) %>% 
  ggplot(aes(x = time, y = SSB, color = area))+geom_line()+theme_classic()+scale_y_log10()

df.plot


# Get the maximum sustainable yield and fishing mortality as a function of MPA size 


protec <- seq( 0  , 1 , length.out = 11) # Fraction of area protected 
F0 <- seq(0,1.5, length.out = 30) # Fishinng mortality outside MPA 


for(i in 1:length(protec)){
  for(j in 1:length(F0)){
    
    df <- load_data_seasons(nseason = 1, # Number of seasons 
                            nspace = 2,
                            myear = 50, # Number of years to run
                            Linf = 30,  # Linfinity
                            K = 1,  # Von Bertalanffy K
                            tau = 3, # age at 50% maturity
                            tau_sel = 3, # Age at 50% selectivity 
                            maxage = 10, # 
                            M = 0.2, # Natural mortality (can be a vector of size myear)
                            b = 1, # bias adjustment 
                            R0 = 1e4, # Unfished recruitment 
                            h = 0.3, # Steepness 
                            nsurvey = 1,  # Number of surveys 
                            SDR = 0,  # Recruitment deviations 
                            catch.true = FALSE, # Use fishing mortality or actual catch
                            F0 = c(F0[j],0),  # Fishing mortality in the two areas 
                            mortality = 'constant', # Mortality type. AR, constant
                            # Movement parameters 
                            movemaxinit = c(.5, .2), # Fraction of fish moving into the other cell
                            R0_dist = c(1-protec[i], protec[i]) # A number between zero and 1 describing the fraction of R0 in a given area 
                            
                            
    )
    
    tmp <- run.agebased.true.catch(df)
    
    
    
    if(i == 1 & j == 1){
      df.out <- data.frame(time = rep(1:df$tEnd,2),
                           SSB = c(tmp$SSB[,1],tmp$SSB[,2]),
                           Catch = c(tmp$Catch, rep(0, df$tEnd)),
                           R = c(tmp$R.save[,1], tmp$R.save[,2]),
                           model = rep(c('fished area', 'mpa') ,each = df$tEnd),
                           protec = rep(protec[i], df$tEnd*2),
                           F0 = rep(F0[j], df$tEnd * 2)
      )
    }else{
      
      df.out <- rbind(data.frame(time = rep(1:df$tEnd,2),
                              SSB = c(tmp$SSB[,1],tmp$SSB[,2]),
                              Catch = c(tmp$Catch, rep(0, df$tEnd)),
                              R = c(tmp$R.save[,1], tmp$R.save[,2]),
                              model = rep(c('fished area', 'mpa') ,each = df$tEnd),
                              protec = rep(protec[i], df$tEnd*2),
                              F0 = rep(F0[j], df$tEnd * 2)
                              ),
                   df.out)
                              
        
      
      
    }
  }
}

p1 <- ggplot(df.out[df.out$time == 50 & df.out$model == 'fished area',], 
       aes(x = F0, y = Catch, color = as.factor(protec)))+geom_line()+theme_classic()
  p1


p1 








#### Se Spawning biomass per areas as a function of protec ### #
  
p2 <- ggplot(df.out[df.out$time == 50,], aes(x = F0, y = SSB,linetype = model ,color = as.factor(protec)))+
  geom_line()+theme_classic()+facet_wrap(~protec)

p2  



p3 <- ggplot(df.out[df.out$time == 50,], aes(x = F0, y = R,linetype = model ,color = as.factor(protec)))+
  geom_line()+theme_classic()+facet_wrap(~protec)

p3  


p1

p2

p3


head(df.out)

  
ggplot(df.out[df.out$model != 'mpa',], aes(x = time, y = Catch, color = as.factor(F0)))+geom_line()+
  facet_wrap(~protec)+scale_y_log10()








