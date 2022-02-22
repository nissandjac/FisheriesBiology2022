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
                        SDR = .0,  # Recruitment deviations 
                        catch.true = FALSE, # Use fishing mortality or actual catch
                        F0 = c(1,0),  # Fishing mortality in the two areas 
                        mortality = 'constant', # Mortality type. AR, constant
                        # Movement parameters 
                        movemaxinit = c(1, .2), # Fraction of fish moving into the other cell
                        R0_dist = c(0.1, .5) # A number between zero and 1 describing the fraction of R0 in a given area 
                        
                        
)




tmp <- run.agebased.true.catch(df)

names(tmp)

df.plot <- data.frame(time = rep(1:df$tEnd,2), 
                      SSB = c(tmp$SSB[,1], tmp$SSB[,2]),
                      R = c(tmp$R.save[,1], tmp$R.save[,2]),
                      area = rep(c('fished','mpa'), each = df$tEnd)) %>% 
  ggplot(aes(x = time, y = SSB, color = area))+geom_line()+theme_classic()

df.plot

