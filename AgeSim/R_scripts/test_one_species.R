source('R/load_files.R')

# Set up parameters with life histories 

df <- load_data_seasons(nseason = 1, # Number of seasons 
                        myear = 50, # Number of years to run
                        Linf = 30,  # Linfinity
                        K = 1,  # Von Bertalanffy K
                        t0 = -0.1, # VonB t0
                        gamma = 1.1, # Slope of maturity function
                        tau = 3, # age at 50% maturity
                        gamma_sel = 1.1, # Slope of selectivity function (trawl)
                        tau_sel = 3, # Age at 50% selectivity 
                        maxage = 10, # 
                        M = 0.2, # Natural mortality (can be a vector of size myear)
                        b = 1, # bias adjustment 
                        R0 = 1e4, # Unfished recruitment 
                        h = 0.8, # Steepness 
                        q = .7, # survey catchability 
                        SDsurv = 0.1, # Survey observation error 
                        SDcatch = 0.01, # Catch observation error 
                        nsurvey = 1,  # Number of surveys 
                        SDR = .6,  # Recruitment deviations 
                        catch.true = FALSE, # Use fishing mortality or actual catch
                        F0 = 0.2,  # Fishing mortality (can be a vector)
                        fishing.type = 'AR', # Fishing type, AR or constant
                        SDF = 0.05, # Fishing mortality deviations 
                        mortality = 'constant', # Mortality type. AR, constant
                        seed = NA # Set seed for replications 
    

)


tmp <- run.agebased.true.catch(df)

### Plot some stuff from the simulation ### 

df.plot <- data.frame(time = 1:df$tEnd, 
                      SS0 = as.numeric(tmp$SSB/tmp$SSB_0),
                      recruits = tmp$R.save,
                      Catch = tmp$Catch,
                      survey = tmp$survey
                      ) %>% pivot_longer(2:5) %>% 
  ggplot(aes(x = time,  y = value, color = name)) + geom_line(size = 1.2)+
  facet_wrap(~name, scales = 'free_y')+theme_classic()+scale_y_continuous('')+theme(legend.position = 'top',
                                                                                    legend.title = element_blank())

df.plot
