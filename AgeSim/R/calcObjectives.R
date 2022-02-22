calcObjectives <- function(MSE.list, 
                           sim.year = 50,
                           projection.year = 99){
  
  

    # 
    # short term catch 
    catch.short <- MSE[MSE$time %in% (sim.year:(sim.year+10)) & MSE$model == 'OM',] %>% 
      group_by(HCR,M, scenario) %>% 
      summarise(catchM = median(catch, na.rm =TRUE),
                catch90 = quantile(catch, probs = 0.9, na.rm =TRUE),
                catch10 = quantile(catch, probs = 0.1, na.rm =TRUE),
                catch75 = quantile(catch, probs = 0.75, na.rm =TRUE),
                catch25 = quantile(catch, probs = 0.25, na.rm =TRUE)
      ) %>% 
      mutate(model = paste(HCR,M, sep = '-'))
    
    # Long term catch 
    catch.long <- MSE[MSE$time %in% ((projection.year-10):projection.year) & MSE$model == 'OM',] %>% 
      group_by(HCR,M) %>% 
      summarise(catchM = median(catch, na.rm =TRUE),
                catch90 = quantile(catch, probs = 0.9, na.rm =TRUE),
                catch10 = quantile(catch, probs = 0.1, na.rm =TRUE),
                catch75 = quantile(catch, probs = 0.75, na.rm =TRUE),
                catch25 = quantile(catch, probs = 0.25, na.rm =TRUE)
      ) %>% 
      mutate(model = paste(HCR,M, sep = '-'))
    
    # Average catch
    catch.med <- MSE[(MSE$time >sim.year) & MSE$model == 'OM',] %>% 
      group_by(HCR,M) %>% 
      summarise(catchM = median(catch, na.rm =TRUE),
                catch90 = quantile(catch, probs = 0.9, na.rm =TRUE),
                catch10 = quantile(catch, probs = 0.1, na.rm =TRUE),
                catch75 = quantile(catch, probs = 0.75, na.rm =TRUE),
                catch25 = quantile(catch, probs = 0.25, na.rm =TRUE)
      ) %>% 
        mutate(model = paste(HCR,M, sep = '-'))
    
    # Long term B/B0
    SSB.long <- MSE[MSE$time %in% ((projection.year-10):projection.year) & MSE$model == 'OM',] %>% 
      group_by(HCR,M) %>% 
      summarise(SSBM = median(biomass/S0, na.rm =TRUE),
                SSB90 = quantile(biomass/S0, probs = 0.9, na.rm =TRUE),
                SSB10 = quantile(biomass/S0, probs = 0.1, na.rm =TRUE),
                SSB75 = quantile(biomass/S0, probs = 0.75, na.rm =TRUE),
                SSB25 = quantile(biomass/S0, probs = 0.25, na.rm =TRUE)
      ) %>% 
      mutate(model = paste(HCR,M, sep = '-'))
    
    # short term B/B0
    SSB.short <- MSE[MSE$time %in% (sim.year:(sim.year+10)) & MSE$model == 'OM',] %>% 
      group_by(HCR,M) %>% 
      summarise(SSBM = median(biomass/S0, na.rm =TRUE),
                SSB90 = quantile(biomass/S0, probs = 0.9, na.rm =TRUE),
                SSB10 = quantile(biomass/S0, probs = 0.1, na.rm =TRUE),
                SSB75 = quantile(biomass/S0, probs = 0.75, na.rm =TRUE),
                SSB25 = quantile(biomass/S0, probs = 0.25, na.rm =TRUE)
      ) %>% 
      mutate(model = paste(HCR,M, sep = '-'))
    
    # Fraction of years overfished 
    SSB.long <- MSE[MSE$time %in% (sim.year:(sim.year+10)) & MSE$model == 'OM',] %>% 
      group_by(HCR,M) %>% 
      summarise(SSBM = median(biomass/S0, na.rm =TRUE),
                SSB90 = quantile(biomass/S0, probs = 0.9, na.rm =TRUE),
                SSB10 = quantile(biomass/S0, probs = 0.1, na.rm =TRUE),
                SSB75 = quantile(biomass/S0, probs = 0.75, na.rm =TRUE),
                SSB25 = quantile(biomass/S0, probs = 0.25, na.rm =TRUE)
      ) %>% 
      mutate(model = paste(HCR,M, sep = '-'))
    
    # Fraction of years with impaired recruitment 
    
    # Interannual catch variability 
    
    # Total catch variability 
    
    
    
    
    
  
}