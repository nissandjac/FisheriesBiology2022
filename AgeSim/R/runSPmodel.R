runSPmodel <- function(r = 1.1, # Run a surplus production model 
                       K = 1000, 
                       q = 1, 
                       nyear = 100,
                       SDrun = 0.05,
                       SDsurv = 0.1,
                       F0 = 0.3){

#'@r growth rate 
#'@K unfished abundance
#'@q survey catchability 
#'@nyear number of years to run simulation 
#'@SDrun stochastic variation 
#'@SDsurv Survey uncertainty 
#'@F0 Fishing mortality 
  
# 
# r <- 1.1
# K <- 100
# q <- 0.001
# SDrun = 0.05
# nyear = 100
# F0 <- 0.3
# SDsurv <- 0.1

  if(length(F0) == 1){  
    Fin <- rep(F0, nyear)
  }
  
  if(length(Fin) != nyear){
    stop('Fishing mortality input has wrong size')
  }
  
B <- S <- C <- rep(NA, nyear) # Biomass 

B[1] <- (K-K*Fin[1])*exp(rnorm(1,mean = 0, sd = SDrun))
C[1] <- Fin[1]*B[1]
for(i in 1:(nyear-1)){
  C[i]<- Fin[i]*B[i] 
  B[i+1] <- (B[i]+B[i]*r*(1-B[i]/K)-C[i])*exp(rnorm(1, mean = 0, sd = SDrun))
#  S[i] <- B[i]*q*exp(rnorm(1, mean = 0, SDsurv))
  # C[i-1] <- B[i-1]*Fin[i-1] 
  # B[i] <- (B[i-1]-C[i-1])*(1+r*(1-(B[i-1]-C[i-1])/K))*exp(rnorm(1, mean = 0, sd = SDrun))
  
}

S = B*q*exp(rnorm(nyear, mean = 0, SDsurv))

return(data.frame(B = B, C = Fin*B, S = S, years = 1:nyear))

}