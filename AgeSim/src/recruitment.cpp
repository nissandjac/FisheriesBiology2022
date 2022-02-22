// Create a file to run the TMB version of sprat
#include <TMB.hpp>
#include <iostream>


template<class Type> // TMB function to run a stock assessment on the MSE
Type objective_function<Type>::operator() ()
{
  // Data input
   DATA_VECTOR(SSBn); // Spawning Biomass
   DATA_VECTOR(Robs); // Observed recruitment
   DATA_INTEGER(nobs); // Number of observations
  // //
  PARAMETER(logh)
  PARAMETER(logSSBzero)
  PARAMETER(logRzero)
  PARAMETER(logSDR)

  Type Rzero = exp(logRzero);
  Type SSBzero = exp(logSSBzero);
  Type steep = exp(logh);
  Type SDR = exp(logSDR);

  // Pre assign stuff
  vector<Type> ml(nobs);
  Type ans = 0.0;

  //  Calculate SR relationship
  // //
  for(int i=0;i<nobs;i++){ // Start time loop
   ml(i) = Type(4.0)*steep*Rzero*SSBn(i)/(SSBzero*(Type(1.0)-steep)+ SSBn(i)*(Type(5.0)*steep-Type(1.0)));
  }

 // Likelihood function
//
   for(int i=0;i<nobs;i++){ // Start time loop
    ans += -dnorm(log(ml(i)+0.001),log(Robs(i)+0.001),SDR,TRUE);
     }

//


// Report calculations
REPORT(ml)
REPORT(ans)
ADREPORT(ml)


  return ans;
}
