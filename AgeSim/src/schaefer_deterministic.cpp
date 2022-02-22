// Create a file to run the TMB version of sprat
#include <TMB.hpp>
#include <iostream>


template <class Type>
vector<Type> cumsum(vector<Type> x) {
  int n = x.size();
  vector<Type> ans(n);
  ans[0] = x[0];
  for (int i = 1; i < n; i++) ans[i] = x[i] + ans[i-1];
  return ans;
}


template<class Type> // TMB function to run a stock assessment on the MSE
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(Catchobs); // Weight in the beginning of the year
  DATA_VECTOR(survey); // Weight in catch
  DATA_INTEGER(nyear); // Number of years to simulate

  PARAMETER_VECTOR(invlogit_Fzero); // Catch estimation (maybe change to exploitation rate)
  PARAMETER(logr); // growth rate
  PARAMETER(logK); // Carrying capacity
  PARAMETER(logq); // survey catchability
  //PARAMETER(logSDcatch); // Standard deviation of catch
  // PARAMETER(logSDF); // fishing mortality error

  // PARAMETER(logSDsurv); // Standard deviation of survey
  DATA_SCALAR(logSDcatch); // Standard deviation of catch
  DATA_SCALAR(logSDsurv); // Standard deviation of survey
  //DATA_SCALAR(logSDF); // fishing mortality error

// Biological parameters
  Type K = exp(logK);
  Type r = exp(logr);
  Type q = exp(logq);
// Uncertainty parameters
  Type SDcatch = exp(logSDcatch);
  Type SDsurv = exp(logSDsurv);
// Catch
  vector<Type> Fzero = invlogit(invlogit_Fzero); // Keep exploitation rate between [0;1]
// Initialize simulation
  vector<Type> B(nyear); // Biomass
  vector<Type> Catch(nyear);
  vector<Type> surv_est(nyear);

  B(0) = K-Fzero(0)*K;

  for(int time=0;time<(nyear-1);time++){
    Catch(time) = Fzero(time)*B(time);
    B(time+1) = B(time)+r*B(time)*(1-B(time)/K)-Catch(time); // Deterministic model
    surv_est(time) = B(time)*q;

  }

  Catch(nyear-1) = B(nyear-1)*Fzero(nyear-1);
  surv_est(nyear-1) = B(nyear-1)*q;

  Type ans_survey=0.0;
  Type ans_catch = 0.0;
  Type ans_prior = 0.0;
  Type ans_rand = 0.0;

  for(int i=0;i<nyear;i++){
     ans_survey += -dnorm(log(surv_est(i)+0.001),log(survey(i)+0.001), SDsurv, TRUE);
     ans_catch += -dnorm(log(Catch(i)+0.001),log(Catchobs(i)+0.001), SDcatch, TRUE);
  }

  // for(int i=1;i<nyear;i++){
  //
  //   ans_prior += -dnorm( log(Fzero(i)), log(Fzero(i-1)), SDF, TRUE ); //
  // }

  //
  // for(int i=0;i<(nyear-1);i++){
  // }

  Type ans = ans_survey+ans_rand+ans_catch+ans_prior;
  vector<Type> ans_vec(4);

  ans_vec(0) = ans_survey;
  ans_vec(1) = ans_catch;
  ans_vec(2) = ans_prior;
  ans_vec(3) = ans_rand;

  REPORT(B)
  REPORT(Catch)
  REPORT(ans_vec)
  REPORT(surv_est)
  REPORT(ans)

  ADREPORT(B)
  ADREPORT(ans_vec)
  ADREPORT(Catch)
  ADREPORT(surv_est)
  ADREPORT(Fzero)

  return ans;
}
