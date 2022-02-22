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
  PARAMETER_VECTOR(logBpred);
  PARAMETER_VECTOR(logr); // growth rate

  PARAMETER(logK); // Carrying capacity
  PARAMETER(logSDB); // Observation error
  PARAMETER(logSDr);
  PARAMETER(logq); // survey catchability
  //DATA_SCALAR(logq); //
  // PARAMETER(logSDcatch); // Standard deviation of catch
  // PARAMETER(logSDsurv); // Standard deviation of survey
  DATA_SCALAR(logSDcatch); // Standard deviation of catch
  DATA_SCALAR(logSDsurv); // Standard deviation of survey
  // PARAMETER(logSDF); // fishing mortality error
  DATA_SCALAR(logSDF); // fishing mortality error

// Biological parameters
  Type K = exp(logK);
  Type q = exp(logq);
// Uncertainty parameters
  Type SDR = exp(logSDr);
  Type SDB = exp(logSDB);
  Type SDcatch = exp(logSDcatch);
  Type SDsurv = exp(logSDsurv);
  Type SDF = exp(logSDF);
// Catch
  vector<Type> Fzero = invlogit(invlogit_Fzero); // Keep exploitation rate between [0;1]
  vector<Type> Bpred = exp(logBpred); // Expected biomass with error
  vector<Type> r = exp(logr);

// Initialize simulation
  vector<Type> B(nyear); // Biomass
  vector<Type> Catch(nyear);
  vector<Type> surv_est(nyear);

  B(0) = Bpred(0)-Bpred(0)*Fzero(0);

  for(int time=0;time<(nyear-1);time++){
    Catch(time) = Fzero(time)*Bpred(time);
    B(time+1) = Bpred(time)+r(time)*Bpred(time)*(1-Bpred(time)/K)-Catch(time); // Deterministic model
    surv_est(time) = Bpred(time)*q;
  }

  Catch(nyear-1) = Bpred(nyear-1)*Fzero(nyear-1);
  surv_est(nyear-1) = Bpred(nyear-1)*q;

  Type ans_survey=0.0;
  Type ans_catch = 0.0;
  Type ans_prior = 0.0;
  Type ans_rand = 0.0;
  Type ans_r = 0.0;

  for(int i=0;i<nyear;i++){
     ans_survey += -dnorm(log(surv_est(i)+0.001),log(survey(i)+0.001), SDsurv, TRUE);
     ans_catch += -dnorm(log(Catch(i)+0.001),log(Catchobs(i)+0.001), SDcatch, TRUE);


  }

  for(int i=1;i<nyear;i++){
    ans_prior += -dnorm( log(Fzero(i)), log(Fzero(i-1)), SDF, TRUE ); //
    ans_r += -dnorm( log(r(i)), log(r(i-1)), SDF, TRUE ); //
    ans_rand += -dnorm(log(Bpred(i)+0.001), log(B(i)+0.001), SDR, TRUE);


  }

  //
  // for(int i=0;i<(nyear-1);i++){
  // }

  Type ans = ans_survey+ans_rand+ans_catch+ans_prior+ans_r;
  vector<Type> ans_vec(5);

  ans_vec(0) = ans_survey;
  ans_vec(1) = ans_catch;
  ans_vec(2) = ans_prior;
  ans_vec(3) = ans_rand;
  ans_vec(4) = ans_r;

  REPORT(B)
  REPORT(Bpred)
  REPORT(Catch)
  REPORT(ans_vec)
  REPORT(surv_est)
  REPORT(ans)
  REPORT(r)

  ADREPORT(B)
  ADREPORT(Bpred)
  ADREPORT(ans_vec)
  ADREPORT(Catch)
  ADREPORT(surv_est)
  ADREPORT(Fzero)
  ADREPORT(r)

  return ans;
}
