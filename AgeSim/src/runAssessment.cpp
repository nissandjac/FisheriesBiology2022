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
  DATA_SCALAR(logSDC); // Catch Uncertainty
  DATA_SCALAR(logSDS); // Survey uncertainty

  PARAMETER_VECTOR(logFzero); // Catch estimation (maybe change to exploitation rate)
  PARAMETER_VECTOR(logBpred);

  PARAMETER(logr); // growth rate
  PARAMETER(logK); // Carrying capacity
  PARAMETER(logSDB); // Observation error
  PARAMETER(logq); // survey catchability
  PARAMETER(logm); // Production shape parameter
  //DATA_SCALAR(logq); //
  //PARAMETER(logSDcatch); // Standard deviation of catch
  // PARAMETER(logSDsurv); // Standard deviation of survey
  //DATA_SCALAR(logSDcatch); // Standard deviation of catch
  //PARAMETER(logSDsurv); // Standard deviation of survey
  PARAMETER(logSDF); // fishing mortality error
  //DATA_SCALAR(logSDF); // fishing mortality error

// Biological parameters
  Type K = exp(logK);
  Type r = exp(logr); // Keep between 0 and 1
  Type q = exp(logq);
  Type m = exp(logm);
// Uncertainty parameters
  Type SDB = exp(logSDB);
  Type SDF = exp(logSDF);
  Type SDC = exp(logSDC);
  Type SDS = exp(logSDS);

  // Catch
  vector<Type> Fzero = logit(logFzero); // Keep exploitation rate between [0;1]
  vector<Type> Bpred = exp(logBpred); // Expected biomass with error
  // Initialize simulation
  vector<Type> B(nyear); // Biomass
  vector<Type> Catch(nyear);
  vector<Type> surv_est(nyear);

  B(0) = Bpred(0)-Fzero(0)*Bpred(0);

  for(int time=0;time<(nyear-1);time++){
    Catch(time) = Fzero(time)*Bpred(time);

    B(time+1) = Bpred(time)+((m+1)/m)*Bpred(time)*r*(1-pow(Bpred(time)/K, m))-Catch(time); // Deterministic model
    surv_est(time) = Bpred(time)*q;

    B(time+1) = CppAD::CondExpLe(B(time+1), Type(0), Type(0.000001), B(time+1)); // Make sure too low K and too high Bpred doesn't destroy the convergence

    // if(B(time+1) < 0){
    //   std::cout << "Biomass next year = "<< B(time+1) << " B/K is :" << B/K << std::endl;
    //  }
  }

  Catch(nyear-1) = Bpred(nyear-1)*Fzero(nyear-1);
  surv_est(nyear-1) = Bpred(nyear-1)*q;

  // Debug the NA values //

  Type ans_survey=0.0;
  Type ans_catch = 0.0;
  Type ans_prior = 0.0;
  Type ans_rand = 0.0;
  Type ans_prior2 = 0.0;

  for(int i=0;i<nyear;i++){
     ans_survey += -dnorm(log(surv_est(i)+0.001),log(survey(i)+0.001), SDB, TRUE); // Assume that survey and biomass have the same uncertainty 

     //if(Catch(i) > 0){
     ans_catch += -dnorm(log(Catch(i)+0.001),log(Catchobs(i)+0.001), SDC, TRUE);
      //}

     ans_rand += -dnorm(log(Bpred(i)+0.001), log(B(i)+0.001), SDB, TRUE);

  }

  for(int i=1;i<nyear;i++){
         ans_prior += -dnorm( log(Fzero(i)+0.0001), log(Fzero(i-1)+0.0001), SDF, TRUE ); //
     }


   ans_prior2 += -dnorm(log(m), log(Type(.188)), Type(.5), TRUE);


  // Print the ans values
  // std::cout << "survey = "<< ans_survey << " catch = " << ans_catch << "prior = " << ans_prior << std::endl;
  // std::cout << "random effect on B = "<< ans_rand << " prior on m = " << ans_catch << "random effect on F = " << ans_prior2 << std::endl;

   // ans_prior2 += -dnorm(log(q), log(Type(.5)), Type(.5), TRUE);
  //
  // if(SDsurv > SDB){
  //   ans_rand += 1e8;
  // }

  //
  // for(int i=0;i<(nyear-1);i++){
  // }

  Type ans = ans_survey+ans_rand+ans_catch+ans_prior+ans_prior2;
  vector<Type> ans_vec(5);

  ans_vec(0) = ans_survey;
  ans_vec(1) = ans_catch;
  ans_vec(2) = ans_prior;
  ans_vec(3) = ans_rand;
  ans_vec(4) = ans_prior2;



  //std::cout << "what value does ans have"<< ans << std::endl;




  REPORT(B)
  REPORT(Bpred)
  REPORT(Catch)
  REPORT(ans_vec)
  REPORT(surv_est)
  REPORT(ans)
  REPORT(Fzero)
  REPORT(r)

  ADREPORT(B)
  ADREPORT(Bpred)
  ADREPORT(ans_vec)
  ADREPORT(Catch)
  ADREPORT(surv_est)
  ADREPORT(Fzero)

  return ans;
}
