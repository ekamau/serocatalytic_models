// Model - elevated death rates due to infection
// stepwise FOI - (i) one fixed lambda; one estimated lambda estimate; 
// (ii) using eqn 41; (iii) estimated mu (iv) fixed IFR/rho = 0.89;

functions{
  real prob_seropos_calc(real lambda, real mu, real rho, int age) {
    vector[age+1] vecS;
    vector[age+1] vecXm;
    vector[age+1] vecXs;
    real prob;
    
    vecS[1] = 1; vecXm[1] = 0; vecXs[1] = 0;
    
    for(i in 1:age){
      int yr = (2017-age)+i; // year at age i
      real foi;
      if(yr < 2014){
        foi = 0.0;
      } else {
        foi = lambda;
      }
        vecS[i+1] = (vecS[i]*exp(-foi));
        vecXm[i+1] = (vecXm[i] + vecS[i]*(1-rho)*(1-exp(-foi)));
        vecXs[i+1] = (vecS[i]*rho*foi*(exp(-mu)-exp(-foi)) + vecXs[i]*(foi-mu)*exp(-mu))/(foi-mu);
    }
    
    prob = (vecXm[age+1] + vecXs[age+1])/(vecXm[age+1] + vecXs[age+1] + vecS[age+1]);
    return prob;
  }

}

data{
  int<lower=0> n_obs; // No. rows in data or no. age classes
  int n_pos[n_obs]; // seropositive
  int n_total[n_obs]; // tested
  int ages[n_obs];
  
}

parameters{
  real<lower=0> lambda; // real<lower=0, upper=gamma>;
  real<lower=0> time_to_die; // death

}

transformed parameters{
  real<lower=0> prop_seropos[n_obs];
  real<lower=0> mu = 1/time_to_die;
  real<lower=0> rho = 0.89; // IFR
  
  for(i in 1:n_obs){
    int age = ages[i];
    prop_seropos[i] = prob_seropos_calc(lambda, mu, rho, age);
  }
    
}

model{
  // priors
  lambda ~ exponential(1);
  time_to_die ~ normal(0.04, 0.02);
   
  // likelihood
  n_pos ~ binomial(n_total, prop_seropos);
}

generated quantities {
  int pos_pred[n_obs] = binomial_rng(n_total, prop_seropos);
  
}
