
data {
  int<lower=1> R;                   //Number of sites
  int<lower=1> T;                   //Number of replications; fixed as 2                                       
  int<lower=-1> y[R, 2, 7];         //  Counts (-1:NA)
  int<lower=1,upper=7> first[R];   //  First occasion
  int<lower=1,upper=7> last[R];   //  Last occasion
  int<lower=0> K;                // Upper bounds of population size             
}

transformed data {
  real<lower=0> lambda_mu;           
  real<lower=0> lambda_sig;        
  lambda_mu = 0;
  lambda_sig = 10;
  int<lower=0> max_y[R, 7];
  int<lower=0> max_y_site[R];
  
  for (i in 1:R) {
    for (k in 1:(first[i] - 1))
      max_y[i, k] = 0;
    for (k in (last[i] +1 ):7)
      max_y[i, k] = 0;
    for (k in first[i]:last[i])
      max_y[i, k] = max(y[i, 1:T, k]);
      max_y_site[i] = max(max_y[i]);
  }
}

parameters {
  matrix [R,7] N;               // latent abundance for each site x day
  vector[7] alphaP;               // alpha of beta dist for det. prob
  vector[7] betaP;               // beta of beta dist for det. prob
}                               // use beta to model det. prob b/c of correlated                                 // behavior in overwintering grassland birds

model {
  // Priors

  // Likelihood
  for (i in 1:R) {
        for (k in first[i]:last[i]) {
        lambda[k] ~ normal(lambda_mu,lambda_sig);
        poisson_lpmf(N[i,k] | lambda[k]) 
        + beta_binomial_lpmf(y[i, 1:T, k] | N[i,k], alphaP[k], betaP[k]);
        target += log_sum_exp(lp);
         
        }
  }
}

generated quantities {
  
  for (k in 1:7) {
    totalN[k] = sum(N[,k])
    p.derived[k] = alphaP[k]/(alphaP[k]+betaP[k]) #derived detection probability
    rho.derived[k] = 1/(alphaP[k]+alphaP[k])   #derived correlation coefficient
    }
}
