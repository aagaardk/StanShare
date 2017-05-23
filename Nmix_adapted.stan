
data {
  int<lower=1> R;                // Number of sites
  int<lower=1> T;                // Number of replications; fixed as 2
  int<lower=-1> y[R, 2, 7];      // Counts (-1:NA)
  int<lower=1,upper=7> first[R]; // First occasion
  int<lower=1,upper=7> last[R];  // Last occasion
  int<lower=0> K;                // Upper bounds of population size
}

transformed data {
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
  real<lower=0,upper=1> omega;   // Suitability
  vector[7]lam;                  // abundance
  vector[7] alpha;               // alpha of beta dist
  vector[7] beta;               // beta of beta dist
}

model {
  // Priors
  // Implicit flat priors [0, 1] are used on omega and p.
  lam ~ normal(0, 10);
  
  // Likelihood
  for (i in 1:R) {
    if (max_y_site[i]) {
      real lp = bernoulli_lpmf(1 | omega);
      
      for (k in first[i]:last[i])
        lp = lp
      + beta_binomial_lpmf(y[i, 1:T, k] | max_y_site[k], alpha[k], beta[k])
      + poisson_lpmf(k | max_y_site[k]);
      target += lp;
    } else {
      real lp[2];
      
      lp[1] = bernoulli_lpmf(0 | omega);
      lp[2] = bernoulli_lpmf(1 | omega);
      for (k in first[i]:last[i])
        lp[2] = lp[2]
      + beta_binomial_lpmf(y[i, 1:T, k] | max_y_site[k], alpha[k], beta[k])
      + poisson_lpmf(k | max_y_site[k]);
      target += log_sum_exp(lp);
    }
  }
}

