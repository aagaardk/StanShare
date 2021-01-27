data {
  int<lower=1> Nobs;
  real Observed[Nobs];
  int<lower=1> Nind;
  int<lower=1, upper = Nind> Individual[Nobs];
  int<lower=0, upper = 1> CV2[Nobs];
  real Batch2[Nobs];
  real Batch3[Nobs];
  real Batch4[Nobs];
  real Batch5[Nobs];
  int<lower = 1> Naocs;
  int<lower = 1> Nsites;
  int<lower = 1, upper = Naocs>  AOC[Nsites];
  int<lower = 1, upper = Nsites> Site[Nind];
  int<lower = 0, upper = 1> Chick[Nind];
}

parameters {
  real indmean[Nind];
  real<lower = 0> indsigma;
  real batch2mean;
  real batch3mean;
  real batch4mean;
  real batch5mean;
  real CV2mean;
  real Chickmean;
  real<lower=0> siteSigma;
  real siteMean[Nsites];
  real AOCmean[Naocs];
  real<lower=0> aocSigma;
}

model {
  siteMean     ~ normal(2, 10);
  siteSigma    ~ cauchy(0, 2.5);
  indsigma     ~ cauchy(0, 2.5);
  aocSigma     ~ cauchy(0, 2.5);
  AOCmean      ~ normal(2, 10);
  Chickmean    ~ normal(0, 10);
  CV2mean      ~ normal(0, 10);
  batch2mean   ~ normal(0, 10);
  batch3mean   ~ normal(0, 10);
  batch4mean   ~ normal(0, 10);
  batch5mean   ~ normal(0, 10);
  
  for(s in 1:Nsites){
    siteMean[s] ~ normal( AOCmean[AOC[s]],
  			  aocSigma);
  }
  
  for(i in 1:Nind){
    indmean[i] ~ normal( siteMean[Site[i]]
			 + Chickmean * Chick[i]
  			 , siteSigma);
  }
  
  for(o in 1:Nobs){
    Observed[o] ~ normal( indmean[Individual[o]] +
  			  Batch2[o] * batch2mean +
  			  Batch3[o] * batch3mean +
  			  Batch4[o] * batch4mean +
			  Batch5[o] * batch5mean +
  			  CV2[o] * CV2mean
			  , indsigma);
  }
}
