data {
	int<lower=0> ncounts;
	int<lower=0> nregions;
	int<lower=0> nyears;
	int<lower=0> nobservers;
	int<lower=0> nunits;
	int<lower=0> nmigday;
	int<lower=0> counts[ncounts];
	int<lower=0,upper=nregions> regions[ncounts];
	int<lower=0,upper=nyears> years[ncounts];
	int<lower=0,upper=nobservers> observers[ncounts];
	int<lower=0,upper=nunits> units[ncounts];
	int<lower=0,upper=nmigday> migday[ncounts];
	real area[ncounts];
	real mean_count;
	vector[ncounts] ScaledEffort;
}

parameters {
	//Regression slopes
	real reg[nregions];
	real yr[nyears];
	real obs[nobservers];
	real unit[nunits];
	real hectares;
	real day;
	
	//Average-count intercept
  real count_intercept;
	
	//Scale parameters
	vector[nregions] B;
	vector[nregions] p;
	real B_mu;
	real p_mu;

  //Standard deviations
  real<lower=0> sigma;
	real<lower=0> B_sigma;
	real<lower=0> p_sigma;
	real<lower=0> count_sigma;
	
	//Linear predictor
	real lambda[ncounts];
}

transformed parameters {
	//Quadratic-day parameter
	real quadday;
  
	//Effort transformed parameters
	vector[ncounts] trans_effort;
	vector[ncounts] effort;
	real dBp_mu[nregions];
	real dBp_zero_one[nregions];

	//Non-linear day effects
	quadday = day^2;
	
  //Effort transformation
	for(i in 1:ncounts){
		trans_effort[i] = ((ScaledEffort[i] ^ p[regions[i]]) - 1) / 
								p[regions[i]];
		effort[i] = B[regions[i]] * trans_effort[i];
	}
	
	for(i in 1:nregions){
	  dBp_mu[i] = sqrt((p[i] - p_mu)^2 + (B[i] - B_mu)^2);
  	dBp_zero_one[i] = sqrt((p[i] - 0)^2 + (B[i] - 1)^2);
	}
}

model {
	//Effort Effects
	B_mu ~ normal(0.0,1000);
	p_mu ~ normal(0.0,1000);
	for(i in 1:nregions){
		B[i] ~ normal(B_mu, B_sigma);
		p[i] ~ normal(p_mu, p_sigma);
  }
	
	//Average-count distribution
	count_intercept ~ normal(mean_count,count_sigma);
	
	//Likelihood function
	for(i in 1:ncounts){
		lambda[i] ~ 
		  normal(
		    (hectares*area[i]) + count_intercept 
		    + (day*migday[i]) + (quadday*migday[i])
		    + yr[years[i]] + reg[regions[i]] 
		    + unit[units[i]] + obs[observers[i]] 
		    + effort[i],sigma);
	}
	counts ~ poisson_log(lambda);
}

generated quantities {
  //Count transformation parameter
  real<lower=0> pred_counts[ncounts];
  
	//Bayesian p-value parameters
  real dBppost[nregions];
	
	//Exponentiated count estimates
	for(i in 1:ncounts){
    pred_counts[i] = exp(lambda[i]);
  }

	//Bayesian p-value calculation
  for(i in 1:nregions){
    dBppost[i] = dBp_mu[i] - dBp_zero_one[i];
  }
}
