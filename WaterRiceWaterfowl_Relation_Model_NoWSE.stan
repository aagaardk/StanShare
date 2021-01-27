data {
	int<lower=0> ndata;
	int years[ndata];
	real rd_mean;
	real wa_mean;
	real TWD_mean;
	real<lower=0> abundance[ndata];
	real<lower=0> waterdepth[ndata];
	real<lower=0> rice[ndata];
}

parameters { 
  real yr[ndata];

	real ricedens_sigma;
	real wfabund_sigma;
	real<lower=0> ricedens_overdisp;
	real<lower=0> wfabund_overdisp;
	
	real ricedens_intercept;
	real wfabund_intercept;
	real TWD_beta;
	real ricedens_wfabund_beta;
	
	real<lower=0> ricedens_est[ndata];
	real<lower=0> wfabund_est[ndata];
}

model {
  //Random Intercepts
	ricedens_intercept ~ normal(rd_mean,10); 
	wfabund_intercept ~ normal(wa_mean,10); 

	//Water depth, and 
	//rice density coefficients
	TWD_beta ~ normal(TWD_mean,0.1);
	ricedens_wfabund_beta ~ normal(rd_mean,0.1);

  //Random rice density year effects
  for(i in 1:ndata){
 		yr[i] ~ normal(0.0,1); 
	}

	//Rice density model
	for(i in 1:ndata){
	  ricedens_est[i] ~ normal(ricedens_intercept 
	                        + yr[i] 
	                        + (TWD_beta * waterdepth[i]), 
	                        ricedens_overdisp);
		rice[i] ~ normal(ricedens_est[i], ricedens_sigma);
	}

	//Waterfowl Abundance Model
	for(i in 1:ndata){
	  wfabund_est[i] ~ normal(wfabund_intercept 
	                    + yr[i]
	                    + (ricedens_wfabund_beta * ricedens_est[i]),
	                    wfabund_overdisp);
		abundance[i] ~ normal(wfabund_est[i], wfabund_sigma);
	}
}

