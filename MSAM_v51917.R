###MSAM using STAN in R
#created by Curtis Burkhalter 2017

#load required packages
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)
set.seed(123)

#set working directory
setwd("~/Personal GHub/STAN-MSAM")

## Read data
## The data file "fritilary.txt" is available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
bdat <- read.table("fritillary.txt", header = TRUE)
y <- array(NA, dim = c(95, 2, 7))	# 95 sites, 2 reps, 7 days

for(k in 1:7) {
  sel.rows <- bdat$day == k
  y[, , k] <- as.matrix(bdat)[sel.rows, 3:4]
}
R = nrow(y)
T = ncol(y)
first <- sapply(1:dim(y)[1], function(i)
  min(grep(FALSE, is.na(y[i, 1, ]))))
last <- sapply(1:dim(y)[1], function(i)
  max(grep(FALSE, is.na(y[i, 1, ]))))
y[is.na(y)] <- -1

## Parameters monitored
params <- c("N","lambda","alpha","beta")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 1

## Initial values
inits <- lapply(1:nc, function(i)
  list(lambda = runif(7,1,1),N = matrix(runif(R*7,1,1),nrow = 95,ncol =7, byrow=TRUE),alphaP = runif(7,1,1), betaP = runif(7,1,1)))

## Call Stan from R
out1 <- stan("Nmix_adapted_noZI.stan",
             data = list(y = y, R = R, T = 2,K = 100),
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             control = list(adapt_delta = 0.9),
             open_progress = FALSE)
print(out1, digits = 3)