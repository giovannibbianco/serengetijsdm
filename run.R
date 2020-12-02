library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("dat.RData")

pars <- c( "b_m", "rho",  "Sigma", "z", "Z", "rhosq", "etasq", "delta")

fits <- stan(file = 'poisson_binomial_dates_pobs_se1.stan',
             data = stan.data,
             #init = init_f,
             pars = pars,
             iter = 2000, thin = 1, chains = 4)