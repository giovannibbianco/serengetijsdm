

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("ser1.RData") # in environment it is read as fits



fit_summary <- summary(fits)$summary 

# responses to environmental variables
responses <- fit_summary[grepl("b_m", rownames(fit_summary)),]
which(sign(responses[,4])==sign(responses[,8])) # which ones do not cross the 0 line with their credible interval 

# traits
traits_effects<-fit_summary[grepl("Z", rownames(fit_summary)),]
which(sign(traits_effects[,4])==sign(traits_effects[,8]))

#r-hat and effective sample size
hist(fit_summary[,10], main = "R-hat Std Model")
hist(fit_summary[,9], main = "n-eff Std Model" )

# phylogenetic signal 
draws <- extract(fits, pars = "rho")
plot(density(draws$rho), main="Rho-Phylogenetic Signal")
