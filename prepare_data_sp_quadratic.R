library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# prepare matrices
# 

X<-read.csv("covariates_quadratic.csv")
X<-as.matrix(X)
X

# scale environmental variables
for(i in 2:ncol(X)){
  X[,i]<-scale(X[,i])
}

# NDVI & Rain

ndvi<-as.matrix(read.csv("ndvi_matrix.csv"))
rain<-as.matrix(read.csv("rain_matrix.csv"))

# scale ndvi and rain
for(i in 1:ncol(ndvi)){
  ndvi[,i]<-(ndvi[,i]-mean(ndvi))/sd(ndvi)
}

for(i in 1:ncol(rain)){
  rain[,i]<-(rain[,i]-mean(rain))/sd(rain)
}


# Observations
obs<-read.csv("obs_matrix.csv")
names(obs)

obs$Animal.sp<-as.factor(obs$Animal.sp)
levels(obs$Animal.sp) #19 species 

# recode species ids progressively
levels(obs$Animal.sp)<-c(1:19)
obs$Month<-as.factor(obs$Month) # range from 1:12 like dims of ndvi and rain


#Species Specific Traits
TT<-read.csv("TT_matrix.csv")
names(TT)

# add intercept for traits
TT[,1]<-1
names(TT)[1]<-"intercept"

TT$browser = as.numeric(TT$guild=="browser")
TT$mix = as.numeric(TT$guild=="mixed")

Traits = cbind(TT[,1], scale(log(TT[,4])), scale(TT[,5]), TT[,6:7])


# Phylogeny
library(ape)
filo<-read.nexus("phylotree.nex")
CC=vcv(filo,corr=T)
CC

# sort species in order
dimnames(CC)

ids<-c(12,9,2,18,17,10,13,19,11,15,7,8,14,3,4,16,1,6,5) # species ids changed

C <- matrix(NA, ncol(CC), ncol(CC))

for(i in 1:ncol(CC)){
  for(j in 1:ncol(CC)){
    C[ids[i],ids[j]] <- CC[i,j]
  }
}

C

# detection probabilities 
p<-c(0.621,0.698,0.624,0.779,0.796,0.791,0.637,0.85,0.803,0.677,0.608,0.6,0.706,0.92,0.72,0.463,0.732,.6,.6)

# model parameters
n_sites<-dim(X)[1] #338 sites
n_sp<-dim(C)[1]
n_pars<-ncol(X) # 11 parameters including intercept
n_dates<-ncol(rain)

table(obs$Animal.sp)

# Species 1 - Wildebeest

sp = 1
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
            data = stan.data,
            pars = pars,
            iter = 1000, thin = 1, chains = 3)


fit_summary <- summary(fit)$summary 

mean_b <- fit_summary[grepl("betas", rownames(fit_summary)),] # creates subset of just regression parameters








# Zebra

sp = 2
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit2 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
             data = stan.data,
             pars = pars,
             iter = 1000, thin = 1, chains = 3)



fit_summary2 <- summary(fit2)$summary 

mean_b2 <- fit_summary2[grepl("betas", rownames(fit_summary2)),] # creates subset of just regression parameters

# Thomson Gazelle

sp = 3
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit3 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
             data = stan.data,
             pars = pars,
             iter = 1000, thin = 1, chains = 3)



fit_summary3 <- summary(fit3)$summary 

mean_b3 <- fit_summary3[grepl("betas", rownames(fit_summary3)),] # creates subset of just regression parameters


# Species 4 Grant's gazelle

sp = 4
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit4 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
             data = stan.data,
             pars = pars,
             iter = 1000, thin = 1, chains = 3)



fit_summary4 <- summary(fit4)$summary 

mean_b4 <- fit_summary4[grepl("betas", rownames(fit_summary4)),] # creates subset of just regression parameters

# Species 5 Topi


sp = 5
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit5 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
             data = stan.data,
             pars = pars,
             iter = 1000, thin = 1, chains = 3)



fit_summary5 <- summary(fit5)$summary 

mean_b5 <- fit_summary5[grepl("betas", rownames(fit_summary5)),] # creates subset of just regression parameters

# species 6 Kongoni

sp = 6
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit6 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
             data = stan.data,
             pars = pars,
             iter = 1000, thin = 1, chains = 3)



fit_summary6 <- summary(fit6)$summary 

mean_b6 <- fit_summary6[grepl("betas", rownames(fit_summary6)),] # creates subset of just regression parameters


# species 7 - Impala

sp = 7
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit7 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
             data = stan.data,
             pars = pars,
             iter = 1000, thin = 1, chains = 3)



fit_summary7 <- summary(fit7)$summary 

mean_b7 <- fit_summary7[grepl("betas", rownames(fit_summary7)),] # creates subset of just regression parameters

# Species 8 Waterbuck

sp = 8
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit8 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
             data = stan.data,
             pars = pars,
             iter = 1000, thin = 1, chains = 3)



fit_summary8 <- summary(fit8)$summary 

mean_b8 <- fit_summary8[grepl("betas", rownames(fit_summary8)),] # creates subset of just regression parameters

# species 9 Warthog

sp = 9
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit9 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
             data = stan.data,
             pars = pars,
             iter = 1000, thin = 1, chains = 3)



fit_summary9 <- summary(fit9)$summary 

mean_b9 <- fit_summary9[grepl("betas", rownames(fit_summary9)),] # creates subset of just regression parameters


# species 10 - Buffalo

sp = 10
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit10 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
              data = stan.data,
              pars = pars,
              iter = 1000, thin = 1, chains = 3)



fit_summary10 <- summary(fit10)$summary 

mean_b10<- fit_summary10[grepl("betas", rownames(fit_summary10)),] # creates subset of just regression parameters




# Species 11 Giraffe


sp = 11
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit11 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
              data = stan.data,
              pars = pars,
              iter = 1000, thin = 1, chains = 3)



fit_summary11 <- summary(fit11)$summary 

mean_b11<- fit_summary11[grepl("betas", rownames(fit_summary11)),] # creates subset of just regression parameters

# Species 12 - Elephant

sp = 12
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit12 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
              data = stan.data,
              pars = pars,
              iter = 1000, thin = 1, chains = 3)



fit_summary12 <- summary(fit12)$summary 

mean_b12<- fit_summary12[grepl("betas", rownames(fit_summary12)),] # creates subset of just regression parameters

# Species 13 Eland

sp = 13
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit13 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
              data = stan.data,
              pars = pars,
              iter = 1000, thin = 1, chains = 3)



fit_summary13 <- summary(fit13)$summary 

mean_b13<- fit_summary13[grepl("betas", rownames(fit_summary13)),] # creates subset of just regression parameters

# Species 16 DikDik


sp = 16
obssp = obs[obs$Animal.sp==sp, ]

stan.data <- list(
  n_obs = dim(obssp)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obssp$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obssp$Month),
  n_max = 50,
  p_obs = p[sp]
)

pars <- c( "betas", "D", "p")


fit16 <- stan(file = 'poisson_binomial_dates_pobs_sp.stan',
              data = stan.data,
              pars = pars,
              iter = 1000, thin = 1, chains = 3)



fit_summary16 <- summary(fit16)$summary 

mean_b16<- fit_summary16[grepl("betas", rownames(fit_summary16)),] # creates subset of just regression parameters

# Make Posterior Table

# Calculate fraction of the posterior same sign as the mean 
B = extract(fit16, pars = "betas")

f = numeric(ncol(B$betas))

for(i in 1: length(f)){
  f[i] = length(which(sign(B$betas[,i]) == sign(mean(B$betas[,i]))))/nrow(B$betas)
}
f

post1<-f
post2<-f
post3<-f
post4<-f
post5<-f
post6<-f
post7<-f
post9<-f
post10<-f
post11<-f
post12<-f
post13<-f
post16<-f

posteriors.signs<-c(post1,post2,post3,post4,post5,post6,post7,post9,post10,post11,post12,post13,post16)# add one species column and 
# names of variables on top

posterior.table<-rbind(mean_b,mean_b2,mean_b3,mean_b4,mean_b5,mean_b6,mean_b7,mean_b9,mean_b10,mean_b11,mean_b12,mean_b13,mean_b16)

full_posterior_table<-cbind(posterior.table,posteriors.signs)

write.csv(full_posterior_table,file="posterior_table_quadratic.csv")

species.name