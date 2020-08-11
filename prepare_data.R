# prepare matrices
# 

X<-read.csv("covariates.csv")
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
p<-c(0.621,0.698,0.624,0.779,0.796,0.791,0.637,0.85,0.803,0.677,0.608,0.6,0.706,0.92,0.72,0.463,0.732,1,1)

# model parameters
n_sites<-dim(X)[1] #338 sites
n_sp<-dim(C)[1]
n_pars<-ncol(X) # 11 parameters including intercept
n_dates<-ncol(rain)

# 

stan.data <- list(
  n_obs = dim(obs)[1],
  n_dates = n_dates,
  n_tcov = 2,
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obs$poly.id),
  K = dim(X)[2], 
  X = X,
  Xt1 = ndvi, # we should add also ndvi
  Xt2 = rain,
  date = as.integer(obs$Month),
  n_max = rep(30, n_sp),
  n_s = as.integer(n_sp),
  n_t =  dim(Traits)[2],
  TT = Traits,
  C = C,
  ones = numeric(n_sp) + 1,
  sp = as.integer(obs$Animal.sp),
  p_obs = p
)

pars <- c( "b_m", "rho",  "Sigma", "z", "Z")


init_f <- function () list(b_m = matrix(0, n_sp, n_pars+2),
                           m = numeric(n_sp * (n_pars+2)))
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fit <- stan(file = 'poisson_binomial_dates_pobs.stan',
            data = stan.data,
            init = init_f,
            pars = pars,
            iter = 100, thin = 1, chains = 1)



