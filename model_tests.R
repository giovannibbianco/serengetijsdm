# testing the model without spatial autocorrelation

# stan model: poisson_binomial_dates_pobs
# iterations ~1000

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# cmdstanr should make things quicker 

# is not on CRAN
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# load it 
library(cmdstanr)
# it requires a shell interface called cmdstan
install_cmdstan(cores = 4)
check_cmdstan_toolchain(fix = TRUE)
# need to set path manually because somehow it doesn't work automatically 
set_cmdstan_path("C:\\Users\\giova\\Documents\\.cmdstanr\\cmdstan-2.25.0")
cmdstan_version()
cmdstan_path()

# cmdstan translates a stan program to c++ and creates and executable file making things faster.
# test it with the example provided 

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file) # compiles the model 

mod$print() # prints the stan program

mod$exe_file() # shows the path to the exe file 

# FITTING THE MODEL
#
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 2,
  refresh = 500
)

# trial model works well - cmdstan was installed successfully



model_script<-("C:\\Users\\giova\\Desktop\\serengetijsdm\\poisson_binomial_dates_pobs.stan")
mod<-cmdstan_model(model_script)
mod$print()# prints the original stan program
mod$exe_file()# shows the path to the executable version of the stan program

# in order to perform mcmc sampling we need a data_list object just like rstan 

fit<-mod$sample(
  data = stan.data,
  seed = 123,
  chains = 3,
  parallel_chains = 3,
  iter_warmup = 1000,
  iter_sampling = 1000,
  thin=1)






# read in env data and scale it 

Xo <- as.matrix(read.csv("covariates_final.csv"))

# scale environmental variables
X = Xo
tmp <- which(dimnames(Xo)[[2]] == "village_dist")
X[, tmp] = scale(Xo[, tmp])
tmp <- which(dimnames(Xo)[[2]] == "tourism_dist")
X[, tmp] = scale(Xo[, tmp])
tmp <- which(dimnames(Xo)[[2]] == "river_dist")
X[, tmp] = scale(log(Xo[, tmp]))
tmp <- which(dimnames(Xo)[[2]] == "woodycover")
X[, tmp] = scale(Xo[, tmp])
tmp <- which(dimnames(Xo)[[2]] == "nitrogen_dist")
X[, tmp] = scale(Xo[, tmp])
tmp <- which(dimnames(Xo)[[2]] == "runoff")
X[, tmp] = scale(Xo[, tmp])
tmp <- which(dimnames(Xo)[[2]] == "fire_freq")
X[, tmp] = scale(Xo[, tmp])
tmp <- which(dimnames(Xo)[[2]] == "water_dist")
X[, tmp] = scale(Xo[, tmp])


# NDVI & Rain time matrices 

rain <- as.matrix(read.csv("rain_matrix.csv"))
grass <- log(as.matrix(read.csv("grass_matrix.csv")))
ndvi <- as.matrix(read.csv("ndvi_matrix.csv"))
a_ndvi <- as.matrix(read.csv("anomaly_ndvi_matrix.csv"))
d_ndvi <- as.matrix(read.csv("Delta_NDVI.csv"))



# scale ndvi and rain
for(i in 1:ncol(ndvi)){
  rain[,i]<-(rain[,i]-mean(rain))/sd(rain)
  grass[,i]<-(grass[,i]-mean(grass))/sd(grass)
  ndvi[,i]<-(ndvi[,i]-mean(ndvi))/sd(ndvi)
  a_ndvi[,i]<-(a_ndvi[,i]-mean(a_ndvi))/sd(a_ndvi)
  d_ndvi[,i]<-(d_ndvi[,i]-mean(d_ndvi))/sd(d_ndvi)
}

dist_mat <- as.matrix(read.csv("dist_matrix.csv"))


# Transect observations
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

# average detection probabilities 
p<-c(0.621,0.698,0.624,0.779,0.796,0.791,0.637,0.85,0.803,0.677,0.608,0.6,0.706,0.92,0.72,0.463,0.732,.6,.6)

# model parameters
n_sites<-dim(X)[1] #338 sites
n_sp<-dim(C)[1] # 19 species
n_pars<-ncol(X) # 11
n_dates<-ncol(rain)

# no spatial autocorrelation 5 time variables 

stan.data <- list(
  n_obs = nrow(obs),
  n_dates = n_dates,
  n_tcov = 5,               # number of covariates changing with time
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(obs$poly.id),
  K = ncol(X), 
  X = X,
  Xt1 = rain, 
  Xt2 = grass,
  Xt3 = ndvi,
  Xt4 = a_ndvi,
  Xt5 = d_ndvi,
  date = as.integer(obs$Month),
  n_max = rep(30, n_sp),
  n_s = as.integer(n_sp),
  n_t =  ncol(Traits),
  TT = Traits,
  C = C,
  ones = numeric(n_sp) + 1,
  sp = as.integer(obs$Animal.sp),
  p_obs = p
  #Dmat = dist_mat[,-1]/10000
  #g_size = obs$Total
  )


pars <- c( "b_m", "rho",  "Sigma", "z", "Z")

fit <- stan(file = 'poisson_binomial_dates_pobs.stan',
             data = stan.data,
             #init = init_f,
             pars = pars,
             iter = 1000, thin = 1, chains = 3)

summary(fit)

# model started at 16:46 19/11
