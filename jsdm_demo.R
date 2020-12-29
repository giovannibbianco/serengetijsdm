# first we generate mock data

# a number of species 
# a phylogenetic tree indicating how these species are related 
# a set of species specific traits for every species 
# a set of environmental covariates (2 changing across space, 2 changing across space and time)


library(ape) # for phylogenetic data
library(mvtnorm) # to compute multivariate probabilities

set.seed(1234)

n_sp=20 # number of species

n_env=2 # environmental covariates 

n_pars=n_env+1 #parameters (environmental cov + intercept)

n_t=2 # number of species specific traits 

# defining species specific  traits 
dgrass=rbeta(n_sp,1,1) # fraction of grass in diet 
log_bm=rnorm(n_sp,0,1) # log of body mass

#arrange species specific traits in a matrix + intercept
TT = as.matrix(cbind(rep(1, n_sp), scale(dgrass), scale(log_bm)))


# simulate phylogeny 
tree=rtree(n_sp)
CC=vcv(tree, corr=T) 


# sort species and re-arrange phylogenetic correlation matrix
tmp = dimnames(CC)
ids = as.numeric(as.factor(tmp[[1]]))

C = matrix(NA, ncol(CC), ncol(CC))
for(i in 1:ncol(CC)){
  for(j in 1:ncol(CC)){
    C[ids[i],ids[j]] = CC[i,j]
  }
}

# Z effect of species-specific trait on environmental responses (betas)
Z = matrix(rnorm((n_t + 1) * n_pars, 0 , 0.5),  (n_t + 1), n_pars)

# expected betas
M=TT%*%Z

# define the variation around the betas
Sigma=diag(n_pars)*0.3
rho=0.7 # defines the role of phylogenetic relatedness

betas=rmvnorm(1,mean=as.vector(M),kronecker(Sigma,rho*C+(1-rho)*diag(n_sp)))
# species responses come from a multivariate distribution where the vector 
# m is the species responses mediated by effect of traits and the variance 
# covariance matrix is adjusted to vary according to phylogenetic variation 
# via a kronecker product

Beta=matrix(betas[1,],n_sp,n_pars)

# animals are observed in transects via distance sampling, we simulate 
# observation like so 
n_sites=100

X=cbind(rep(1,n_sites),matrix(rnorm(n_sites*n_env),n_sites,n_env))
# X matrix contains all the environmental data (plus a column of 1s for the intercepts)

B=5 # max sampling distance

# simulate real abundances
N = matrix(NA, n_sites, n_sp)
for(i in 1: n_sp){
  N[,i] = rpois(n_sites, lambda = exp(X %*% Beta[i, ]))
}

p = rbeta(n_sp, 5,3) # probability of detection is unique for every species

# simulate observed dataset
data = NULL
for (i in 1:n_sites) {
  for(j in 1:n_sp){
    if(N[i,j] > 0){
      #d = runif(N[i,j], 0, B)
      #gs = rpois(N[i,j], lambda.group[j]) + 1
      #sigma = exp(Beta[j,1] + gs * Beta[j,2] + TT[j,2] * Beta[j,3])
      #p = exp(-d * d/(2 * (sigma^2)))
      y = rbinom(N[i,j], 1, p[j])
      #d = d[y == 1]
      #gs = gs[y == 1]
      y = y[y == 1]
      if (sum(y) > 0){
        data = rbind(data, cbind(rep(i, sum(y)), rep(j, sum(y)), y))
      }
    }
  }
}



colnames(data) = c("site","sp", "y")
datos = as.data.frame(data)


# run the stan model to estimate parameters

# library(cmdstanr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


stan_dat <- list(
  n_obs = nrow(datos),
  area = rep(1.0, n_sites),
  n_sites = n_sites,
  site = as.integer(datos$site), # (c(s, rep(1:n_sites, each = nzs ))),
  K = ncol(X), 
  X = X,
  n_max = rep(100, n_sp),
  n_s = as.integer(n_sp),
  n_t = dim(TT)[2],
  TT = TT,
  C = C,
  ones = numeric(n_sp) + 1,
  sp = datos$sp,
  p_obs = p
)

# BASE MODEL

fit <- stan(file = 'poisson_binomial_pobs.stan',
            data = stan_dat,
            #init = init_f,
            #pars = pars,
            iter = 1000, thin = 1, chains = 3)
# 578

# 738 gio's pc

# NO POISSON-BIN LOOP

fit_2 <- stan(file = 'jsdm_stan_demo_noloop.stan',
            data = stan_dat,
            #init = init_f,
            #pars = pars,
            iter = 1000, thin = 1, chains = 3)

# 371 seconds gio's pc 

# PETERSON MODEL

init_f <- function () list(beta_std = matrix(0, n_sp, n_pars))

fit_p <- stan(file = 'jsdm_stan_demo_noloop_cholesky_p.stan',
              data = stan_dat,
              init = init_f,
              #pars = pars,
              iter = 1000, thin = 1, chains = 3)


# 1900 s gio's pc

# NHUURRE MODEL

init_f <- function () list(beta_std = matrix(0, n_sp, n_pars))

fit_3<- stan(file = 'jsdm_stan_demo_noloop_cholesky_n.stan',
             data = stan_dat,
             init = init_f,
             #pars = pars,
             iter = 1000, thin = 1, chains = 3)

# 400 seconds gio's pc


pars <- c( "b_m", "rho",  "Sigma", "z", "D")


# some plots
fit_summary<-summary(fit)$summary
fit_summary_2<-summary(fit_2)$summary
fit_summary_3<-summary(fit_3)$summary

# compare r-hat and effective sample size across models

par(mfrow = c(3,2))

hist(fit_summary[,10], main = "R-hat Std Model")
hist(fit_summary[,9], main = "n-eff Std Model" )

hist(fit_summary_2[,10], main = "R-hat Noloop")
hist(fit_summary_2[,9], main = "n-eff Noloop" )

hist(fit_summary_3[,10], main = "R-hat Cholesky Noloop")
hist(fit_summary_3[,9], main = "n-eff Cholesky Noloop" )

# Phylogenetic signals across models 

par(mfrow = c(3,1))

draws <- extract(fit, pars = "rho")
draws_2 <- extract(fit_2, pars = "rho")
draws_3<-extract(fit_3, pars = "rho")


plot(density(draws$rho), main="Rho Std Model")
abline(v=rho, col="red")

plot(density(draws_2$rho), main="Rho Noloop")
abline(v=rho, col="red")

plot(density(draws_3$rho), main="Rho Cholesky Noloop")
abline(v=rho, col="red")


# plot trait level parameters

zs <- fit_summary[grepl("z", rownames(fit_summary)),]
zs_2 <- fit_summary_2[grepl("z", rownames(fit_summary_2)),]
zs_3 <- fit_summary_3[grepl("Z", rownames(fit_summary_3)),] #this model estimates matrix Z directly


#plot(c(Z) - zs[,1])
df<- data.frame(x = 1:dim(zs)[1],
                 tz = c(Z),
                 fz = zs[,1],
                 L = zs[,4],
                 U = zs[,8])


df_2 <- data.frame(x = 1:dim(zs_2)[1],
                 tz = c(Z),
                 fz = zs_2[,1],
                 L = zs_2[,4],
                 U = zs_2[,8])

df_3 <- data.frame(x = 1:dim(zs_3)[1],
                 tz = c(Z),
                 fz = zs_3[,1],
                 L = zs_3[,4],
                 U = zs_3[,8])

par(mfrow = c(3,1))

zplot1<-ggplot(df, aes(x = x, y = tz)) +
  geom_point(size = 3, color="red") +
  geom_point(aes(y = fz), size = 2) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  theme_classic()+
  ggtitle("Trait Parameters Std Model")


zplot2<-ggplot(df_2, aes(x = x, y = tz)) +
  geom_point(size = 3, color="red") +
  geom_point(aes(y = fz), size = 2) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  theme_classic()+
  ggtitle("Trait Parameters Noloop")


zplot3<-ggplot(df_3, aes(x = x, y = tz)) +
  geom_point(size = 3, color="red") +
  geom_point(aes(y = fz), size = 2) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  theme_classic()+
  ggtitle("Trait Parameters Cholesky Noloop")

library(ggpubr)

ggarrange(zplot1,zplot2,zplot3)



# plot intercepts and slopes

# parameters estimates 

bs <- fit_summary[grepl("b_m", rownames(fit_summary)),]
bs_2 <- fit_summary_2[grepl("b_m", rownames(fit_summary_2)),]
bs_3 <- fit_summary_3[grepl("b_m", rownames(fit_summary_3)),]

nf = layout(matrix(c(1,2,3,4,0,0),3,2,byrow=TRUE), widths=c(1,1), heights=c(1,1,0.1))

# STD MODEL BETA PLOTS

par(mfrow = c(2,2))


#layout.show(nf)
op <- par( mar = c(3, 3, 2, 2) + 0.1, mgp = c(3.5, 1, 0), las = 1, bty = "n", cex = 1.2)
plot(scale(dgrass), Beta[,2],  ylab = "", xlab = "", main = "intercept std mod", ylim=c(-3,3))
points(scale(dgrass), bs[seq(1, (n_sp*3) ,by=3) + 1,1], pch = 19, col = 2)

plot(scale(dgrass), Beta[,3], ylab = "", xlab = "", main = "slope std mod", ylim=c(-3,3) )
points(scale(dgrass), bs[seq(1, (n_sp*3) ,by=3) + 2,1], pch = 19, col = 2)

plot(scale(log_bm),Beta[,2], ylab = "", xlab = "", ylim=c(-3,3))
points(scale(log_bm), bs[seq(1, (n_sp*3) ,by=3) + 1,1], pch = 19, col = 2)

plot(scale(log_bm),Beta[,3], ylab = "", xlab = "", ylim=c(-3,3))
points(scale(log_bm),  bs[seq(1, (n_sp*3) ,by=3) + 2,1], pch = 19, col = 2)
mtext("         scaled grass           scaled log body mass", side = 1, line = -1, outer = TRUE, cex=1.3)


# NOLOOP BETA PLOTS 

#layout.show(nf)
op <- par( mar = c(3, 3, 2, 2) + 0.1, mgp = c(3.5, 1, 0), las = 1, bty = "n", cex = 1.2)
plot(scale(dgrass), Beta[,2],  ylab = "", xlab = "", main = "intercept noloop", ylim=c(-3,3))
points(scale(dgrass), bs_2[seq(1, (n_sp*3) ,by=3) + 1,1], pch = 19, col = 2)

plot(scale(dgrass), Beta[,3], ylab = "", xlab = "", main = "slope noloop", ylim=c(-3,3) )
points(scale(dgrass), bs_2[seq(1, (n_sp*3) ,by=3) + 2,1], pch = 19, col = 2)

plot(scale(log_bm),Beta[,2], ylab = "", xlab = "", ylim=c(-3,3))
points(scale(log_bm), bs_2[seq(1, (n_sp*3) ,by=3) + 1,1], pch = 19, col = 2)

plot(scale(log_bm),Beta[,3], ylab = "", xlab = "", ylim=c(-3,3))
points(scale(log_bm),  bs_2[seq(1, (n_sp*3) ,by=3) + 2,1], pch = 19, col = 2)
mtext("         scaled grass           scaled log body mass", side = 1, line = -1, outer = TRUE, cex=1.3)


# NOLOOP CHOLESKY BETA PLOTS 
op <- par( mar = c(3, 3, 2, 2) + 0.1, mgp = c(3.5, 1, 0), las = 1, bty = "n", cex = 1.2)
plot(scale(dgrass), Beta[,2],  ylab = "", xlab = "", main = "intercept noloop chol.", ylim=c(-3,3))
points(scale(dgrass), bs_3[seq(1, (n_sp*3) ,by=3) + 1,1], pch = 19, col = 2)

plot(scale(dgrass), Beta[,3], ylab = "", xlab = "", main = "slope noloop chol.", ylim=c(-3,3) )
points(scale(dgrass), bs_3[seq(1, (n_sp*3) ,by=3) + 2,1], pch = 19, col = 2)

plot(scale(log_bm),Beta[,2], ylab = "", xlab = "", ylim=c(-3,3))
points(scale(log_bm), bs_3[seq(1, (n_sp*3) ,by=3) + 1,1], pch = 19, col = 2)

plot(scale(log_bm),Beta[,3], ylab = "", xlab = "", ylim=c(-3,3))
points(scale(log_bm),  bs_3[seq(1, (n_sp*3) ,by=3) + 2,1], pch = 19, col = 2)
mtext("         scaled grass           scaled log body mass", side = 1, line = -1, outer = TRUE, cex=1.3)





par(mfrow = c(2,2))

# BETAS POSTERIOR VS TRUE VALUES 

plot(c(Beta), 
     c(bs[seq(1, (n_sp*3) ,by=3) ,1],
       bs[seq(1, (n_sp*3) ,by=3) + 1,1],
       bs[seq(1, (n_sp*3) ,by=3) + 2,1]), 
     xlab = "true value", 
     ylab = "posterior mean",
     main="Posterior Mean ~ True Values Std Model")
abline(0,1, col="red")

plot(c(Beta), 
     c(bs_2[seq(1, (n_sp*3) ,by=3) ,1],
       bs_2[seq(1, (n_sp*3) ,by=3) + 1,1],
       bs_2[seq(1, (n_sp*3) ,by=3) + 2,1]), 
     xlab = "true value", 
     ylab = "posterior mean",
     main="Posterior Mean ~ True Values Noloop Model")
abline(0,1, col="red")

plot(c(Beta), 
     c(bs_3[seq(1, (n_sp*3) ,by=3) ,1],
       bs_3[seq(1, (n_sp*3) ,by=3) + 1,1],
       bs_3[seq(1, (n_sp*3) ,by=3) + 2,1]), 
     xlab = "true value", 
     ylab = "posterior mean",
     main="Posterior Mean ~ True Values Noloop-Cholesky Model")
abline(0,1, col="red")




# yet another one
df <- data.frame(x = 1:dim(bs)[1],
                 tb = c(Beta),
                 fb = c(bs[seq(1, (n_sp*3) ,by=3) ,1], bs[seq(1, (n_sp*3) ,by=3) + 1,1], bs[seq(1, (n_sp*3) ,by=3) + 2,1]),
                 L =  c(bs[seq(1, (n_sp*3) ,by=3) ,4], bs[seq(1, (n_sp*3) ,by=3) + 1,4], bs[seq(1, (n_sp*3) ,by=3) + 2,4]),
                 U =  c(bs[seq(1, (n_sp*3) ,by=3) ,8], bs[seq(1, (n_sp*3) ,by=3) + 1,8], bs[seq(1, (n_sp*3) ,by=3) + 2,8])
)

df_2<- data.frame(x = 1:dim(bs_2)[1],
                  tb = c(Beta),
                  fb = c(bs_2[seq(1, (n_sp*3) ,by=3) ,1], bs_2[seq(1, (n_sp*3) ,by=3) + 1,1], bs_2[seq(1, (n_sp*3) ,by=3) + 2,1]),
                  L =  c(bs_2[seq(1, (n_sp*3) ,by=3) ,4], bs_2[seq(1, (n_sp*3) ,by=3) + 1,4], bs_2[seq(1, (n_sp*3) ,by=3) + 2,4]),
                  U =  c(bs_2[seq(1, (n_sp*3) ,by=3) ,8], bs_2[seq(1, (n_sp*3) ,by=3) + 1,8], bs_2[seq(1, (n_sp*3) ,by=3) + 2,8])
)

df_3<- data.frame(x = 1:dim(bs_3)[1],
                  tb = c(Beta),
                  fb = c(bs_3[seq(1, (n_sp*3) ,by=3) ,1], bs_3[seq(1, (n_sp*3) ,by=3) + 1,1], bs_3[seq(1, (n_sp*3) ,by=3) + 2,1]),
                  L =  c(bs_3[seq(1, (n_sp*3) ,by=3) ,4], bs_3[seq(1, (n_sp*3) ,by=3) + 1,4], bs_3[seq(1, (n_sp*3) ,by=3) + 2,4]),
                  U =  c(bs_3[seq(1, (n_sp*3) ,by=3) ,8], bs_3[seq(1, (n_sp*3) ,by=3) + 1,8], bs_3[seq(1, (n_sp*3) ,by=3) + 2,8])
)

# true vs fitted betas with credible intervals 

tbfb1<-ggplot(df, aes(x = x, y = tb)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fb), size = 1) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  ggtitle("True vs Estimated ± CI (Std Mod)")+
  theme_classic()



tbfb2<-ggplot(df_2, aes(x = x, y = tb)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fb), size = 1) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  ggtitle("True vs Estimated ± CI (Noloop Mod)")+
  theme_classic()



tbfb3<-ggplot(df_3, aes(x = x, y = tb)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fb), size = 1) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  ggtitle("True vs Estimated ± CI (Noloop-Cholesky Mod)")+
  theme_classic()



ggarrange(tbfb1,tbfb2,tbfb3)



# plot density estimates
D <- fit_summary[grepl("D", rownames(fit_summary)),]
D2 <- fit_summary_2[grepl("D", rownames(fit_summary_2)),]
D3 <- fit_summary_3[grepl("D", rownames(fit_summary_3)),]




df <- data.frame(x = 1:dim(D)[1],
                 td = colSums(N)/n_sites,
                 fd = D[,1],
                 L =  D[,4],
                 U =  D[,8]
)


df_2 <- data.frame(x = 1:dim(D2)[1],
                 td = colSums(N)/n_sites,
                 fd = D2[,1],
                 L =  D2[,4],
                 U =  D2[,8]
)

df_3 <- data.frame(x = 1:dim(D3)[1],
                 td = colSums(N)/n_sites,
                 fd = D3[,1],
                 L =  D3[,4],
                 U =  D3[,8]
)



dplot<-ggplot(df, aes(x = x, y = td)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fd), size = 1) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  ggtitle("Posterior density Std Model")+
  theme_classic()

dplot2<-ggplot(df_2, aes(x = x, y = td)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fd), size = 1) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  ggtitle("Posterior density Noloop Model")+
  theme_classic()

dplot3<-ggplot(df_3, aes(x = x, y = td)) +
  geom_point(size = 2, color="red") +
  geom_point(aes(y = fd), size = 1) +
  geom_linerange(aes(ymin = L, ymax = U)) +
  ggtitle("Posterior density Noloop-Cholesky Model")+
  theme_classic()

ggarrange(dplot,dplot2,dplot3)
