library(rstan)

# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(cmdstanr)
# install_cmdstan(cores = 2)

# set_cmdstan_path("/Users/juanmanuelmorales/cmdstan")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# prepare matrices
# 

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

# for(i in 2:ncol(X)){
#   X[,i] <- scale(Xo[ ,i])
# }

# NDVI & Rain

rain <- as.matrix(read.csv("rain_matrix.csv"))
grass <- log(as.matrix(read.csv("grass_matrix.csv")))
ndvi <- as.matrix(read.csv("ndvi_matrix.csv"))
a_ndvi <- as.matrix(read.csv("anomaly_ndvi_matrix.csv"))
d_ndvi <- as.matrix(read.csv("Delta_NDVI.csv"))


dist_mat <- as.matrix(read.csv("dist_matrix.csv"))

# scale ndvi and rain
for(i in 1:ncol(ndvi)){
  rain[,i]<-(rain[,i]-mean(rain))/sd(rain)
  grass[,i]<-(grass[,i]-mean(grass))/sd(grass)
  ndvi[,i]<-(ndvi[,i]-mean(ndvi))/sd(ndvi)
  a_ndvi[,i]<-(a_ndvi[,i]-mean(a_ndvi))/sd(a_ndvi)
  d_ndvi[,i]<-(d_ndvi[,i]-mean(d_ndvi))/sd(d_ndvi)
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

stan.data <- list(
  n_obs = nrow(obs),
  n_dates = n_dates,
  n_tcov = 5,
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
  p_obs = p,
  Dmat = dist_mat[,-1]/10000
  #g_size = obs$Total
)

pars <- c( "b_m", "rho",  "Sigma", "z", "Z")


init_f <- function () list(b_m = matrix(0, n_sp, n_pars+5))


modelName <- "poisson_binomial_dates_pobs_se"
# modelName <- "hmm-2-state_classic"
if (TRUE) {
  model_file <- file.path(paste(getwd(),"/", modelName, ".stan", sep = ""))
  mod <- cmdstan_model(model_file)
}

fit <- mod$sample(
  data = stan.data, 
  init = init_f,
  chains = 3, 
  parallel_chains = 3,
  iter_warmup = 1000, 
  iter_sampling = 1000, 
  thin = 1, 
  seed = 123)




# fit <- stan(file = 'poisson_binomial_dates_pobs.stan',
#             data = stan.data,
#             init = init_f,
#             pars = pars,
#             iter = 100, thin = 1, chains = 1)


Y = matrix(NA, n_sites, n_sp)
for(i in 1:n_sites){
  for(j in 1:n_sp){
    tmp = which(obs$Animal.sp == j & obs$poly.id == i)
    Y[i,j] = length(tmp)
  }
}

fit_summary <- summary(fit)$summary

# check R_hat
hist(fit_summary[,10], 100)

# check n_eff
hist(fit_summary[,9], 100)
# tmp <- which(fit_summary[,9] < 100)

print(fit, pars = c( "rho"))

# plot of regression lines by spp
mean_b <- fit_summary[grepl("b_m", rownames(fit_summary)),]
mean_z <- fit_summary[grepl("Z", rownames(fit_summary)),]

# matrix with regression coefs in rows by spp in cols
mb = matrix(mean_b[,1], n_pars+2, n_sp)
# matrix with coefs for traits in rows and covariates in cols
mz = matrix(mean_z[,1], ncol(Traits), n_pars+2, byrow = TRUE)

op = par(mfrow= c(4,3), las = 1)
for(j in 2:dim(X)[2]){
  #f = 1 # (j-1)*n_sp - n_sp + 1
  plot(Xo[,j], exp(mb[1,1]+mb[j,1]*X[,j]), type = "l", ylim = c(0,2),
       xlab = colnames(X)[j], ylab = "lambda")
  for(i in 2:n_sp) lines(Xo[,j], exp(mb[1,i]+mb[j,i]*X[,j]))
}
par(op)

# find coefficients different from zero
tmp = which(sign(mean_z[,1]) == sign(mean_z[,4]) & sign(mean_z[,1])==sign(mean_z[,8]))
mean_z[tmp,c(1,4,8:10)]

tmp = which(sign(mean_b[,1]) == sign(mean_b[,4]) & sign(mean_b[,1])==sign(mean_b[,8]))
mean_b[tmp,c(1,4,8:10)]



plot(X[,3], mean_b[1,1]+mean_b[1,2]*X[,2], type = "l", ylim = c(-1,1))
for(i in 2:n_sp) lines(X[,2], mean_b[i,1]+mean_b[i,2]*X[,2])

zs <- extract(fit, pars = "Z")

f = matrix(NA, 5,13)
for(i in 1:5){
  for(j in 1:13){
    f[i,j] = length(which(sign(zs$Z[,i,j]) == sign(mean(zs$Z[,i,j]))))/length(zs$Z[,i,j])
  }
}

mean_zs <- fit_summary[grepl("z", rownames(fit_summary)),]

library(coda)
nsim <- 1500
xp <- seq(min(Traits[,2]), max(Traits[,2]), length.out = 101)
xg <- seq(min(Traits[,3]), max(Traits[,3]), length.out = 100)
z1 <- matrix(NA, nsim, length(xp))
z2 <- matrix(NA, nsim, length(xp))
z3 <- matrix(NA, nsim, length(xg))
z4 <- matrix(NA, nsim, length(xg))


# coeffs and body mass

#tiff(filename="zetas.tif", res=res, unit="cm", height=17, width=17)
nf = layout(matrix(1:16,4,4,byrow=TRUE), widths=c(1,1), heights=rep(1,16))
#layout.show(nf)
op <- par( mar = c(3, 3, 2, 2) + 0.1, las = 1, bty = "n", cex = 0.8)

plot(TT[,4] , mb[5,], pch = 19, col=rgb(0,0,0, 0.4), 
     ylab = "", xlab = "", main = "",  log = "x")
lines(TT[,4], mz[1,5] + mz[2,5]*Traits[,2], lwd = 3)


for(i in 1:nsim) z1[i,] <- zs$Z[i,1,1] + zs$Z[i,2,1]*xp
ci1 <- HPDinterval(as.mcmc(z1)) 
lines(seq(min(TT[,4]), max(TT[,4]), length.out = 101), ci1[,1])
lines(seq(min(TT[,4]), max(TT[,4]), length.out = 101), ci1[,2])


plot(TT[,5] , mb[5,], pch = 19, col=rgb(0,0,0, 0.4), 
     ylab = "", xlab = "", main = "")
lines(TT[,5], mz[1,2] + mz[3,5]*Traits[,3], lwd = 3)
for(i in 1:nsim) z1[i,] <- zs$Z[i,1,1] + zs$Z[i,2,1]*xp
ci1 <- HPDinterval(as.mcmc(z1)) 
lines(seq(min(TT[,4]), max(TT[,4]), length.out = 101), ci1[,1])
lines(seq(min(TT[,4]), max(TT[,4]), length.out = 101), ci1[,2])



plot(TT[,4] , mb[1,], pch = 19, col=rgb(0,0,0, 0.4), 
     ylab = "", xlab = "", main = "",  log = "x")
lines(TT[,4], mz[1,1] + mz[2,1]*Traits[,2], lwd = 3)
for(i in 1:nsim) z1[i,] <- zs$Z[i,1,1] + zs$Z[i,1,2]*xp
ci1 <- HPDinterval(as.mcmc(z1)) 
lines(seq(min(TT[,3]), max(TT[,3]), length.out = 101), ci1[,1])
lines(seq(min(TT[,3]), max(TT[,3]), length.out = 101), ci1[,2])



mean_b <- fit_summary[grepl("b_m", rownames(fit_summary)),]
mean_z <- fit_summary[grepl("Z", rownames(fit_summary)),]
# plot coeffs against body mass
op = par(mfrow = c(2,5), mar = c(3, 3, 2, 2) + 0.1, las = 1, bty = "n")
for(i in 1:10){
  idx = seq()
  plot(TT[,3] , mean_b[seq(1, (n_sp*n_pars), by=n_pars) + i -1, 1], 
       pch = 19, col=rgb(0,0,0, 0.4), 
       ylab = "", main = colnames(X)[i], xlab = "", log = "x")
  lines(TT[,3], mean_z[i,1] + mean_z[i+n_pars,1]*Traits[,2], lwd = 3)
}
par(op)

# now against water dependence
op = par(mfrow = c(2,5), mar = c(3, 3, 2, 2) + 0.1, las = 1, bty = "n")
for(i in 1:10){
  idx = seq()
  plot(TT[,4] , mean_b[seq(1, (n_sp*n_pars), by=n_pars) + i -1, 1], 
       pch = 19, col=rgb(0,0,0, 0.4), 
       ylab = "", main = colnames(X)[i], xlab = "")
  lines(TT[,4], mean_z[i,1] + mean_z[i+n_pars*2,1]*Traits[,3], lwd = 3)
}
par(op)




plot(TT[,3], mean_b[seq(1, (n_sp*n_pars) ,by=n_pars) + 1, 1], pch = 19, col=rgb(0,0,0, 0.4),
     ylab = "", xlab = "", main = "", log="x")

plot(TT[,3], mean_b[seq(1, (n_sp*n_pars) ,by=n_pars) + 2, 1], pch = 19, col=rgb(0,0,0, 0.4),
     ylab = "", xlab = "", main = "", log="x")



lines(TT[,3], mean_z[2,1] + mean_z[4,1]*Traits[,2], lwd = 3)
for(i in 1:nsim) z2[i,] <- zs$z[i,4] + zs$z[i,5]*xp
ci2 <- HPDinterval(as.mcmc(z2)) 
lines(0:100, ci2[,1])
lines(0:100, ci2[,2])

#for(i in 1:100) lines(fd, zs$z[i,4] + zs$z[i,5]*TT[,2], col=rgb(0,0,0, 0.2) )
plot(TT[,4], mean_b[seq(1, (n_sp*2),by=2),1], pch = 19, col=rgb(0,0,0, 0.4),
     ylab = "", xlab = "", main = "")
lines(bs, mean_z[1,1] + mean_z[5,1]*TT[,3], lwd = 3)
for(i in 1:nsim) z3[i,] <- zs$z[i,1] + zs$z[i,3]*xg
ci3 <- HPDinterval(as.mcmc(z3)) 
lines(seq(min(bs), max(bs), length.out = 100), ci3[,1])
lines(seq(min(bs), max(bs), length.out = 100), ci3[,2])


plot(TT[,4], mean_b[seq(1, (n_sp*2) ,by=2) + 1, 1], pch = 19, col=rgb(0,0,0, 0.4),
     ylab = "", xlab = "",  main = "")

for(i in 1:nsim) z4[i,] <- zs$z[i,4] + zs$z[i,6]*xg
ci4 <- HPDinterval(as.mcmc(z4)) 
#lines(seq(min(bs), max(bs), length.out = 100), ci4[,1])
#lines(seq(min(bs), max(bs), length.out = 100), ci4[,2])

#curve(mean_z[2,1] + mean_z[6,1]*x, add = TRUE, lwd = 3)
#mtext("log body mass", 
#      side = 3, line = -31.5, outer = TRUE, cex=1.2)
#mtext("frugivory", 
#      side = 3, line = -17.5, outer = TRUE, cex=1.2)
#mtext("intercept", side = 2, line = -1, las = 0, outer = TRUE, cex=1.2)
par(op)
par(mfrow=c(1,1))
dev.off()

#----
# prepare tables
gens <- sort(unique(datos$Bird_genus))
slopes <- mean_b[seq(2, (n_sp*2)+1 ,by=2),]
rownames(slopes) = gens
write.csv(slopes, file = "slopes.csv")

intercepts <-  mean_b[seq(1, (n_sp*2),by=2),]
rownames(intercepts) = gens
write.csv(intercepts, file = "intercepts.csv")

write.csv(mean_z, file = "zetas.csv")

eta <- fit_summary[grepl("etasq", rownames(fit_summary)),]
delta <- fit_summary[grepl("delta", rownames(fit_summary)),]
rhos <- fit_summary[grepl("rhosq", rownames(fit_summary)),]

write.csv(rbind(eta, delta, rhos), file = "etas.csv") 

rho <- fit_summary[grepl("rho", rownames(fit_summary)),]

sigma <- fit_summary[grepl("Sigma", rownames(fit_summary)),]
write.csv(rbind(sigma, rho[1,]), file = "sigma.csv") 

#-----

# histogram of slopes  
tiff(filename="hist_slopes.tif", res=res, unit="cm", height=8, width=8)
op <- par(las = 1, bty = "n", cex = 1)
hist(mean_b[seq(2, (n_sp*2)+1 ,by=2),1], 60, 
     main = "", xlab = "slope", col=rgb(0,0,0, 0.5))
par(op)
dev.off()

gens <- sort(unique(datos$Bird_genus))
slopes <- mean_b[seq(2, 310, by = 2), 1]
toplot = gens[which(slopes < 0)]

plot(lipid, rep(-10, length(lipid)), ylim=c(0,100), ylab = "visits", xlab = "")
for(i in 1: length(toplot)){
  tmp <- which(datos$Bird_genus == toplot[i])
  idx <- which(gens ==  toplot[i])
  points(lipid[tmp], Y[tmp], pch = 19, col=rgb(0,0,0, 0.4))
  curve(exp( mean_b[2*idx-1, 1]  + mean_b[2*idx, 1] * x), add = TRUE, lwd = 2)
}



print(fit, pars = c("rhosq", "rhosqp", "etasq", "etasqp", "delta", "deltap"))

plot(fit, plotfun = "rhat")
plot(fit, pars = list("Z"), 
     show_density = TRUE, ci_level = 0.5, fill_color = "gray")


# plot response against lipids -------------------------------------------------

mz <- matrix(mean_z[,1], 3, 2, byrow = TRUE)

xl <- seq(-0.9, 2.4, by = 0.1 )
fro <- seq(0, 100, by = 20)
fr <- (fro - mean(fd))/ sd(fd)

xlt <- xl * sd(datos$lipid) + mean(datos$lipid)
# grs <- viridis(length(fr))

rwb <- colorRampPalette(colors = c("grey", "blue"))
grs <- rwb(length(fr))
#ev <- matrix(NA, length(xl), length(fr))

tiff(filename="lipids.tif", res=res, unit="cm", height=10, width=11)
op <- par(las = 1, bty = "n", cex = 1)
plot(xlt, rep(-1, length(xl)), ylim = c(0,8),
     xlab = "lipids", ylab = "expected visits")
legend("bottomright", legend = paste( fro, "%"), 
       col = grs, pch = 19, bty = "n", cex = 0.6)
for(i in 1:length(fr)){
  lines(xlt, exp( mz[1,1] + mz[2,1] * fr[i]  + (mz[2,1] + mz[2,2] *fr[i]) * xl), 
        lwd = 3, col = grs[i])
  #ev[,i] = exp( mz[1,1] + mz[2,1] * fr[i]  + (mz[2,1] + mz[2,2] *fr[i]) * xl)
}
par(op)
dev.off()
#--------------

#  plot by gens ---------------------------------------------------------------

#shadesOfGrey <- colorRampPalette(c("grey0", "grey100"))
#grs <- shadesOfGrey(100)
grs <- viridis(100)

gens <- sort(unique(datos$Bird_genus))

toplot <- c("Trichothraupis", "Tachyphonus", "Ramphocelus", "Dacnis","Tersina", "Tangara")

ids <- NULL
for(i in 1:length(toplot)){
  ids <- c(ids, round(fd[which(gens ==  toplot[i])]))
}

sl <- sort(datos$lipid)
ssl <- sort(lipid)

op <- par(mfrow = c(1,2), las = 1, bty = "n", cex = 1.2, mar = c(4, 3, 2, 2) + 0.1)

plot(datos$lipid, rep(-10, length(lipid)), ylim=c(0,60), ylab = "visits", xlab = "lipid")
legend("topleft", legend = paste(toplot, ", frug = ", ids), col = grs[ids], pch = 19, bty = "n", cex = 0.9)
for(i in 1: length(toplot)){
  tmp <- which(datos$Bird_genus == toplot[i])
  idx <- which(gens ==  toplot[i])
  #points(datos$lipid[tmp], Y[tmp], pch = 19, col=rgb(0,0,0, 0.4))
  lines(sl, exp( mean_b[2*idx-1, 1]  + mean_b[2*idx, 1] * ssl), lwd = 3, col = grs[round(fd[idx])])
}

toplot <- c("Phyllomyias","Tyrannus", "Elaenia","Pitangus", "Empidonomus")

ids <- NULL
for(i in 1:length(toplot)){
  ids <- c(ids, ceiling(fd[which(gens ==  toplot[i])]))
}
# Empidonomus 40
# Myiarchus 30
# Myiodynastes 20
# Tyrannus 10
# Pitangus 30 
# Phyllomyias 0
# Elaenia 20

plot(datos$lipid, rep(-10, length(lipid)), ylim=c(0,60), ylab = "visits", xlab = "lipid")
legend("topleft", legend = paste(toplot, ", frug = ", ids), col = grs[ids], pch = 19, bty = "n", cex = 0.9)

for(i in 1: length(toplot)){
  tmp <- which(datos$Bird_genus == toplot[i])
  idx <- which(gens ==  toplot[i])
  #points(datos$lipid[tmp], Y[tmp], pch = 19, col=rgb(0,0,0, 0.4))
  lines(sl, exp( mean_b[2*idx-1, 1]  + mean_b[2*idx, 1] * ssl), lwd = 3, col = grs[ceiling(fd[idx]+0.1)])
}

par(op)

#-------------------------------------------------------------------------------

# sum of residuals from data and predicted data
sr <- extract(fit, pars = "sr")
sr_rep <- extract(fit, pars = "sr_rep")
plot(sr$sr, sr_rep$sr_rep)
abline(0,1)


# lipids vs counts

library(coda)


# residuals a la DHARMa -------------------------------------------------------
y_rep <- extract(fit, pars = "y_pred")
err <- numeric(length(Y))
for(i in 1:length(Y)){
  ec <- ecdf(y_rep$y_pred[,i])
  err[i] <- ec(Y[i])
}
hist(err, 100) # we would like these to be uniformly distributed

plot(Y, colMeans(y_rep$y_pred), xlim = c(0,500), ylim=c(0,500))
abline(0,1)
# ---------



y_rep <- extract(fit, pars = "y_pred")
ul <- unique(lipid)
frug <- datos$Fruit_Diet

uf <- unique(datos$Fruit_Diet)

ff <- seq(0, 100, by = 10)
ll <- seq(0, 100, by = 10)

mm <- matrix(NA, length(ll), length(ff))
mr <- matrix(NA, length(ll), length(ff))
for(i in 1:(length(ll) -1)){
  for(j in (1:length(ff) -1)){
    tmp <- which(datos$lipid >= ll[i] & datos$lipid < ll[i+1] & datos$Fruit_Diet >= ff[j] & datos$Fruit_Diet < ff[j+1])
    if(! is.null(tmp)){
      mm[i,j] <- mean(Y[tmp])
      mr[i,j] <- mean(y_rep$y_pred[,tmp])
    } 
  }
}

plot(mm,mr)
abline(0,1)

ym <- numeric(length(ul))

for (i in 1: length(ul)){
  tmp <- which(lipid == ul[i])
  ym[i] <- mean(Y[tmp])
  #HPDinterval(as.mcmc(y_rep$y_pred[ , tmp]))
  
}


eta <- fit_summary[grepl("etasq", rownames(fit_summary)),]
delta <- fit_summary[grepl("delta", rownames(fit_summary)),]
rhos <- fit_summary[grepl("rhosq", rownames(fit_summary)),]

hist(Dmat, 100, freq = FALSE, main = "", xlab = "distance (thousand km)", ylim = c(0,1.2))
xtmp <- seq(0, 10, by = 0.1)
lines(xtmp, eta[1,1]*exp(-rhos[1,1]*xtmp^2), lwd = 3)
#lines(xtmp, eta[1,1]*exp(-rhos[1,1]*xtmp^2)+delta[1,1], lwd = 3)

hist(DistP, 30, freq = FALSE, main = "", xlab = "cophenetic distance", ylim = c(0,0.1))
xtmp <- seq(0, 400, by = 1)
lines(xtmp, eta[2,1]*exp(-rhos[2,1]*(xtmp/100)^2), lwd = 3)
#lines(xtmp, eta[2,1]*exp(-rhos[2,1]*(xtmp/100)^2)+delta[2,1], lwd = 3)


# save(fit, file = "fit_pois.RData")
#, control = list(adapt_delta = 0.9, max_treedepth = 15))

#plot(fit, pars = "rho")


