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
  ndvi[,i]<-scale(ndvi[,i])
}

for(i in 1:ncol(rain)){
  rain[,i]<-scale(rain[,i])
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

# scale values and log mass
TT$bodymass<-scale(log(TT$bodymass))
TT$water_dependence<-scale(TT$water_dependence)

# recode guild as numeric
# 0 is grazer
# 1 is browser
# 2 is mixed
TT$guild<-c(0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 2, 2, 0, 2, 1, 0, 2, 2)

TT<-as.matrix(TT)
TT

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
n_sites<-dim(X)[1]
n_sp<-dim(C)[1]
n_pars<-ncol(X)





