# variables for plots

#  13 environmental vars + intercept 

Xo <- as.matrix(read.csv("covariates_final.csv"))
varnames<-c(colnames(Xo),"rain","grass","ndvi","a_ndvi","d_ndvi")

#traits
TT<-read.csv("TT_matrix.csv")
names(TT)

# select coefs for each of the 5 traits (Intercept (Grazer), Body Mass, Water Dependence,
# Browser Correction, Mixed feeders Corr)

#bodymass

bodymass_coefs<-as.data.frame(traits_effects[15:28,])
waterdep_coefs<-as.data.frame(traits_effects[29:42,])
browser_coefs<-as.data.frame(traits_effects[43:56,])
mixedfeed_coefs<-as.data.frame(traits_effects[57:70,])




# select different coefs for environmental variables 

intercepts<-as.data.frame(responses[seq(length.out=19, from=1, by=14),])
village_coefs<-as.data.frame(responses[seq(length.out=19, from=2, by=14),])
tourism_coefs<-as.data.frame(responses[seq(length.out=19, from=3, by=14),])
drainages_coefs<-as.data.frame(responses[seq(length.out=19, from=4, by=14),])
woodycover_coefs<-as.data.frame(responses[seq(length.out=19, from=5, by=14),])
nitrogen_coefs<-as.data.frame(responses[seq(length.out=19, from=6, by=14),])
runoff_coefs<-as.data.frame(responses[seq(length.out=19, from=7, by=14),])
fire_coefs<-as.data.frame(responses[seq(length.out=19, from=8, by=14),])
water_coefs<-as.data.frame(responses[seq(length.out=19, from=9, by=14),])
rain_coefs<-as.data.frame(responses[seq(length.out=19, from=10, by=14),])
grass_coefs<-as.data.frame(responses[seq(length.out=19, from=11, by=14),])
ndvi_coefs<-as.data.frame(responses[seq(length.out=19, from=12, by=14),])
andvi_coefs<-as.data.frame(responses[seq(length.out=19, from=13, by=14),])
dndvi_coefs<-as.data.frame(responses[seq(length.out=19, from=14, by=14),])

sp.names<-as.vector(TT[,2])


# Intercept Plot

ggplot(intercepts, aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(intercepts[,4]),
                ymax=as.numeric(intercepts[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Intercepts")+
  ylab("Parameter estimate")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# Village Distance Plot

ggplot(village_coefs, aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(village_coefs[,4]),
                ymax=as.numeric(village_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to Distance from Villages")+
  ylab("Parameter estimate")+
  ylim(-1,1.5)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



# Tourism Plot 

ggplot(tourism_coefs, aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(tourism_coefs[,4]),
                ymax=as.numeric(tourism_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to Tourism")+
  ylab("Parameter estimate")+
  ylim(-1,1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# Drainages plot

ggplot(drainages_coefs, aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(drainages_coefs[,4]),
                ymax=as.numeric(drainages_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to distance from drainages")+
  ylab("Parameter estimate")+
  ylim(-1,1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#woodycover plot

ggplot(woodycover_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(woodycover_coefs[,4]),
                ymax=as.numeric(woodycover_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to density of tree cover")+
  ylab("Parameter estimate")+
  ylim(-1,1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#nitrogen plot


ggplot(nitrogen_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(nitrogen_coefs[,4]),
                ymax=as.numeric(nitrogen_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to distance from N hotspot")+
  ylab("Parameter estimate")+
  ylim(-1,1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# runoff plot

ggplot(runoff_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(runoff_coefs[,4]),
                ymax=as.numeric(runoff_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to water runoff")+
  ylab("Parameter estimate")+
  ylim(-1,1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# fire plot

ggplot(fire_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(fire_coefs[,4]),
                ymax=as.numeric(fire_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to frequency of fires")+
  ylab("Parameter estimate")+
  ylim(-1,1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# water plot

ggplot(water_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(water_coefs[,4]),
                ymax=as.numeric(water_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to distance from permanent water sources")+
  ylab("Parameter estimate")+
  ylim(-4,4)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# rain plot

ggplot(rain_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(rain_coefs[,4]),
                ymax=as.numeric(rain_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to precipitations levels")+
  ylab("Parameter estimate")+
  ylim(-1,1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# grass abundance


ggplot(grass_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(grass_coefs[,4]),
                ymax=as.numeric(grass_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to grass height")+
  ylab("Parameter estimate")+
  ylim(-1,1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# ndvi plot

ggplot(ndvi_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(ndvi_coefs[,4]),
                ymax=as.numeric(ndvi_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to NDVI")+
  ylab("Parameter estimate")+
  ylim(-100,100)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# andvi plot

ggplot(andvi_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(andvi_coefs[,4]),
                ymax=as.numeric(andvi_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to NDVI")+
  ylab("Parameter estimate")+
  ylim(-100,100)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# dndvi plot

ggplot(dndvi_coefs,aes(x=sp.names,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(dndvi_coefs[,4]),
                ymax=as.numeric(dndvi_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Herbivore Species")+
  ggtitle("Species Responses to NDVI")+
  ylab("Parameter estimate")+
  ylim(-100,100)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+

    theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# TRAITS PLOTS 

#bodymass plot
ggplot(bodymass_coefs,aes(x=varnames,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(bodymass_coefs[,4]),
                ymax=as.numeric(bodymass_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Environmental Covariate")+
  ggtitle("Effects of bodymass on responses to landscape")+
  ylab("Parameter estimate")+
  ylim(-1,1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# water dependence plot 

ggplot(waterdep_coefs,aes(x=varnames,y=mean))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(waterdep_coefs[,4]),
                ymax=as.numeric(waterdep_coefs[,8]),width=0.4, size=0.8)+
  theme_bw()+
  xlab("Environmental Covariate")+
  ggtitle("Effects of water dependence on responses to landscape")+
  ylab("Parameter estimate")+
  ylim(-2,2)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# 
