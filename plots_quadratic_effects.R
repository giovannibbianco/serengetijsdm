# plot quadratic effects single species models

library(dotwhisker)
library(dplyr)
library(tidyr)

sp.names<-c("wildebeest","zebra","thomson","grant","topi","kongoni","impala","warthog","buffalo","giraffe","elephant","eland","dikdik")

parnames<-c("intercept","village_dist","tourism_dist","river_dist","woodycover","nitrogen_dist","runoff","fire","slope","sodium_dist","water_dist","water_quad","runoff_quad","nitrogen_quad","sodium_quad","ndvi","rain")

# mean of parameters is first column of mean_b summaries
# bind them all together to get df with all mean values 

par_df<-cbind(parnames,mean_b[,1],mean_b2[,1],mean_b3[,1],mean_b4[,1],mean_b5[,1],mean_b6[,1],mean_b7[,1],mean_b9[,1],mean_b10[,1],mean_b11[,1],mean_b12[,1],mean_b13[,1],mean_b16[,1])
colnames(par_df)[2:14]<-sp.names

par_df<-data.frame(par_df)
write.csv(par_df,file="regression_pars_quadratic.csv")

CI_low<-cbind(parnames,mean_b[,4],mean_b2[,4],mean_b3[,4],mean_b4[,4],mean_b5[,4],mean_b6[,4],mean_b7[,4],mean_b9[,4],mean_b10[,4],mean_b11[,4],mean_b12[,4],mean_b13[,4],mean_b16[,4])
colnames(CI_low)[2:14]<-sp.names
CI_low<-data.frame(CI_low)
write.csv(CI_low,file = "CI_low_quadratic.csv")


CI_hi<-cbind(parnames,mean_b[,8],mean_b2[,8],mean_b3[,8],mean_b4[,8],mean_b5[,8],mean_b6[,8],mean_b7[,8],mean_b9[,8],mean_b10[,8],mean_b11[,8],mean_b12[,8],mean_b13[,8],mean_b16[,8])
colnames(CI_hi)[2:14]<-sp.names
CI_hi<-data.frame(CI_hi)
write.csv(CI_hi,file = "CI_hi_quadratic.csv")

# plot parameters by variable for all species by row
# add credible intervals with geom_errorbar 

library(ggplot2)
library(dplyr)
library(forcats)

# load regression parameters data frame

par_df<-read.csv("regression_pars_quadratic.csv")


str(par_df)

# transpose it so that every column is a variable and every row is a species

par_df_t <- data.frame(t(par_df[-1]))

species.name<-row.names(par_df_t)

par_df_t<-cbind(species.name,par_df_t)

row.names(par_df_t)<-NULL

colnames(par_df_t)<-c("species.name",parnames)


par_df_t$species.name<-as.factor(par_df_t$species.name)

par_df_t<-as.data.frame(par_df_t)

str(par_df_t)

# add grazing guild column
par_df_t$species.name
guild<-as.factor(c("G","G","G","B","G","G","M","G","G","B","M","M","B"))

par_df_t<-cbind(par_df_t,guild)

# Quadratic Effect of Distance to Nitrogen Hotspots By Species

nitrogenplot<-par_df_t %>%
  mutate(species.name=fct_relevel(species.name,"dikdik","thomson","grant","impala",
                                  "topi","kongoni","wildebeest","giraffe","eland","buffalo",
                                  "warthog","zebra","elephant"))%>%
  ggplot(aes(x=species.name,y=nitrogen_quad,color=guild))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(CI_low[14,][2:14]),ymax=as.numeric(CI_hi[14,][2:14]),width=0.4, size=0.8)+
  theme_bw()+
  scale_x_discrete(labels=c("Dikdik","Thomson's Gazelle","Grant's Gazelle", "Impala",
                            "Topi","Kongoni","Wildebeest","Giraffe","Eland", "Buffalo", "Warthog",
                            "Zebra","Elephant"))+ 
  
  scale_color_discrete(name="Guild",breaks=c("B","G","M"),
                       labels=c("Browser", "Grazer", "Mixed-Feeder"))+
  xlab("")+
  ggtitle("Nitrogen Distance Quadratic effect")+
  ylab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(text=element_text(size=14,  family="Arial"))+
  ylim(-2,2)


# Quadratic Effect of Distance to Permanent Water By Species


waterplot<-par_df_t %>%
  mutate(species.name=fct_relevel(species.name,"dikdik","thomson","grant","impala",
                                  "topi","kongoni","wildebeest","giraffe","eland","buffalo",
                                  "warthog","zebra","elephant"))%>%
  ggplot(aes(x=species.name,y=water_quad,color=guild))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(CI_low[12,][2:14]),ymax=as.numeric(CI_hi[12,][2:14]),width=0.4, size=0.8)+
  theme_bw()+
  scale_x_discrete(labels=c("Dikdik","Thomson's Gazelle","Grant's Gazelle", "Impala",
                            "Topi","Kongoni","Wildebeest","Giraffe","Eland", "Buffalo", "Warthog",
                            "Zebra","Elephant"))+ 
  scale_color_discrete(name="Guild",breaks=c("B","G","M"),
                       labels=c("Browser", "Grazer", "Mixed-Feeder"))+
  xlab("")+
  ggtitle("Water quadratic effect")+
  ylab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(text=element_text(size=14,  family="Arial"))+
  ylim(-2,5)

# Quadratic Effect of Runoff By Species

runoff<-par_df_t %>%
  mutate(species.name=fct_relevel(species.name,"dikdik","thomson","grant","impala",
                                  "topi","kongoni","wildebeest","giraffe","eland","buffalo",
                                  "warthog","zebra","elephant"))%>%
  ggplot(aes(x=species.name,y=runoff_quad,color=guild))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(CI_low[13,][2:14]),ymax=as.numeric(CI_hi[13,][2:14]),width=0.4, size=0.8)+
  theme_bw()+
  scale_x_discrete(labels=c("Dikdik","Thomson's Gazelle","Grant's Gazelle", "Impala",
                            "Topi","Kongoni","Wildebeest","Giraffe","Eland", "Buffalo", "Warthog",
                            "Zebra","Elephant"))+ 
  
  scale_color_discrete(name="Guild",breaks=c("B","G","M"),
                       labels=c("Browser", "Grazer", "Mixed-Feeder"))+
  xlab("")+
  ggtitle("Runoff quadratic effect")+
  ylab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(text=element_text(size=14,  family="Arial"))+
  ylim(-2,5)

# Quadratic Effect of Distance to Sodium Hotspots By Species


sodium<-par_df_t %>%
  mutate(species.name=fct_relevel(species.name,"dikdik","thomson","grant","impala",
                                  "topi","kongoni","wildebeest","giraffe","eland","buffalo",
                                  "warthog","zebra","elephant"))%>%
  ggplot(aes(x=species.name,y=sodium_quad,color=guild))+
  geom_point(size=2)+
  geom_hline(yintercept = 0, colour="grey",linetype=2, lwd=1)+
  geom_errorbar(ymin=as.numeric(CI_low[15,][2:14]),ymax=as.numeric(CI_hi[15,][2:14]),width=0.4, size=0.8)+
  theme_bw()+
  scale_x_discrete(labels=c("Dikdik","Thomson's Gazelle","Grant's Gazelle", "Impala",
                            "Topi","Kongoni","Wildebeest","Giraffe","Eland", "Buffalo", "Warthog",
                            "Zebra","Elephant"))+ 
  
  scale_color_discrete(name="Guild",breaks=c("B","G","M"),
                       labels=c("Browser", "Grazer", "Mixed-Feeder"))+
  xlab("")+
  ggtitle("Sodium quadratic effect")+
  ylab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(text=element_text(size=14,  family="Arial"))+
  ylim(-4,5)

# single species non linear responses

# Distance from permanent water sources


par_df

# intercept of thomson is -0.62

# beta for water distance is 0.56

# beta for water distance squared is -0.25

# vector of water distances from the normalised x matrix

water_distances_squared<-X[,12]
water_distances<-X[,11]

gazelle_abundance=-0.62+0.56*water_distances+(-0.25*water_distances_squared)

plot(gazelle_abundance~water_distances,main="Thomson's gazelle responses to distance from water")

# dik dik 

# intercept is -5.02
# beta for water distance -2.45
# beta for water distance squared 2.49

dikdik_abundance=-5.02+(-2.45*water_distances)+2.49*water_distances_squared
plot(dikdik_abundance~water_distances, main="Dikdik responses to distance from water")

# impala 

# intercept -3.34
# beta for water distance is -2.65
# beta for water distance squared is 1.25


impala_abundance=-3.34+(-2.65*water_distances)+1.25*water_distances_squared
plot(impala_abundance~water_distances, main="Impala responses to distance from water")

# grant's gazelle

# intercept -2
# beta for water distance 1.49
# beta for water distance squared -0.91

grants_abundance=-2+1.49*water_distances+(-0.91*water_distances_squared)

plot(grants_abundance~water_distances,main="Grant's gazelle response to distance from water")

# Giraffe 

# intercept -3.31
# beta for water distance  -1.82
# beta for water distance squared 1.65


giraffe_abundance=-3.31+(-1.82*water_distances)+1.65*water_distances_squared

plot(giraffe_abundance~water_distances,main="Giraffe's response to distance from water")

# nitrogen hotspots distance 

nitrogen_distances<-X[,6]

nitrogen_distances_squared<-X[,14]

# species that show non linear relationships are thomson's gazelle and Kongoni

# Thomson's gazelle

# intercept -0.62
# response to nitrogen 0.66
# response to nitrogen squared -0.28

thomson_abundance=-0.62+0.66*nitrogen_distances+(-0.28*nitrogen_distances_squared)

plot(thomson_abundance~nitrogen_distances, main="Thomson's gazelle abundance as a function of nitrogen distance")

# Kongoni 

# intercept -2.74
# response to distance from nitrogen 1.73
# response to distance from nitrogen squared -1.09

kongoni_abundance=-2.74+1.73*nitrogen_distances+(-1.09*nitrogen_distances_squared)
plot(kongoni_abundance~nitrogen_distances, main="Kongoni's abundance as a function of nitrogen distance")


# Runoff (not a distance)

runoff<-X[,7]
runoff_squared<-X[,13]

# kongoni 

# intercept -2.74
# response to runoff 1.3
# response to runoff squared -1.65  

kongoni_abundance_runoff=-2.74 + 1.3*runoff + (-1.65*runoff_squared)

plot(kongoni_abundance_runoff~runoff, main="Kongoni abundance as a function of runoff")

# grant's gazelle
# intercept -2
# response to runoff 1.36
# response to runoff squared -1.34  

grants_abundance_runoff=-2 +1.36*runoff + (-1.34*runoff_squared)

plot(grants_abundance_runoff~runoff,main="Grant's gazelle abundance as a function of runoff")


