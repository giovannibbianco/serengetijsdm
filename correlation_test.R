# correlation test time series matrices 

delta<-read.csv("Delta_NDVI.csv")
ndvi<-read.csv("ndvi_matrix.csv")
anomaly<-read.csv("anomaly_ndvi_matrix.csv")
grass<-read.csv("grass_matrix.csv")
rain<-read.csv("rain_matrix.csv")

cor(ndvi,grass)
cor(ndvi,delta)
cor(ndvi,anomaly)
cor(anomaly,delta)
cor(rain,grass)
cor(rain,ndvi)
cor(rain,anomaly)
