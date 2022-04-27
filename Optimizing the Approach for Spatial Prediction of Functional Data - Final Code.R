#65 Fourier basis functions used to smooth the average daily temperature raw data curve for Halifax
##Graphs present the raw data curve for Halifax (light blue curve) against a smoothed curve (dark blue). Smoothing is done via basis expansion with 5 (left) and 75 (right) Fourier basis functions.
smooth = c(5,65,75)
for(j in 1:3){
  nbasis = smooth[j]
  yearRng = c(0,365)
  daybasis = create.fourier.basis(yearRng, nbasis)
  tempav = CanadianWeather$dailyAv[, , 'Temperature.C']
  daytempfd <- with(CanadianWeather, smooth.basis(day.5,
                                                  tempav, daybasis,
                                                  fdnames=list("Day", "Station", "degrees C"))$fd )
  for(i in 2:2){
    plot(tempav[,i],axes=TRUE,xlab="Day",ylab="degrees C",
         main=CanadianWeather$place[i])
    lines(daytempfd[i],col=2)
    axisIntervals(1)
    axis(2)
    readline("Press <return to continue")
  }
}

#Penalized smoothing of the Canadian Weather daily average temperature data set with λ = 10^(−2)
daybasis365 = create.fourier.basis(c(0,365),365)

harmLfd = vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))

tempfdPar = fdPar(daybasis365,harmLfd,1e-2)
tempfd = smooth.basis(1:365,daily$tempav,tempfdPar)

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd,xlab='day',ylab='temperature',cex.lab=1.5,cex.axis=1.5,col = "plum3")

#Penalized smoothing of the Canadian Weather daily average temperature data set with λ = 10^(-7)

tempfdPar2 = fdPar(daybasis365,harmLfd,1e7)
tempfd2 = smooth.basis(1:365,daily$tempav,tempfdPar3)

quartz()
mar = c(8, 8, 4, 2)
plot(tempfd2$fd,xlab='Day',ylab='Temperature.C',cex.lab=1.5,cex.axis=1.5,col=4)


#The values of the generalized cross-validation (GCV) criterion for the log precipitation data (left) and daily average temperature data for Halifax (right). The roughness penalty was defined by harmonic acceleration.

logprecav <- CanadianWeather$dailyAv[dayOfYearShifted,,3] #Replace with 1 for daily average temperature
CanadianWeather
dayrange = c(0,365)
daybasis365 = create.fourier.basis(dayrange,365)
Lcoef = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, dayrange)
loglam = seq(2,9,0.25) #change boundaries accordingly
nlam = length(loglam)
dfsave = rep(NA,nlam)
gcvsave = rep(NA,nlam)
for (ilam in 1:nlam){
  cat(paste('log10 lambda =',loglam[ilam],'\n'))
  lambda = 10^loglam[ilam]
  fdParobj = fdPar(daybasis365, harmaccelLfd, lambda)
  smoothlist = smooth.basis(day.5, logprecav,fdParobj)
  dfsave[ilam] = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
} 
plot(loglam, gcvsave,xlab="log10lambda", ylab="GCV value", lwd=2, pch=1, type="o")
loglam[which.min(gcvsave)] 

#Time series signal displaying Gaussian noise and its complementary ACF.
t = 0:365
y_stationary <- rnorm(length(t),mean=1,sd=1)
y_stationary <- y_stationary/max(y_stationary) #normalize each for simplicity
plot(t,y_stationary)
abline(lm(y_trend~t))
frame()
par(mfcol=c(2,1))
#the stationary signal and ACF
plot(t,y_stationary,
     type='l',col='red',
     xlab="time (t)",
     ylab="Y(t)",
     main="Stationary signal")
acf(y_stationary,lag.max = length(y_stationary),
    xlab="lag #", ylab='ACF', main=' ')

#Time series signal displaying the cumulative sum of Gaussian noise and its complementary ACF
y_trend <- cumsum(rnorm(length(t),mean=1,sd=4))+t/100
y_stationary <- y_stationary/max(y_stationary)
y_trend <- y_trend/max(y_trend)
#the trend signal and ACF
plot(t,y_trend,
     type='l',col='red',
     xlab="time (t)",
     ylab="Y(t)",
     main = "Trend signal")
acf(y_trend,lag.max = length(y_trend),
    xlab = "lag #", ylab = 'ACF', main=' ')

#Left Halifax daily average temperature raw data (black crosses) and corresponding smoothed curve (purple line). (Right) The complementary ACF.
#Left) Halifax daily average precipitation raw data (black crosses) and corresponding smoothed curve (purple curve). (Right) The complementary ACF
nbasis = 365
yearRng = c(0,365)
daybasis = create.fourier.basis(yearRng, nbasis)
tempav = CanadianWeather$dailyAv[, , 1] #change to 2 for precipitation
daytempfd <- with(CanadianWeather, smooth.basis(day.5,tempav, daybasis, fdnames=list("Day", "Station", "degrees C"))$fd )
#daytempfd[1]
plot(tempav[,2],axes=TRUE,xlab="Day",ylab="Precipitation (mm)", main=CanadianWeather$place[2], pch =3)
lines(daytempfd[2],col='violet', lwd=4)
axisIntervals(1)
axis(2)
acf(tempav[,2],lag.max = length(tempav[,1]), xlab="lag #", ylab='ACF', main=' ',)

#Visualization of the Meuse river floodplain  and the relative lead concentrations
#Meuse dataset comprises of 4 heavy metals measured in the top soil in a flood plain along the river meuse. Apparently polluted sediment is carried by the river and mostly deposited close to the river bank
library(sp)
library(gstat)
data(meuse)
class(meuse)
names(meuse)
coordinates(meuse) = ~x+y
class(meuse)
summary(meuse)
coordinates(meuse) #164 spatial coordinates
spplot(meuse, "lead", do.log = T, colorkey = TRUE, main = "Lead concentrations (ppm)") #the x and y axis are the spatial coordinates
bubble(meuse, "zinc", col = c("#00ff0088"), main = "zinc concentrations (ppm)")
data(meuse.grid)
summary(meuse.grid)
str(meuse.grid) #convert from data frame to SpatialPixelsDataFrame
coordinates(meuse.grid) = ~x+y
gridded(meuse.grid) = TRUE
meuse$loglead <- log(meuse$lead)
image(meuse.grid["dist"])
title("distance to river (red=0)")

#Point pairs that have specific separation distances are plotted
library(gstat)
hscat(log(lead) ~ 1, meuse, (0:9) * 100)

#Variogram cloud for log(lead) values
#Variogram - we assume there is a constant trend for the variable log(zinc)
vgm1 = variogram(log(lead)~1, meuse)
#plot variogram cloud
meuse.vc <- variogram(log(lead) ~ 1, meuse, cloud = TRUE)
plot(meuse.vc, ylab=bquote(gamma), xlab=c("h (separation distance in m)"))

#Sample variogram fitted with theoretical variogram model for log(lead) values
vgm1
plot(vgm1)
#When comparing two variogram models, spherical vs Gaussian, using universal kriging, it seems using the spherical model is slightly better than the gaussian model
#Mean error should be closer to 0, RMSE should be smaller, MSDR should be closer to 1
vgm1.fit = fit.variogram(vgm1, model = vgm(0.05156318, "Sph", 965.1573, 0.51530771))
vgm1.fit
#plot the fitted variogram and the observed variogram on the same graph
plot(vgm1, vgm1.fit, ylab=bquote(gamma), xlab=c("h (separation distance in m)"))

#Representation of kriged lead concentrations for river Meuse floodplain, generated by Ordinary Kriging.
vgm1.kriged = krige(log(lead)~1, meuse, meuse.grid, model = vgm1.fit)
spplot(vgm1.kriged["var1.pred"])
vgm1.kriged

#Locations of 35 Weather Stations in Canadian Weather data set (left). Locations of 35 weather stations in Canadian Maritime Provinces data set (right).
install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
ggplot(data = world) +
  geom_sf()
ggplot(data = world) +
  geom_sf(color = "black", fill = "lightgrey") +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-150, -50), ylim = c(40, 80), expand = FALSE)+
  geom_point(data=data.frame(CanadianWeather$coordinates), aes(x=-(CanadianWeather$coordinates[,2]), y=CanadianWeather$coordinates[,1]), colour="red", fill="pink",pch=21, size=5, alpha=I(0.7))
ggplot(data = world) +
  geom_sf(color = "black", fill = "lightgrey") +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-120, -50), ylim = c(40, 70), expand = FALSE)+
  geom_point(data=data.frame(maritimes.coords), aes(x=maritimes.coords[,1], y=maritimes.coords[,2]), colour="red", fill="pink",pch=21, size=5, alpha=I(0.7))


#Both panels display the locations at which temperature measurements are taken (black dots). (Left) 10 randomly selected unobserved locations (red dots). (Right) A single fixed unobserved site located in Moncton (red dot)
library(geofd)
data(maritimes.data)
data(maritimes.coords)
#coords for Moncton 46.8199° N, 68.4766
n<-dim(maritimes.data)[1]
data(maritimes.newcoords)
argvals<-seq(1,n, by=1)
data(maritimes.avg)
col1<-sample( (min(maritimes.coords[,1])*100):(max(maritimes.coords[,1])*100),
              10, replace=TRUE)/100
col2<-sample( (min(maritimes.coords[,2])*100):(max(maritimes.coords[,2])*100),
              10, replace=TRUE)/100
new.coords <- cbind(col1,col2)
# Replace col1 with -64.6839 and col2 with 46.1053 for fixed prediction site at Moncton
zoom_to <- c(-64.6839, 46.1053)  #METAR/Synop Information for CYQM (71705) in Moncton, N. B., Canada
zoom_level <- 5
lon_span <- 360 / 2^zoom_level
lat_span <- 180 / 2^zoom_level
lon_bounds <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
lat_bounds <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)
ggplot() + geom_sf(data = worldmap) +
  geom_sf(data = st_sfc(st_point(zoom_to), crs = 4326),
          color = 'red', size = 3) +
  coord_sf(xlim = lon_bounds, ylim = lat_bounds) +
  geom_point(data=data.frame(maritimes.coords), aes(x=maritimes.coords[,1], y=maritimes.coords[,2]), color="black", fill="black",pch=19, size=2, alpha=I(0.7))+
  #geom_point(data=data.frame(new.coords), aes(x=new.coords[,1], y=new.coords[,2]), color="red",pch=19, size=3)+
  labs(x = "Longitude") +
  labs(y = "Latitude") +
  theme_bw()

#The smoothed temperature curves for the 35 visited locations (left) and the prediction curve for the unvisited site at Moncton (right)
# Prediction by okfd
okfd.res<-okfd(new.coords=new.coords, coords=maritimes.coords,
               data=maritimes.data, smooth.type="fourier",
               nbasis=65, argvals=argvals)
okfd.res$datafd
par(mfrow=c(1,1))
new.coords
plot(okfd.res$datafd, lty=1, col=8,
     main="Smoothed", xlab="Day", ylab="Temperature (Degrees C)")

plot(okfd.res$argvals, okfd.res$krig.new.data, col=1, lwd=2,type="l", lty=1, main="Predictions", xlab="Day",ylab="Temperature (Degrees C)")
lines(maritimes.avg,  type="p", pch=20,cex=0.5, col=2, lwd=1)

#OKFD predictions at ten randomly selected sites from Canadian Maritimes Provinces
matplot(okfd.res$argvals, okfd.res$krig.new.data, col=1, lwd=1, type="l", lty=1,
        main="Predictions", xlab="Day", ylab="Temperature (Degrees C)")
