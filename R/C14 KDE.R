#-----------------------------------------------------------------------------------
# Eyeball the density of archaeological 14C through time slices
# Initially use euclidean space, but eventually use anisotropic kernels
# Initially use arbitaray bandwidth, but eventually use cross-validation
#-----------------------------------------------------------------------------------
# generate
load('hyd.cal.dates.RData')
dim(hyd)

plot(hyd$Longitude, hyd$Latitude)
rgl::plot3d(hyd$Longitude, hyd$Latitude, hyd$wm.date)

hist(hyd$wm.date, breaks=50)
range(hyd$wm.date)

breaks <- seq(2000,14000,by=1000)
N <- length(breaks)-1

plot(sub$Longitude, sub$Latitude, pch=16)

require(maps)
xlim <- c(-20,45)
ylim <- c(10,35)

n <- 12
sub <- subset(hyd, wm.date>=breaks[n] & wm.date<=breaks[n+1])
plot(NULL, xlim=xlim, ylim=ylim, frame.plot=FALSE, axes=FALSE, xlab='', ylab='', main=paste(breaks[n],'to',breaks[n+1],'BP'))
map('world', xlim=xlim, ylim=ylim, col='grey90', add=TRUE, fill=TRUE, border='grey')
points(sub$Longitude, sub$Latitude, pch=16)