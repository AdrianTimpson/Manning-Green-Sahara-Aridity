#--------------------------------------------------------------------
library(sf)
require(maps)
source('functions.R')
cols <- readLines('../tools/colours.txt')
d1 <- read.csv('../tools/animals.summary.csv')
load('R objects/hyd.cal.dates.RData')
data <- hyd

# assign data top time slices
breaks <- seq(2000,14000,by=500)
d1$timebin <- cut(d1$meanBP, breaks, 1:(length(breaks)-1))
N <- length(breaks)-1
slice.names <- c()
for(n in 1:N)slice.names[n] <- paste0(breaks[n]/1000, '-', breaks[n+1]/1000, 'kyrBP')

# plot region
plot.xlim <- c(-20,45)
plot.ylim <- c(13,40)

par(mfrow=c(4,6),mar=c(2,2,3,1))
for(n in 1:N){
	plot(NULL, xlim=plot.xlim, ylim=plot.ylim, bty='n', xaxs='i',yaxs='i',axes=TRUE, xlab='', ylab='',main=slice.names[n])
	polygon(x=c(plot.xlim,rev(plot.xlim)), y=rep(plot.ylim,each=2),col='steelblue',border='steelblue')
	map('world', xlim=plot.xlim, ylim=plot.ylim, col='grey90', add=TRUE, fill=TRUE, border='grey90')
	sub <- subset(d1, timebin==n)
	points(sub$Longitude, sub$Latitude,pch=16, col=sub$WDI.mean.col)
	}
#--------------------------------------------------------------------




