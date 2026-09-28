#------------------------------------------------------------------------------------------------------------
# Fig 5 hydrology and WDI through continuous time
#------------------------------------------------------------------------------------------------------------
source('functions.R')

load('R objects/hyd.cal.dates.RData')
faunal <- read.csv('../tools/faunal phases.csv')
kd.pal <- stat.pal <- readLines('../tools/colours.txt')

minDate <- 2000
maxDate <- 14000

dim(hyd)
length(unique(hyd$SiteID))
#------------------------------------------------------------------------------------------------------------
# kernel density of hydrological dates
#------------------------------------------------------------------------------------------------------------
# kernel density
kd <- MASS::kde2d(x=hyd$wm.date, y=hyd$Latitude, h=c(1000,2),n =100, lims=c(minDate,maxDate,13,37)) # 17 to 32

# rescale
kd$z <- sqrt(kd$z)
#------------------------------------------------------------------------------------------------------------
# add kd value at each faunal phase, for later correlation
#------------------------------------------------------------------------------------------------------------
N <- nrow(faunal)
for(n in 1:N){
	x.dist <- abs(kd$x-faunal$meanBP[n])
	x.i <- which(x.dist==min(x.dist))
	y.dist <- abs(kd$y-faunal$Latitude[n])
	y.i <- which(y.dist==min(y.dist))
	faunal$kd.value[n] <- kd$z[x.i,y.i]
	}
#------------------------------------------------------------------------------------------------------------
# various objects
stat <- faunal$WDI.mean
cols <- faunal$WDI.mean.col
stat.inc <- round(seq(min(stat),max(stat),length.out=5),1)
stat.cols <- convert.stat.to.colours(stat.inc, col=stat.pal, v=stat.inc)
kd.inc <- seq(0,1,length.out=20)
kd.cols <- convert.stat.to.colours(kd.inc, col=kd.pal, v=kd.inc)
cor <- cor(faunal$WDI.mean, faunal$kd.value)
cor.test(faunal$WDI.mean, faunal$kd.value)
#------------------------------------------------------------------------------------------------------------
# plot
#------------------------------------------------------------------------------------------------------------
pdf(file='../plots/fig 4 kernel density.pdf',width=14,height=10)
par(mar=c(4,4,1,4))

image(kd, col=colorRampPalette(rev(kd.pal))(30), xlab='yrsBP', ylab='Latitude')
points(hyd$wm.date, hyd$Latitude, pch=3, col='grey')
points(faunal$meanBP, faunal$Latitude, col='black', pch=16, cex=1.3)
points(faunal$meanBP, faunal$Latitude, col=cols, pch=16)

# legend is very fiddly
X <- c(10500,12700,13000)
Y <- c(31,29,28,26)

text(x=X[2],y=Y[1], 'Animal WDI (mean)',cex=0.9, pos=2)
for(i in 1:length(stat.cols)){
	y <- Y[1] - (i-1)*0.3
	points(X[2],y,pch=16,col='black', cex=1.3)
	points(X[2],y,pch=16,col=stat.cols[i])
	text(X[3], y, labels=stat.inc[i], cex=0.7, pos=4, offset=0) 
	}

text(x=X[2],y=Y[2], 'Hydrological 14C dates',cex=0.9, pos=2)
points(x=X[2], y=Y[2], pch=3, col='grey',cex=1.5)

text(x=X[2],y=Y[3], 'Hydrological 14C density',cex=0.9, pos=2)
y.scale <- 15
for (i in 1:(length(kd.inc)-1)) {
	y <- Y[3] - (i-1)/y.scale
	rect(X[2],y,X[3],y-1/y.scale, col=kd.cols[i], border=NA)
	}
rect(X[2],y-1/y.scale,X[3],Y[3],col=NA, border='grey')

text(x=X[3], y=Y[3], 'high density', cex=0.7, pos=4, offset=0)
text(x=X[3], y=y-1/y.scale, 'low density', cex=0.7, pos=4, offset=0)
text(x=X[3],y=Y[4], paste("Pearson's R =",round(cor,2)),cex=0.7, pos=2)

dev.off()

#------------------------------------------------------------------------------------------------------------
