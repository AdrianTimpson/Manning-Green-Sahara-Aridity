#------------------------------------------------------------------------------------------------------------
# overheads
#------------------------------------------------------------------------------------------------------------
source('functions.R')
source('get hydrology.R')
dim(d)# 1323
#------------------------------------------------------------------------------------------------------------
# generate bands
#------------------------------------------------------------------------------------------------------------
N <- 15
minDate <- 2000
maxDate <- 15000
lat.band.key <- create.band.key(d=d, N=N, direction='Latitude')
long.band.key <- create.band.key(d=d, N=N, direction='Longitude')
#--------------------------------------------------------------------------------------------------
# generate SPDs for bands
#--------------------------------------------------------------------------------------------------
SPD.lat <- make.band.spds(d=d, minDate=minDate, maxDate=maxDate, inc=20, band.key=lat.band.key, direction='Latitude')
SPD.long <- make.band.spds(d=d, minDate=minDate, maxDate=maxDate, inc=20, band.key=long.band.key, direction='Longitude')
#--------------------------------------------------------------------------------------------------
# calculate 75% quantile on the SPDs
#--------------------------------------------------------------------------------------------------
lat.band.key <- cbind(lat.band.key, SPD.quantiles(SPD.lat,prob=c(0.125,0.875)))
long.band.key <- cbind(long.band.key, SPD.quantiles(SPD.long,prob=c(0.125,0.875)))
range(lat.band.key$count)
#--------------------------------------------------------------------------------------------------
# Fit models to onset/decay
#--------------------------------------------------------------------------------------------------
lat.mod.onset <- lm(onset~poly(mean,2,raw=TRUE), data=lat.band.key) 
lat.mod.decay <- lm(decay~poly(mean,2,raw=TRUE), data=lat.band.key) 
long.mod.onset <- lm(onset~poly(mean,3,raw=TRUE), data=long.band.key) # 3rd order
long.mod.decay <- lm(decay~poly(mean,1,raw=TRUE), data=long.band.key) # linear
summary(lat.mod.onset)
summary(lat.mod.decay)
summary(long.mod.onset)
summary(long.mod.decay)

lat <- seq(min(lat.band.key$mean),max(lat.band.key$mean),length.out=50)
lat.onset <- predict(lat.mod.onset, data.frame(mean=lat))
lat.decay <- predict(lat.mod.decay, data.frame(mean=lat))
lat.model <- data.frame(mean=lat, onset=lat.onset, decay=lat.decay)

long <- seq(min(long.band.key$mean),max(long.band.key$mean),length.out=50)
long.onset <- predict(long.mod.onset, data.frame(mean=long))
long.decay <- predict(long.mod.decay, data.frame(mean=long))
long.model <- data.frame(mean=long, onset=long.onset, decay=long.decay)
#--------------------------------------------------------------------------------------------------
# Plots!!
#--------------------------------------------------------------------------------------------------
for(ext in c('png','pdf')){
	for(choice in c('latitude','longitude')){
if(ext=='pdf'){
	plot.fnc <- pdf
	cex <- 1
	width <- 13
	height <- 7
	}
if(ext=='png'){
	plot.fnc <- png
	cex <- 6
	width <- 5000
	height <- 2500
	}
if(choice=='latitude'){
	file <- paste('../plots/hydro ',choice,' bands.',ext,sep='')
	band.key <- lat.band.key
	SPDs <- SPD.lat
	model <- lat.model
	}
if(choice=='latitude'){
	band.key <- lat.band.key
	SPDs <- SPD.lat
	model <- lat.model
	yscale <- 300
	extra <- c(-0.3,1.5)
	}
if(choice=='longitude'){
	band.key <- long.band.key
	SPDs <- SPD.long
	model <- long.model
	yscale <- 600
	extra <- c(-0.6,5)
	}
file <- paste('../plots/hydro ',choice,' bands.',ext,sep='')
colfunc <- colorRampPalette(c("seagreen", "gold"))
cols <- colfunc(N)

legend.shift <- rep(4:0,3)
years <- as.numeric(colnames(SPDs))

plot.fnc(file=file,width=width,height=height)
par(mar=c(4,4,1,4), cex=cex)
plot(NULL, xlim=c(maxDate+1000,minDate),ylim=range(model$mean)+extra,xlab='yrs BP', ylab=choice,main="")

for(n in N:1){
	lines(years,SPDs[n,]*yscale+band.key$mean[n],col=cols[n],type='l',lwd=1*cex,xlim=c(maxDate,minDate),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
	i <- years<=band.key$onset[n] & years>=band.key$decay[n]
	poly.x <- c(years[i],rev(years[i]))
	poly.y <- c(SPDs[n,i],rep(0,sum(i)))
	polygon(x=poly.x,y=poly.y*yscale+band.key$mean[n],col=cols[n],border='grey',lwd=1)

	# add the legends
	K <- 65
	X1 <- c(minDate, minDate) - legend.shift[n]*K
	X2 <- c(minDate, minDate-30) - legend.shift[n]*K
	X3 <- minDate-40 - legend.shift[n]*K
	X4 <- maxDate+1500

	lines(x=X1,y=c(0,1)+band.key$mean[n],col=cols[n],lwd=1)
	lines(x=X2,y=c(0,0)+band.key$mean[n],col=cols[n],lwd=1)
	lines(x=X2,y=c(1,1)+band.key$mean[n],col=cols[n],lwd=1)
	text(x=X3, y=band.key$mean[n], labels='0',col='grey',cex=0.5,pos=4,offset=0)
	text(x=X3, y=1+band.key$mean[n], labels='1',col='grey',cex=0.5,pos=4,offset=0)

	text(x=X4, y=band.key$mean[n], paste(band.key$label[n],'deg'),cex=0.5,pos=4,offset=0)
	text(x=X4-1100, y=band.key$mean[n], paste('N=',band.key$count[n]),cex=0.5,pos=4,offset=0)
	}

lines(model$onset, model$mean,col='seagreen3', lwd=3*cex)
lines(model$decay, model$mean,col='orange', lwd=3*cex)

points(band.key$onset,band.key$mean,col='grey',pch=16,cex=1.2)
points(band.key$decay,band.key$mean,col='grey',pch=16,cex=1.2)
points(band.key$onset,band.key$mean,col='seagreen3',pch=16,cex=0.8)
points(band.key$decay,band.key$mean,col='orange',pch=16,cex=0.8)

legend(x=X4-1600, y=max(band.key$mean), cex=0.5,
	bg='white', box.col='white',
	pch=c(NA,15,16,16,NA),
	col=c('grey','grey','seagreen3','orange','grey'),
	lwd=c(1,NA,NA,NA,3)*cex,
	legend=c('Summed probability distribution (SPD)','75% quantile of SPD','Onset','Termination','Best fit curve'))
text(x=X4, y=max(band.key$mean)+2, paste('Total hydrologically active associated radiocarbon dates: N=',nrow(d)),cex=0.5,pos=4,offset=0)
dev.off()
}}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
