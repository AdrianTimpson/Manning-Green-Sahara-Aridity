#---------------------------------------------------------------
# aggregate worldclim over several years to provide average rainfall
#---------------------------------------------------------------
folder <- 'wc2.1_cruts4.09_10m_prec_1950-1959'
files <- list.files(folder, full.names=TRUE, pattern='.tif')

library(terra)
area <- ext(-20, 38, -20, 38)

# mean precipitation for each month
months <- sprintf("%02d", 1:12)
seasons <- c('Dec, Jan, Feb','Mar, Apr, May','Jun, Jul, Aug','Sep, Oct, Nov')
months <- c('12.tif|01.tif|02.tif','03.tif|04.tif|05.tif','06.tif|07.tif|08.tif','09.tif|10.tif|11.tif')
cols <- rev(terrain.colors(5)[1:4])
pdf('worldclim.pdf', width=10,height=8)
par(mfrow=c(2,2),mar=c(1,1,1,1))
for(s in 1:length(months)){
	yrs <- files[grepl(months[s], files)]
	N <- length(yrs)
	d <- rast(yrs[1])
	d <- crop(d, area)
	for(n in 2:N){
		d1 <- rast(yrs[n])
		d1 <- crop(d1,area)
		d <- d + d1
		}
	d <- d/N

	# thresholds 
	t1 <- 100/12
	t2 <- 200/12
	t3 <- 500/12
	d[d<=t1] <- 1	
	d[d>t1 & d<=t2] <- 2
	d[d>t2 & d<=t3] <- 3
	d[d>t3] <- 4
	plot(d, main=seasons[s], col=cols,colNA='skyblue', legend=FALSE)
	legend(x=-19,y=10,fill=cols, legend=c('<100','<200','<500','>500'),bty='n',cex=0.7, title='Rainfall rate (mm/yr)')
	}

#---------------------------------------------------------------
# look at rainfall by latitude within study area
#-------------------------------------------------------------------------
par(mfrow=c(1,1),mar=c(4,4,1,1))
P <- 100
cols <- colorRampPalette(c('red','blue','green','blue','red'))(12) 
months <- sprintf("%02d", 1:12)
lats <- seq(-20,38,length.out=P)
plot(NULL, xlim=c(-20,38), ylim=c(0,350), xlab='latitude', ylab='rainfall (mm)')
for(s in 1:length(months)){
	yrs <- files[grepl(months[s], files)]
	N <- length(yrs)
	d <- rast(yrs[1])
	d <- crop(d, area)
	for(n in 2:N){
		d1 <- rast(yrs[n])
		d1 <- crop(d1,area)
		d <- d + d1
		}
	d <- d/N
	coord <- crds(d, na.rm=FALSE)
	res <- matrix(,P-1,2)
	for(p in 2:P){
		min <- lats[p-1]
		max <- lats[p]
		mid <- (min+max)/2
		i <- coord[,'y']>=min & coord[,'y']<max
		x <- d[i]
		m <- mean(x[,1],na.rm=TRUE)
		res[p-1,] <- c(mid,m)
		}
	res <- as.data.frame(res); names(res) <- c('latitude','rainfall')
	lines(res, lwd=2, col=cols[s])
	}
legend(bty='n', x=28, y=250, legend=months, lwd=2, col=cols, cex=0.8)
#dev.off()
#-------------------------------------------------------------------------


