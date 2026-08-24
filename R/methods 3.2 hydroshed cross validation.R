#-----------------------------------------------------------------------------------------------------------------------------
# first look at hydrosheds data
#-----------------------------------------------------------------------------------------------------------------------------
# shape files
files <- list.files('../hydrosheds/hybas_af_lev01-12_v1c',pattern='.shp', full.names=TRUE)
files <- files[!grepl('.xml',files)]

# hydrology
load('hyd.cal.dates.RData')


# unique sites only
sites <- unique(hyd[,c('Latitude','Longitude')])

library(sf)
xlim <- c(-20,45)
ylim <- c(10,40)

sf_use_s2(FALSE) # stops errors due to curvature shit
#-----------------------------------------------------------------------------------------------------------------------------
# find best catchments
#-----------------------------------------------------------------------------------------------------------------------------
levels <- 1:6
N <- length(levels)
match <- mismatch <- numeric(N)
	for(n in 1:N){

		sf <- read_sf(files[n])  
		# filter the catchments to study area
		sf <- st_crop(sf, xmin=min(xlim), ymin=min(ylim), xmax=max(xlim), ymax=max(ylim))


		# CV. randomly split data into two  sets
		i <- sample(c(TRUE, FALSE),size = nrow(sites), replace=TRUE)
		d1 <- sites[i,]
		d1.sf <- st_as_sf(d1, coords = c("Longitude", "Latitude"), crs = 4326)
		d2 <- sites[!i,]
		d2.sf <- st_as_sf(d2, coords = c("Longitude", "Latitude"), crs = 4326)

	
		# allocate data points to their catchments
		dd1 <- as.data.frame(st_join(d1.sf, sf))[,'HYBAS_ID']
		dd2 <- as.data.frame(st_join(d2.sf, sf))[,'HYBAS_ID']
	
		# count how many catchment-timebins are active in both d1 and d2
		dd1 <- unique(dd1[!is.na(dd1)])
		dd2 <- unique(dd2[!is.na(dd2)])
		match[n] <- sum(dd1%in%dd2)
	
		# count how many are in one but not the other
		mismatch[n] <- sum(!dd1%in%dd2) + sum(!dd2%in%dd1)
		}

data.frame(level=levels, matches=match, mismatches=mismatch, profit=match-mismatch)


plot(sf$geometry)
points(d1$Longitude, d1$Latitude, pch=16, col=1)
points(d2$Longitude, d2$Latitude, pch=16, col=2)

library(maps)
plot(NULL, xlim=xlim, ylim=ylim, frame.plot=FALSE, axes=FALSE, xlab='', ylab='')
map('world', xlim=xlim, ylim=ylim, col='grey90', add=TRUE, fill=TRUE, border='grey')
points(d1$Longitude, d1$Latitude, pch=16, col=1)
points(d2$Longitude, d2$Latitude, pch=16, col=2)
#-----------------------------------------------------------------------------------------------------------------------------
# find best time bins
#-----------------------------------------------------------------------------------------------------------------------------
bw <- seq(100,2000,by=100)
data <- hyd
N <- length(bw)
for(w in bw){
	breaks <- unique(c(seq(2000,14000,by=w),14000)) 
	name <- paste0('w',w)
	data[,name] <- cut(data$wm.date, breaks, 1:(length(breaks)-1))
	}

# look at level 4
level <- 4
R <- 50 
all.matches <- all.mismatches <- matrix(,N,R)
for(r in 1:R){
	
	sf <- read_sf(files[level])  

	# filter the catchments to study area
	sf <- st_crop(sf, xmin=min(xlim), ymin=min(ylim), xmax=max(xlim), ymax=max(ylim))

	# CV. randomly split data into two  sets
	i <- sample(c(TRUE, FALSE),size = nrow(hyd), replace=TRUE)
	d1 <- data[i,]
	d1.sf <- st_as_sf(d1, coords = c("Longitude", "Latitude"), crs = 4326)
	d2 <- data[!i,]
	d2.sf <- st_as_sf(d2, coords = c("Longitude", "Latitude"), crs = 4326)

	names <- paste0('w',bw)
	summary <- as.data.frame(matrix(,N,2)); names(summary) <- c('binwidth','matches')
	for(n in 1:N){
		# allocate data points to their catchments
		d1 <- as.data.frame(st_join(d1.sf, sf))[,c(names[n],'HYBAS_ID')]
		d2 <- as.data.frame(st_join(d2.sf, sf))[,c(names[n],'HYBAS_ID')]
	
		# count how many catchment-timebins are active in both d1 and d2
		d1 <- unique(d1[!is.na(d1$HYBAS_ID),])
		d2 <- unique(d2[!is.na(d2$HYBAS_ID),])
		d1 <- paste(d1[,1],d1[,2])
		d2 <- paste(d2[,1],d2[,2])

		matches <- sum(d1%in%d2)
		mismatches <- sum(!d1%in%d2) + sum(!d2%in%d1)
		
		all.matches[n,r] <- matches
		all.mismatches[n,r] <- mismatches
		}
	}

plot(NULL, xlim=range(bw), ylim=range(all.matches), xlab='time bins (yrs)', ylab='CV matches in same catchment')
for(r in 1:R)points(bw,all.matches[,r])
#-----------------------------------------------------------------------------------------------------------------------------
# plot level 4 for some arbitary time bins
#-----------------------------------------------------------------------------------------------------------------------------
data <- hyd
binwidth <- 500
level <- 4
breaks <- unique(c(seq(2000,14000,by=binwidth),14000)) 
data$timebin <- cut(data$wm.date, breaks, 1:(length(breaks)-1))
N <- length(breaks)-1
bin.names <- c()
for(n in 1:N)bin.names[n] <- paste0(breaks[n]/1000, '-', breaks[n+1]/1000, 'kyrBP')
sf <- read_sf(files[level])  
	
# filter the catchments to study area
sf <- st_crop(sf, xmin=min(xlim), ymin=min(ylim), xmax=max(xlim), ymax=max(ylim))

pdf(file='../plots/catchments 0.5kyr.pdf', width=14, height=7)
par(mfrow=c(4,6), mar=c(1,1,1,1))
for(n in 1:N){

	# plot the overall catchments
	plot(sf$geometry, main=bin.names[n])

	# get data points in the first time bin
	sub <- subset(data, timebin==n)
	d.sf <- st_as_sf(sub, coords = c("Longitude", "Latitude"), crs = 4326)

	# allocate data points to their catchments
	d <- as.data.frame(st_join(d.sf, sf))[,c('timebin','HYBAS_ID')]

	# colour the catchments
	d.sub <- subset(sf, HYBAS_ID%in%unique(d$HYBAS_ID))
	plot(d.sub$geometry, add=TRUE, col='darkseagreen')

	# add the points
	points(sub$Longitude,sub$Latitude, pch=16, col=1)
	}
dev.off()
#-----------------------------------------------------------------------------------------------------------------------------