#-----------------------------------------------------------------------------------------------------------------------------
# Cross validation method to choosing catchments
#-----------------------------------------------------------------------------------------------------------------------------
# shape files
files <- list.files('../../hydrosheds/hybas_af_lev01-12_v1c',pattern='.shp', full.names=TRUE)
files <- files[!grepl('.xml',files)]

# hydrology
load('R objects/hyd.cal.dates.RData')

# unique sites only
sites <- unique(hyd[,c('Latitude','Longitude')])
xlim <- c(-19,35)
ylim <- c(13,37)

library(sf)
sf_use_s2(FALSE) # stops errors due to curvature shit
#-----------------------------------------------------------------------------------------------------------------------------
# find best catchments
#-----------------------------------------------------------------------------------------------------------------------------
levels <- 1:6
R <- 1000
N <- length(levels)
match.025 <- match.500 <- match.975 <- mismatch.025 <- mismatch.500 <- mismatch.975 <- profit.025 <- profit.500 <- profit.975 <- numeric(N)
for(n in 1:N){
	match <- mismatch <- profit <- numeric(R)
	sf <- read_sf(files[n])  
	# filter the catchments to study area
	sf <- st_crop(sf, xmin=min(xlim), ymin=min(ylim), xmax=max(xlim), ymax=max(ylim))

	for(r in 1:R){
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
		match[r] <- sum(dd1%in%dd2)
	
		# count how many are in one but not the other
		mismatch[r] <- sum(!dd1%in%dd2) + sum(!dd2%in%dd1)

		# count profit
		profit[r] <- match[r] - mismatch[r]
		}
	
	# store quantiles
	match.025[n] <- quantile(match, prob=0.025)
	match.500[n] <- quantile(match, prob=0.5)
	match.975[n] <- quantile(match, prob=0.975)
	mismatch.025[n] <- quantile(mismatch, prob=0.025)
	mismatch.500[n] <- quantile(mismatch, prob=0.5)
	mismatch.975[n] <- quantile(mismatch, prob=0.975)
	profit.025[n] <- quantile(profit, prob=0.025)
	profit.500[n] <- quantile(profit, prob=0.5)
	profit.975[n] <- quantile(profit, prob=0.975)
	}

summary <- data.frame(levels, match.025, match.500, match.975, mismatch.025, mismatch.500, mismatch.975, profit.025, profit.500, profit.975)
save(summary, file='R objects/catchment.summary.RData')
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
R <- 1000 
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

# generate a summary of quantiles
q.025 <- q.500 <- q.975 <- numeric(N)
for(n in 1:N){
	q.025[n] <- quantile(all.matches[n,], prob=0.025)
	q.500[n] <- quantile(all.matches[n,], prob=0.500)
	q.975[n] <- quantile(all.matches[n,], prob=0.975)
	}

summary <- data.frame(bw, q.025, q.500, q.975)
save(summary, file='R objects/timebin.summary.RData')

plot(NULL, xlim=range(bw), ylim=range(all.matches), xlab='time bins (yrs)', ylab='CV matches in same catchment',xaxt='n')
axis(1,bw,bw,las=2)
lines(bw,q.025,col='grey')
lines(bw,q.500,lwd=2)
lines(bw,q.975,col='grey')
points(bw,q.500,pch=16)
points(bw,q.025,pch=16,col='grey')
points(bw,q.975,pch=16,col='grey')
legend(x=1800,y=140, pch=c(16,16),col=c('black','grey'),legend=c('median','95% CI'),bty='n')
#-----------------------------------------------------------------------------------------------------------------------------

