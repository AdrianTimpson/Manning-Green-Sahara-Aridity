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
#-----------------------------------------------------------------------------------------------------------------------------
# find best catchments
#-----------------------------------------------------------------------------------------------------------------------------
levels <- 1:6
R <- 1000
N <- length(levels)
match.025 <- match.500 <- match.975 <- mismatch.025 <- mismatch.500 <- mismatch.975 <- profit.025 <- profit.500 <- profit.975 <- n.catch <- numeric(N)
for(n in 1:N){
	print(paste('level',n))
	match <- mismatch <- profit <- numeric(R)
	sf <- read_sf(files[n]) 
	sf <- st_make_valid(sf)
	sf <- st_transform(sf, crs = 4326)
	st_agr(sf) <- 'constant'
 
	# filter the catchments to study area
	sf <- st_crop(sf, xmin=min(xlim), ymin=min(ylim), xmax=max(xlim), ymax=max(ylim))
	
	# allocate data points to their catchments
	d.sf <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326)
	d <- as.data.frame(st_join(d.sf, sf))
	d <- d[!is.na(d$HYBAS_ID),]

	# store total number of catchments
	n.catch[n] <- nrow(sf)

	for(r in 1:R){
		# CV. randomly split data into two  sets
		i <- sample(c(TRUE, FALSE),size = nrow(d), replace=TRUE)
		d1 <- unique(d[i,'HYBAS_ID'])
		d2 <- unique(d[!i,'HYBAS_ID'])
	
		# count how many catchment-timebins are active in both d1 and d2
		match[r] <- sum(d1%in%d2)
	
		# count how many are in one but not the other
		mismatch[r] <- sum(!d1%in%d2) + sum(!d2%in%d1)

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

summary <- data.frame(levels, n.catch, match.025, match.500, match.975, mismatch.025, mismatch.500, mismatch.975, profit.025, profit.500, profit.975)
save(summary, file='R objects/catchment.summary.RData')
#-----------------------------------------------------------------------------------------------------------------------------
# find best time bins
#-----------------------------------------------------------------------------------------------------------------------------
bw <- seq(100,2000,by=100)
names <- paste0('w',bw)
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

sf <- read_sf(files[level])  
sf <- st_make_valid(sf)
sf <- st_transform(sf, crs = 4326)
st_agr(sf) <- 'constant'

# filter the catchments to study area
sf <- st_crop(sf, xmin=min(xlim), ymin=min(ylim), xmax=max(xlim), ymax=max(ylim))

# allocate data points to their catchments
d.sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)
d <- as.data.frame(st_join(d.sf, sf))
d <- d[!is.na(d$HYBAS_ID),]

for(r in 1:R){

	if(r%in%seq(1,R,length.out=10))print(paste('r=',r))# housekeeping

	# CV. randomly split data into two  sets
	i <- sample(c(TRUE, FALSE),size = nrow(d), replace=TRUE)
	d1 <- d[i,]
	d2 <- d[!i,]

	# for each timebin (n) count how many catchments are active in both d1 and d2
	for(n in 1:N){

		dd1 <- unique(d1[,c(names[n],'HYBAS_ID')])
		dd2 <- unique(d2[,c(names[n],'HYBAS_ID')])

		dd1 <- paste(dd1[,1],dd1[,2])
		dd2 <- paste(dd2[,1],dd2[,2])

		matches <- sum(dd1%in%dd2)
		mismatches <- sum(!dd1%in%dd2) + sum(!dd2%in%dd1)
		
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

