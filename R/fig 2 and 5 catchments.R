#-----------------------------------------------------------------------------------------------------------------------------
# plot level 4 catchments with hydrological and archaeological data
#-----------------------------------------------------------------------------------------------------------------------------
load('R objects/hyd.cal.dates.RData')
load('R objects/sites.RData')
load('R objects/arch.phases.RData')
faunal <- read.csv('../tools/faunal.csv')

library(sf)
library(maps)
library(jaccard)

# assign time bin to each hydrological date
binwidth <- 500
breaks <- unique(c(seq(2000,14000,by=binwidth),14000)) 
hyd$timebin <- cut(hyd$wm.date, breaks, 1:(length(breaks)-1))

# catchments as sf objects
catchments <- read_sf('../../hydrosheds/hybas_af_lev01-12_v1c/hybas_af_lev04_v1c.shp')  
catchments <- st_make_valid(catchments)
catchments <- st_transform(catchments, crs = 4326)
st_agr(catchments) <- 'constant'

# hydrology, faunal,arch and sites as sf objects
hyd <- st_as_sf(hyd, coords = c("Longitude", "Latitude"), crs = 4326)
st_agr(hyd) <- 'constant'

faunal <- st_as_sf(faunal, coords = c("Longitude", "Latitude"), crs = 4326)
st_agr(faunal) <- 'constant'

arch <- st_as_sf(arch, coords = c("Longitude", "Latitude"), crs = 4326)
st_agr(arch) <- 'constant'

sit <- st_as_sf(sit, coords = c("Longitude", "Latitude"), crs = 4326)
st_agr(arch) <- 'constant'

# reduce catchments to those that have some part above 13 degrees.
# # find catchments to grey out if partly below 13 and no hydro data above 13
N <- nrow(catchments)
decision <- c()
for(n in 1:N){

	# reset conditions to be safe
	cond1 <- cond2 <- col <- NA

	catchment <- catchments[n,]
	points <- st_intersection(hyd,catchment)

	limits <- st_bbox(catchment)
	if(limits$ymin >=13)decision[n] <- 'colour' # colour if entirely above 13
	if(limits$ymax <13)decision[n] <- 'remove' # remove if entirely below 13

	# grey out if partly below 13 and no hydro data above 13
	cond1 <- limits$ymin < 13 & limits$ymax>=13
	if(cond1){
		if(nrow(points)==0)decision[n] <- 'grey'
		if(nrow(points)>0){
			# grey out catchments partly south of study area with no data points
			cond2 <- limits$ymax < 13
			if(cond2)decision[n] <- 'grey'
			if(!cond2)decision[n] <- 'colour'
			}
		}	
	}
catchments$decision <- decision
catchments <- subset(catchments, decision!='remove')

# assign data types for each site
sit$Faunal <- 0
sit$Faunal[sit$SiteID%in%faunal$SiteID] <- 1
sit$class <- NA
sit$class[sit$Archaeology ==1 & sit$Hydrology==0 & sit$Faunal==0] <- 'archaeological'
sit$class[sit$Archaeology ==0 & sit$Hydrology==1 & sit$Faunal==0] <- 'hydrological'
sit$class[sit$Archaeology ==1 & sit$Hydrology==1 & sit$Faunal==0] <- 'mixed'
sit$class[sit$Archaeology ==1 & sit$Hydrology==0 & sit$Faunal==1] <- 'archaeological & faunal'
sit$class[sit$Archaeology ==0 & sit$Hydrology==1 & sit$Faunal==1] <- 'mixed'
sit$class[sit$Archaeology ==1 & sit$Hydrology==1 & sit$Faunal==1] <- 'mixed'

# assign data types for each site
arch$class <- 'archaeological'
arch$class[arch$PhaseID%in%faunal$PhaseID] <- 'archaeological & faunal'

# colours for legend
leg.cols <- c('firebrick','tomato','blue','purple')
classes <- names(table(sit$class))
class.names <- paste0(classes, ' n=', as.numeric(table(sit$class)))

# plot region
plot.xlim <- c(-20,45)
plot.ylim <- c(13,40)
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
# plot Fig 2A
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
# count number timeslices each catchment has points in.
d <- as.data.frame(st_join(hyd, catchments))[,c('timebin','HYBAS_ID')]
tb <- table(unique(d)$HYBAS_ID)
tb1 <- data.frame(ID=names(tb), count=as.numeric(tb))
tb2 <- data.frame(ID=catchments$HYBAS_ID[!catchments$HYBAS_ID%in%tb1$ID], count=0)
tb <- rbind(tb1, tb2)

# convert counts to colours
q <- unique(c(0,round(as.numeric(quantile(tb1$count,prob=seq(0,1,length.out=20))))))
nq <- length(q)
cols <- rev(readLines('../tools/colours.txt'))[1:3]
mypal <- colorRampPalette( cols )(nq)
counts <- min(q):max(q)
count.col <- c()
for(i in 1:nq)count.col[counts<=rev(q)[i]] <- rev(mypal)[i]

pdf(file='../plots/fig 2A catchments.pdf', width=14, height=7)

par(mfrow=c(1,1), mar=c(3,3,1,1))
plot(NULL, xlim=plot.xlim, ylim=plot.ylim, bty='n', xaxs='i',yaxs='i',axes=TRUE, xlab='', ylab='')
axis(1,at=min(plot.xlim):max(plot.xlim),labels=NA)
axis(2,at=min(plot.ylim):max(plot.ylim),labels=NA)
polygon(x=c(plot.xlim,rev(plot.xlim)), y=rep(plot.ylim,each=2),col='steelblue',border='steelblue')
map('world', xlim=plot.xlim, ylim=plot.ylim, col='grey90', add=TRUE, fill=TRUE, border='grey90')

# colour active catchments
N <- nrow(tb)
for(n in 1:N){
	id <- tb$ID[n]
	count <- tb$count[n]
	catchment <- subset(catchments, HYBAS_ID==id)
	if(catchment$decision=='colour')col <- count.col[counts==count]	
	if(catchment$decision=='grey')col <- 'grey90'	
	plot(catchment$geometry, add=TRUE, col=col, border='grey')
	}

# add hydrology, archaeology and faunal
cex <- 0.6
for(n in 1:length(classes)){
	sub <- subset(sit, class==classes[n])
	plot(sub$geometry, add=TRUE,  pch=16, col=leg.cols[n],cex=cex)
	}

# add legends
text(x=39, y=39, labels='Number of 0.5kyr timebins activated',cex=0.7)
phytools::add.color.bar(leg=10, lims=range(counts),cols=count.col, outline=FALSE, prompt=FALSE, digits=2, x=34,y=38.4,title='',subtitle='',lwd=10) 
legend(x=36, y=38, pch=16,col=leg.cols, legend=class.names, bty='n',cex=0.7)

dev.off()
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
# plot Fig 2B: six timeslices
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
arch$midBP <- (arch$StartBP + arch$EndBP)/2
slicewidth <- 2000
breaks <- unique(c(seq(2000,14000,by=slicewidth),14000)) 
N <- length(breaks)-1
slice.names <- c()
for(n in 1:N)slice.names[n] <- paste0(breaks[n]/1000, '-', breaks[n+1]/1000, 'kyrBP')

nq <- 4
count.col <- colorRampPalette( cols )(nq+1)
counts <- 0:nq

pdf(file='../plots/fig 2B catchments slices.pdf', width=14, height=6)
par(mfrow=c(2,3), mar=c(3,3,4,1))
N <- length(breaks)-1
for(n in 1:N){

	# base plot
	plot(NULL, xlim=plot.xlim, ylim=plot.ylim, bty='n', xaxs='i',yaxs='i',axes=TRUE, xlab='', ylab='',main=slice.names[n])
	polygon(x=c(plot.xlim,rev(plot.xlim)), y=rep(plot.ylim,each=2),col='steelblue',border='steelblue')
	map('world', xlim=plot.xlim, ylim=plot.ylim, col='grey90', add=TRUE, fill=TRUE, border='grey90')

	# plot all catchments in grey
	plot(catchments$geometry, add=TRUE, col='grey90', border='grey')

	# get data points in slice
	hyd.sub <- subset(hyd, wm.date>=breaks[n] & wm.date<=breaks[n+1])
	arch.sub <- subset(arch, midBP>=breaks[n] & midBP<=breaks[n+1])

	# count number timeslices each catchment has hydrological points in.
	d <- as.data.frame(st_join(hyd.sub, catchments))[,c('timebin','HYBAS_ID')]
	tb <- table(unique(d)$HYBAS_ID)
	tb1 <- data.frame(ID=names(tb), count=as.numeric(tb))
	tb2 <- data.frame(ID=catchments$HYBAS_ID[!catchments$HYBAS_ID%in%tb1$ID], count=0)
	tb <- rbind(tb1, tb2)

	# colour active catchments
	C <- nrow(tb)
	for(c in 1:C){
		id <- tb$ID[c]
		count <- tb$count[c]
		catchment <- subset(catchments, HYBAS_ID==id)
		if(catchment$decision=='colour')col <- count.col[counts==count]	
		if(catchment$decision=='grey')col <- 'grey90'	
		plot(catchment$geometry, add=TRUE, col=col, border='grey')
		}

	# add hydrology, arch and faunal
	cex <- 0.9
	points(hyd.sub, pch=16, col=leg.cols[3], cex=cex)
	points(subset(arch.sub, class=='archaeological') , pch=16, col=leg.cols[1], cex=cex)
	points(subset(arch.sub, class=='archaeological & faunal') , pch=16, col=leg.cols[2], cex=cex)
	}


# add legends to last plot
text(x=28, y=39, labels='Number of 0.5kyr timebins activated',cex=0.9)
phytools::add.color.bar(leg=30, lims=c(0,nq),cols=count.col, outline=FALSE, prompt=FALSE, digits=2, x=13,y=38,title='',subtitle='',lwd=8) 
legend(x=18, y=38, pch=16,col=leg.cols[1:3], legend=classes[1:3], bty='n',cex=0.9)

dev.off()
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
# plot Fig 5: 24 timeslices with faunal and non-hydro sites, and jaccard
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
slicewidth <- 500
breaks <- unique(c(seq(2000,14000,by=slicewidth),14000)) 
N <- length(breaks)-1
slice.names <- c()
for(n in 1:N)slice.names[n] <- paste0(breaks[n]/1000, '-', breaks[n+1]/1000, 'kyrBP')

# store stats for jaccard distance
jac <- p <- numeric(N-1)

pdf(file='../plots/fig 5 catchments slices.pdf', width=25, height=10)
par(mfrow=c(4,6), mar=c(2,2,2,1))

for(n in 1:(N-1)){

	# get data points in slice
	hyd.sub <- subset(hyd, wm.date>=breaks[n] & wm.date<=breaks[n+1])
	arch.sub <- subset(arch, midBP>=breaks[n] & midBP<=breaks[n+1])

	# identify which catchments have hydro data, arch data
	hyd.d <- unique(as.data.frame(st_join(hyd.sub, catchments))[,c('HYBAS_ID')])
	arch.d <- unique(as.data.frame(st_join(arch.sub, catchments))[,c('HYBAS_ID')])

	# calculate jaccard distance in each time slice
	hyd.vec <- catchments$HYBAS_ID%in%hyd.d
	arch.vec <- catchments$HYBAS_ID%in%arch.d
	jt <- jaccard.test(hyd.vec,arch.vec,method='exact')
	jac[n] <- jt$statistics
	p[n] <- jt$pvalue 

	# base plot
	plot(NULL, xlim=plot.xlim, ylim=plot.ylim, bty='n', xaxs='i',yaxs='i',axes=TRUE, xlab='', ylab='',main=paste0(slice.names[n],' p=',round(p[n],3)),cex.axis=0.7,cex.main=0.9,mgp=c(3,0.2,0),tck=-0.05)
	polygon(x=c(plot.xlim,rev(plot.xlim)), y=rep(plot.ylim,each=2),col='steelblue',border='steelblue')
	map('world', xlim=plot.xlim, ylim=plot.ylim, col='grey90', add=TRUE, fill=TRUE, border='grey90')

	# plot all catchments in grey
	plot(catchments$geometry, add=TRUE, col='grey90', border='grey')

	# colour active catchments
	C <- length(hyd.d)
	for(c in 1:C){
		if(length(hyd.d)>0){
			id <- hyd.d[c]
			catchment <- subset(catchments, HYBAS_ID==id)	
			plot(catchment$geometry, add=TRUE, col=cols[3], border='grey')
			}
		}
	# add hydrology, arch and faunal
	cex <- 0.9
	points(hyd.sub, pch=16, col=leg.cols[3], cex=cex)
	points(subset(arch.sub, class=='archaeological') , pch=16, col=leg.cols[1], cex=cex)
	points(subset(arch.sub, class=='archaeological & faunal') , pch=16, col=leg.cols[2], cex=cex)
	}

# add legends to last tile
plot(NULL, xlim=c(0,1),ylim=c(0,1),bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
legend(x=0, y=1, pch=16,col=leg.cols[1:3], legend=classes[1:3], bty='n',cex=0.9)
legend(x=0, y=0.5, pch=15,col=cols[3], legend='hydrologically active catchment', bty='n',cex=0.9)
dev.off()

# jaccard summary
jac
round(p,3)
sum(p<0.02)
#-----------------------------------------------------------------------------------------------------------------------------
