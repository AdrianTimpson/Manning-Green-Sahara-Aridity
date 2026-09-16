
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Functions and overheads
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
get.hydrology <- function(){

	# This function pulls hydrology data direct from database, 
	# cleaning etc, various merges etc
	require(BIADconnect)

	# get data from database
	conn  <-  init.conn(dbname ='greensahara')
	d.hyd <- query.database("SELECT * FROM greensahara.Hydrology")
	d.dat <- query.database("SELECT LabCode, SiteID, Method, C14Age, C14SD FROM greensahara.DatingSamples")
	d.sit <- query.database("SELECT SiteID, SiteName, Latitude, Longitude FROM greensahara.CommonSites")
	disconnect()

	# merge tables and keep only sites with hydrological data
	# only use C14 ... for now 'only' bother with the 96% of dates that are C14
	d <- merge(d.hyd,d.dat,by='LabCode')
	d <- merge(d,d.sit,by='SiteID')
	d <- subset(d, Method=='C14')
	i <- which(d$ProxyTypeLezine=='Arid' | d$ProxyTypeSP=='Arid' | d$ProxyTypeOther=='Arid')
	d <- d[-i,]
return(d)}
#--------------------------------------------------------------------------------------------------
get.animals <- function(){

	# This function hydrology pulls data direct from database, 
	# cleaning etc, various merges etc
	require(BIADconnect)

	# get data from database
	d.eco <- query.database("SELECT TaxonCode,WaterCode FROM greensahara.FaunalTaxaList")
	d.spe <- query.database("SELECT TaxonCode,NISP,PhaseID FROM greensahara.FaunalSpecies")
	d.pha <- query.database("SELECT SiteID,PhaseID,StartBP,EndBP FROM greensahara.CommonPhases")
	d.sit <- query.database("SELECT SiteID,Latitude,Longitude FROM greensahara.CommonSites")
	disconnect()

	d.ani <- merge(d.eco,d.spe,by='TaxonCode')
	d.comb <- merge(d.pha,d.ani, by='PhaseID')
	d.comb <- merge(d.comb,d.sit, by='SiteID')
	d.comb <- subset(d.comb, !is.na(WaterCode))

	d <- d.comb
return(d)}
#--------------------------------------------------------------------------------------------------
make.animal.summary <- function(d){

# MOST OF THIS SHOULD COME OUT!!
# NOTE ALSO THE PRESENCE HAS DIFFERENT LEVELS|!!
	# various summary statistics etc
	phases <- unique(d$PhaseID)
	N <- length(phases)
	names <- c('SiteID','PhaseID','StartBP','EndBP','Latitude','Longitude','WDI.mean','WDI.min','WDI.max','Shannon','Richness')
	summary <- as.data.frame(matrix(,N,length(names))) 
	names(summary) <- names

	for(n in 1:N){
		sub <- subset(d,PhaseID==phases[n])
		sub$NISP[is.na(sub$NISP)] <- 0
		summary$SiteID[n] <- unique(sub$SiteID)
		summary$PhaseID[n] <- unique(sub$PhaseID)
		summary$StartBP[n] <- unique(sub$StartBP)
		summary$EndBP[n] <- unique(sub$EndBP)
		summary$Latitude[n] <- unique(sub$Latitude)
		summary$Longitude[n] <- unique(sub$Longitude)
		summary$WDI.mean[n] <- mean(as.numeric(sub$WaterCode),na.rm=TRUE) 
		summary$WDI.min[n] <- min(as.numeric(sub$WaterCode),na.rm=TRUE) 
		summary$WDI.max[n] <- max(as.numeric(sub$WaterCode),na.rm=TRUE) 
		summary$Shannon[n] <- shannon(sub$NISP + 1) # add one to deal with presence only
		summary$Richness[n] <- nrow(sub)
		}
	summary$meanBP <- (summary$StartBP + summary$EndBP)/2
	summary <- subset(summary, !is.na(meanBP))

	# add colour for various stats

	cols <- readLines('../tools/colours.txt')

	summary$WDI.mean.col <- convert.stat.to.colours(summary$WDI.mean, cols, v=1:5)
	summary$Shannon.col <- convert.stat.to.colours(summary$Shannon, rev(cols), v=seq(min(summary$Shannon,na.rm=TRUE),max(summary$Shannon,na.rm=TRUE),length.out=5))
	summary$Richness.col <- convert.stat.to.colours(summary$Richness, rev(cols), v=seq(0,max(summary$Richness),length.out=5))

	# store summary and raw data	
	write.csv(summary, file='../tools/animals.summary.csv',row.names=FALSE)
	write.csv(d, file='../tools/animals.raw.csv',row.names=FALSE)

return(NULL)}

#--------------------------------------------------------------------------------------------------
create.band.key <- function(d, N, direction){

	# d: data including 'Latitude', 'Longitude', 'age', 'sd'
	# N: number of bands
	# direction either 'Latitude' or 'Longitude

	band.key <- as.data.frame(matrix(,N,6))
	names(band.key) <- c('band','count','min','max','mean','label')

	rank <- rank(d[,direction])
	posts <- seq(1,nrow(d),length.out=N+1)
	for(n in 1:N){
		sub <- d[rank>=posts[n] & rank<posts[n+1],]
		band.key$band[n] <- n
		band.key$count[n] <- nrow(sub)
		band.key$min[n] <- min(sub[,direction],na.rm=TRUE)
		band.key$max[n] <- max(sub[,direction],na.rm=TRUE)
		band.key$mean[n] <- mean(sub[,direction],na.rm=TRUE)
		band.key$label[n] <- paste(round(range(sub[,direction], na.rm=TRUE),1),collapse=' to ')
		}
return(band.key)}
#--------------------------------------------------------------------------------------------------
make.band.spds <- function(d, minDate, maxDate, inc, band.key, direction){

	require(ADMUR)
	d$datingType <- '14C'
	CalArray <- makeCalArray( calcurve=intcal20, calrange=c(minDate,maxDate), inc=inc )
	years <- CalArray$cal
	N <- max(band.key$band)
	SPDs <- matrix(,N,length(years))

	for(n in 1:N){
		sub <- d[d[,direction]>=band.key$min[n] & d[,direction]<=band.key$max[n],]
		sub$datingType <- '14C'
		spd <- summedCalibrator(data=sub,CalArray)[,1]
		SPDs[n,] <- spd/sum(spd)
		}
	colnames(SPDs) <- years
return(SPDs)}
#--------------------------------------------------------------------------------------------------
SPD.quantiles <- function(SPD,prob){
	N <- nrow(SPD)
	SPD.cumsum <- SPD
	for(n in 1:N)SPD.cumsum[n,] <- 1 - cumsum(SPD[n,]) 
	quant <- as.data.frame(matrix(,N,length(prob)))
	names(quant) <- c('onset','decay')
	years <- as.numeric(colnames(SPD))
	for(n in 1:N)quant[n,] <- approx(y=years, x = SPD.cumsum[n,], xout=prob)$y 
return(quant)}
#--------------------------------------------------------------------------------------------------
convert.stat.to.colours <- function(stat, cols, v){
	pal <- colorRampPalette(cols)
	i <- findInterval(stat,v)
	x <- pal(length(v))[i]
return(x)}
#--------------------------------------------------------------------------------------------------
data.in.kml <- function(data, kml.path, latcolname, longcolname, index=NULL){
	# returns the data in a particular polygon, use the full kml file path of the polygon
	# data: the dataframe containing lats and longs
	# latcolname, longcolname: the column names in the data object that specify latitude and longitude
	# index: can be used to select particular polygons only in the kml

	require(sf)
	pts <- st_as_sf(data, coords = c(longcolname,latcolname), crs = 4326)
	poly <- st_read(kml.path, crs = 4326)
	st_agr(pts) <- st_agr(poly) <- 'constant'
	pts.in.poly <- st_intersection(pts, poly)
	data.in.poly <- data[row.names(data)%in%row.names(pts.in.poly),]
return(data.in.poly)}
#--------------------------------------------------------------------------------------------------
getKMLnames <- function(file){	
	names <- readLines(file)
	names <- names[grep('\t\t\t<name>',names)]
	names <- gsub('\t\t\t<name>','',names)
	names <- gsub('</name>','',names)
return(names)}
#-----------------------------------------------------------------------------------------------
shannon <- function(x){
	p <- x/sum(x)
	res <- sum(p * log(p)) *(-1)
return(res)}
#-----------------------------------------------------------------------------------------------
KEplot <- function(kd, hyd, stat, date, lat, cols,  cor, file, X, Y, txt, stat.pal, kd.pal){

	pdf(file=file,width=13,height=7)
	par(mar=c(4,4,1,4))

	image(kd, col=colorRampPalette(rev(kd.pal))(30), xlab='yrsBP', ylab='Latitude')
	points(hyd$wm.date, hyd$Latitude, pch=3, col='grey')
	points(date, lat, col='black', pch=16, cex=1.3)
	points(date, lat, col=cols, pch=16)

	# legend is very fiddly
	stat.inc <- round(seq(min(stat),max(stat),length.out=5),1)
	stat.cols <- convert.stat.to.colours(stat.inc, col=stat.pal, v=stat.inc)

	kd.inc <- seq(0,1,length.out=20)
	kd.cols <- convert.stat.to.colours(kd.inc, col=kd.pal, v=kd.inc)

	text(x=X[2],y=Y[1], txt,cex=0.9, pos=2)
	for (i in 1:length(stat.cols)) {
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
	}
#--------------------------------------------------------------------------------------------------

