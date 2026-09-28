
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
get.faunal <- function(){

	# This function pulls faunal data direct from database, 
	# cleaning etc, various merges etc
	require(BIADconnect)

	# get data from database
	d.eco <- query.database("SELECT TaxonCode,WaterCode FROM greensahara.FaunalTaxaList")
	d.spe <- query.database("SELECT TaxonCode,NISP,PhaseID,PresenceOnly FROM greensahara.FaunalSpecies")
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

#--------------------------------------------------------------------------------------------------

