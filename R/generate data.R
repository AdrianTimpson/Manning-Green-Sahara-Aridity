#---------------------------------------------------------------------
# Not for public
# Script generates various downstream files direct from database
# Will not work without access to database!! 
#---------------------------------------------------------------------
source('functions.R')
require(BIADconnect)
#---------------------------------------------------------------------
# hydrological 14C dates
#---------------------------------------------------------------------
# get all from database
hyd <- get.hydrology() 

# remove data in Nile polygon and south of 13 degrees Latitude
nile <- data.in.kml(data=hyd,kml.path='../shapes/Nile polygon.kml',latcolname='Latitude', longcolname='Longitude',index=1)
hyd <- hyd[!row.names(hyd)%in%row.names(nile),]
hyd <- hyd[which(hyd$Latitude>=13),]

# create a kml of points
file <- "../shapes/hydrologically_active_sites.kml"
file.remove(file)
hyd.sites <- unique(hyd[,c('SiteID','Latitude','Longitude')])
hyd.sites <- sf::st_as_sf(hyd.sites, coords = c('Longitude','Latitude'), crs=4326)
sf::st_write(hyd.sites, file, layer='hydrology', driver="KML") 
tmp <- readLines(file)
tmp <- gsub('<name>.*</name>','<name></name>',tmp)
writeLines(tmp, con=file)

# Calibrate dates, store weighted mean. Takes 5 mins
library(ADMUR)
names(hyd)[names(hyd)=='C14Age'] <- 'age'
names(hyd)[names(hyd)=='C14SD'] <- 'sd'
minDate <- 2000
maxDate <- 14000
CalArray <- makeCalArray( calcurve=intcal20, calrange=c(minDate,maxDate) )
cal <- dateCalibrator(hyd, CalArray)

years <- as.numeric(row.names(cal))
hyd$wm.date <- round(colSums(cal*as.numeric(row.names(cal))))

# retain those within the date range of interest
i <- colSums(cal)>0.999999
hyd <- hyd[i,]
save(hyd, file='R objects/hyd.cal.dates.RData')
#---------------------------------------------------------------------
# faunal data, study area north of 13
#---------------------------------------------------------------------
# get all from database
faunal <- get.faunal()

# subset to study area and period
i <- which(faunal$StartBP<=14000 & faunal$EndBP>=2000 & faunal$Latitude>=13)
faunal <- faunal[i,]
nile <- data.in.kml(data=faunal,kml.path='../shapes/Nile polygon.kml',latcolname='Latitude', longcolname='Longitude',index=1)
faunal <- faunal[!row.names(faunal)%in%row.names(nile),]

# mid date
faunal$midBP <- (faunal$StartBP + faunal$EndBP)/2

# create a kml of points
file <- "../shapes/faunal_sites.kml"
file.remove(file)
faunal.sites <- unique(faunal[,c('SiteID','Latitude','Longitude')])
faunal.sites <- sf::st_as_sf(faunal.sites, coords = c('Longitude','Latitude'), crs=4326)
sf::st_write(faunal.sites, file, layer='faunal', driver="KML") 
tmp <- readLines(file)
tmp <- gsub('<name>.*</name>','<name></name>',tmp)
writeLines(tmp, con=file)

# create faunal summary data in the tools folder
write.csv(faunal, file='../tools/faunal.csv',row.names=FALSE)
#---------------------------------------------------------------------
# Archaeological sites (e.g.containing 14C but no faunal), study area north of 13
# Hydro has no phases, therefore all phases are archaeological
#---------------------------------------------------------------------
sit <- query.database("SELECT SiteID,Latitude,Longitude FROM greensahara.CommonSites WHERE Archaeology='1'")
pha <- query.database("SELECT SiteID,PhaseID,StartBP,EndBP FROM greensahara.CommonPhases")
pha <- merge(pha,sit,all.x=TRUE,all.y=FALSE)
disconnect()
# subset to study area and period
i <- which(pha$StartBP<=14000 & pha$EndBP>=2000 & pha$Latitude>=13)
pha <- pha[i,]
nile <- data.in.kml(data=pha,kml.path='../shapes/Nile polygon.kml',latcolname='Latitude', longcolname='Longitude',index=1)
arch <- pha[!row.names(pha)%in%row.names(nile),]
save(arch, file='R objects/arch.phases.RData')

# create a kml of points
file <- "../shapes/arch_sites.kml"
file.remove(file)
arch.sites <- unique(sit[,c('SiteID','Latitude','Longitude')])
arch.sites <- sf::st_as_sf(arch.sites, coords = c('Longitude','Latitude'), crs=4326)
sf::st_write(arch.sites, file, layer='archaeological', driver="KML") 
tmp <- readLines(file)
tmp <- gsub('<name>.*</name>','<name></name>',tmp)
writeLines(tmp, con=file)
#---------------------------------------------------------------------
# all sites
#---------------------------------------------------------------------
sit <- query.database("SELECT SiteID, Latitude, Longitude, Archaeology, Hydrology FROM greensahara.CommonSites")
# reduce to study area
nile <- data.in.kml(data=sit,kml.path='../shapes/Nile polygon.kml',latcolname='Latitude', longcolname='Longitude',index=1)
sit <- sit[!row.names(sit)%in%row.names(nile),]
sit <- sit[which(sit$Latitude>=13),]
save(sit, file='R objects/sites.RData')
#---------------------------------------------------------------------




