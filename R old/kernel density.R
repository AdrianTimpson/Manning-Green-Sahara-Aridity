#------------------------------------------------------------------------------------------------------------
# various plots of ecological variables against KD of hydro dates
#------------------------------------------------------------------------------------------------------------
source('functions.R')

hyd <- get.hydrology()
ani <- get.animals()

minDate <- 2000
maxDate <- 15000

dim(hyd)
length(unique(hyd$SiteID))
#------------------------------------------------------------------------------------------------------------
# Calibration takes a while, do once and save
#------------------------------------------------------------------------------------------------------------
# library(ADMUR)

# CalArray <- makeCalArray( calcurve=intcal20, calrange=c(minDate,maxDate) )
# cal <- dateCalibrator(hyd, CalArray)

# retain those within the date range of interest
# i <- colSums(cal)>0.999999

# years <- as.numeric(row.names(cal))
# hyd$wm.date <- colSums(cal*as.numeric(row.names(cal)))
# hyd <- hyd[i,]
#save(hyd, file='hyd.cal.dates.RData')

load('hyd.cal.dates.RData')
dim(hyd)
#------------------------------------------------------------------------------------------------------------
# kernel density
kd <- MASS::kde2d(x=hyd$wm.date, y=hyd$Latitude, h=c(1000,2),n =100, lims=c(minDate,maxDate,17,32))

# rescale
kd$z <- sqrt(kd$z)

# animal summary
summary <- make.animal.summary(ani)

# add kd value to summary
N <- nrow(summary)
for(n in 1:N){
	x.dist <- abs(kd$x-summary$meanBP[n])
	x.i <- which(x.dist==min(x.dist))
	y.dist <- abs(kd$y-summary$Latitude[n])
	y.i <- which(y.dist==min(y.dist))
	summary$kd.value[n] <- kd$z[x.i,y.i]
	}
#------------------------------------------------------------------------------------------------------------
# plots
#------------------------------------------------------------------------------------------------------------
date <- summary$meanBP
lat <- summary$Latitude
X <- c(11500,13700,14000)
Y <- c(31,29,28,26)
kd.pal <- readLines('../tools/colours.txt')
#------------------------------------------------------------------------------------------------------------
# WDI mean

stat <- summary$WDI.mean
cols <- summary$WDI.mean.col
cor <- cor(stat, summary$kd.value)
file <- '../plots/KD vs WDI.mean.pdf'
txt <- 'Animal WDI (mean)'
stat.pal <- readLines('../tools/colours.txt')
KEplot(kd, hyd, stat, date, lat, cols, cor, file, X, Y, txt, stat.pal, kd.pal)
#------------------------------------------------------------------------------------------------------------
# WDI min

#stat <- summary$WDI.min
#cols <- summary$WDI.min.col
#cor <- cor(stat, summary$kd.value)
#file <- '../plots/KD vs WDI.min.pdf'
#txt <- 'Animal WDI (min)'
#stat.pal <- readLines('../tools/colours.txt')
#KEplot(kd, hyd, stat, date, lat, cols, cor, file, X, Y, txt, stat.pal, kd.pal)
#------------------------------------------------------------------------------------------------------------
# Richness

#stat <- summary$Richness
#cols <- summary$Richness.col
#cor <- cor(stat, summary$kd.value)
#file <- '../plots/KD vs Richness.pdf'
#txt <- 'Species Richness'
#stat.pal <- rev(readLines('../tools/colours.txt'))
#KEplot(kd, hyd, stat, date, lat, cols, cor, file, X, Y, txt, stat.pal, kd.pal)
#------------------------------------------------------------------------------------------------------------
# Shannon

#stat <- summary$Shannon
#cols <- summary$Shannon.col
#cor <- cor(stat, summary$kd.value)
#file <- '../plots/KD vs Shannon.pdf'
#txt <- 'Species Shannon'
#stat.pal <- rev(readLines('../tools/colours.txt'))
#KEplot(kd, hyd, stat, date, lat, cols, cor, file, X, Y, txt, stat.pal, kd.pal)
#------------------------------------------------------------------------------------------------------------
