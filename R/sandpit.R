#-----------------------------------------------------------------------------------
# Investigate some issues. 
# Why is the below object now bigger than the one I generate in the 'generate data' script?
# at what point and why is there anything stored as .RData files?

# OK I get it now. Latest is correct.
#-----------------------------------------------------------------------------------
# get 14C hydrology dates
load('R objects/hyd.cal.dates.RData')
hyd1 <- hyd

e1 <- hyd1[!hyd1$SiteID%in%hyd$SiteID,]
e <- hyd[!hyd$SiteID%in%hyd1$SiteID,]

subset(hyd1, SiteID%in%c('S1016'))


#-----------------------------------------------------------------------------------
# map them all
require(maps)
xlim <- c(-20,45)
ylim <- c(10,35)
plot(NULL, xlim=xlim, ylim=ylim, frame.plot=FALSE, axes=FALSE, xlab='', ylab='', main='all hydro 14C dates')
map('world', xlim=xlim, ylim=ylim, col='grey90', add=TRUE, fill=TRUE, border='grey')
points(hyd$Longitude, hyd$Latitude, pch=16)
points(e1$Longitude, e1$Latitude, pch=16,col='red')
#-----------------------------------------------------------------------------------
# summary
dim(hyd)
hist(hyd$wm.date, breaks=50)
range(hyd$wm.date)
#-----------------------------------------------------------------------------------


