#------------------------------------------------------------------------------------------------------------
# Explore other aspects of hydrology
#------------------------------------------------------------------------------------------------------------
load('hyd.cal.dates.RData')

# investigate lake levels for those with multiple samples for the same lake
d <- hyd[!is.na(hyd$LakeLevel),]

d$col[d$LakeLevel=='LOW'] <- 'orange'
d$col[d$LakeLevel=='MED'] <- 'green'
d$col[d$LakeLevel=='HIGH'] <- 'blue'

plot(d$wm.date, d$Latitude, pch=16, col=d$col)

rgl::plot3d(d$Latitude,  d$Longitude, d$wm.date, col=d$col, size=5)
#------------------------------------------------------------------------------------------------------------
unique(subset(d, Latitude>29)$SiteID)