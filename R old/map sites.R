#------------------------------------------------------------------------------------------------------------
# overheads
#------------------------------------------------------------------------------------------------------------
source('functions.R')
source('get hydrology.R')
source('get animals.R')
#------------------------------------------------------------------------------------------------------------
# reduce to unique sites
hydro <- unique(d[,c('site','Latitude','Longitude')]) # 515
animal <- unique(summary[,c('SiteID','Latitude','Longitude')]) # 85
#------------------------------------------------------------------------------------------------------------
# create a kml of points, open in google earth, edit icons, save and crop
#------------------------------------------------------------------------------------------------------------
hydro <- sf::st_as_sf(hydro, coords = c('Longitude','Latitude'), crs=4326)
sf::st_write(hydro, "../shapes/hydrologically_active_sites.kml", layer='hydrology', driver="KML") 
tmp <- readLines("../shapes/hydrologically_active_sites.kml")
tmp <- gsub('<name>.*</name>','<name></name>',tmp)
writeLines(tmp, con="../shapes/hydrologically_active_sites.kml")

animal <- sf::st_as_sf(animal, coords = c('Longitude','Latitude'), crs=4326)
sf::st_write(animal, "../shapes/animal_sites.kml", layer='animals', driver="KML") 
tmp <- readLines("../shapes/animal_sites.kml")
tmp <- gsub('<name>.*</name>','<name></name>',tmp)
writeLines(tmp, con="../shapes/animal_sites.kml")
#------------------------------------------------------------------------------------------------------------






