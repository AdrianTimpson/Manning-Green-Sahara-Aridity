#--------------------------------------------------------------------
# All faunal analysis and Fig 3
#--------------------------------------------------------------------
data <- read.csv(file='../tools/faunal.csv')
source('functions.R')
lat.posts <- seq(13,37,by=2)
lat.mid <- seq(14,36,by=2)
cols <- readLines('../tools/colours.txt')
#--------------------------------------------------------------------
# generate proportions of each WDI in each phase
N <- length(lat.posts)-1
N.phases <- N.spec <- numeric(N)
lat.names <- c(); for(n in 1:N)lat.names[n] <- paste(lat.posts[n],lat.posts[n+1],sep='-')
wc <- as.data.frame(matrix(,N,5)); names(wc) <- c('WDI.1','WDI.2','WDI.3','WDI.4','WDI.5')
for(n in 1:N){
	slice <- subset(data, Latitude>=lat.posts[n] & Latitude<=lat.posts[n+1])
	total.phases <- length(unique(slice$PhaseID))
	N.phases[n] <- total.phases
	N.spec[n] <- length(unique(slice$TaxonCode))
	for(code in 1:5){
		n.phases <- length(unique(subset(slice, WaterCode==code)$PhaseID))
		wc[n,code] <- n.phases/total.phases
		}
	}
row.names(wc) <- lat.names
#--------------------------------------------------------------------
############
# problem with duplications. Katie fixing it
N <- length(phases)
duplicates <- c()
for(n in 1:N){
	sub <- subset(data,PhaseID==phases[n])
	dup <- unique(sub$TaxonCode[duplicated(sub$TaxonCode)])
	if(length(dup)>0)duplicates <- c(duplicates, paste(phases[n], dup))
	}
duplicates

# pre-process presence only entries-
# needed for weighting the WDI code by NISP
several issues - if all species have the same level (e.g. rare) and there are no NISP, they can be weighted equally etc
table(data$PresenceOnly)

#################
#--------------------------------------------------------------------
# generate phase summary of WDI
phases <- unique(data$PhaseID)

N <- length(phases)
names <- c('SiteID','PhaseID','StartBP','EndBP','Latitude','Longitude','WDI.mean','WDI.min','WDI.max','Richness')
summary <- as.data.frame(matrix(,N,length(names))) 
names(summary) <- names
for(n in 1:N){
	sub <- subset(data,PhaseID==phases[n])
	summary$SiteID[n] <- unique(sub$SiteID)
	summary$PhaseID[n] <- unique(sub$PhaseID)
	summary$StartBP[n] <- unique(sub$StartBP)
	summary$EndBP[n] <- unique(sub$EndBP)
	summary$Latitude[n] <- unique(sub$Latitude)
	summary$Longitude[n] <- unique(sub$Longitude)
	summary$WDI.mean[n] <- mean(as.numeric(sub$WaterCode),na.rm=TRUE) # should this be weighted by NISP?? 
#	summary$WDI.mean[n] <- sum(sub$WaterCode * log(sub$NISP))/sum(log(sub$NISP)) # should this be weighted by NISP?? 
	summary$WDI.min[n] <- min(as.numeric(sub$WaterCode),na.rm=TRUE) 
	summary$WDI.max[n] <- max(as.numeric(sub$WaterCode),na.rm=TRUE) 
	summary$Richness[n] <- nrow(sub)
	}
summary$meanBP <- (summary$StartBP + summary$EndBP)/2
summary <- subset(summary, !is.na(meanBP))

# generate continuous colours for WDI.mean
summary$WDI.mean.col <- convert.stat.to.colours(summary$WDI.mean, cols, v=seq(1,5,length.out=20))
#--------------------------------------------------------------------
# plot
#--------------------------------------------------------------------

#pdf(file='../plots/WDI circles.pdf', width=14, height=7)

N <- length(lat.posts)-1
x.text <- c(6,7)
par(mar=c(5,3,1,1))
plot(NULL, xlim=c(0.5,11), ylim=c(13,39), bty='n', xlab=paste('WDI',paste(rep(' ',120),collapse=''),'Richness'), ylab='Latitude',xaxt='n',yaxt='n')
axis(1, 1:5)
axis(2, lat.posts, las=1)
for(n in 1:(N+1))lines(x=c(0,11),y=rep(lat.posts[n],2),col='grey',lty=2)
abline(v=mean(x.text), col='grey', lty=2)
text(x=x.text[1],y=38,'N phases')
text(x=x.text[2],y=38,'N species')
scale <- 4
for(n in 1:N){
	points(x=1:5, y=rep(lat.mid[n],5), pch=16, cex=as.numeric(wc[n,])*scale,col=cols)
	text(x.text[1],y=lat.mid[n], N.phases[n])
	text(x.text[2],y=lat.mid[n], N.spec[n])
	}

x.shift <- 7.5
x.scale <- 0.08
points(x.scale*summary$Richness+x.shift, summary$Latitude, pch=16, col=summary$WDI.mean.col, cex=1.1)
xticks <- seq(0,45,by=5)
axis(1, at=x.scale*xticks+x.shift, labels=xticks, cex.axis=0.8)

legend(bty='n',x=0,y=40,pch=16,pt.cex=c(0.1,0.3,0.5,1)*scale,legend=c('10%','30%','50%','100%'),col='grey',y.intersp=2, horiz=TRUE)
#dev.off()
#--------------------------------------------------------------------



