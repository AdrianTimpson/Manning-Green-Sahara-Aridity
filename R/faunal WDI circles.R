#--------------------------------------------------------------------
# calculate proportion of phases that contain a WDI
#--------------------------------------------------------------------
data <- read.csv('../tools/animals.raw.csv')
lat.posts <- seq(13,37,by=2)
lat.mid <- seq(14,36,by=2)
cols <- readLines('../tools/colours.txt')

# generate summary data
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
# plot
#--------------------------------------------------------------------
pdf(file='../plots/WDI circles.pdf', width=14, height=7)
x.text <- c(6,7)
plot(NULL, xlim=c(0.5,9), ylim=c(13,38), bty='n', xlab='WDI', ylab='Latitude',xaxt='n',yaxt='n')
axis(1, 1:5)
axis(2, lat.posts, las=1)
for(n in 1:(N+1))lines(x=c(0,7.5),y=rep(lat.posts[n],2),col='grey',lty=2)
abline(v=mean(x.text), col='grey', lty=2)
text(x=x.text[1],y=38,'N phases')
text(x=x.text[2],y=38,'N species')
scale <- 4
for(n in 1:N){
	points(x=1:5, y=rep(lat.mid[n],5), pch=16, cex=as.numeric(wc[n,])*scale,col=cols)
	text(x.text[1],y=lat.mid[n], N.phases[n])
	text(x.text[2],y=lat.mid[n], N.spec[n])
	}
legend(bty='n',x=8,y=37,pch=16,pt.cex=c(0.1,0.3,0.5,1)*scale,legend=c('10%','30%','50%','100%'),col='grey',y.intersp=2)
dev.off()
#--------------------------------------------------------------------










range(data$StartBP, na.rm=TRUE)
range(data$EndBP, na.rm=TRUE)
range(data$Latitude, na.rm=TRUE)
length(unique(data$PhaseID))
length(unique(data$SiteID))
sum(data$NISP,na.rm=TRUE)
#--------------------------------------------------------------------
# generate summary data
#--------------------------------------------------------------------
# choose some posts
lat.posts <- seq(13,37,by=2)
cols <- readLines('../tools/colours.txt')

N <- length(lat.posts)-1
lat.names <- c(); for(n in 1:N)lat.names[n] <- paste(lat.posts[n],lat.posts[n+1],sep='-')
wc <- as.data.frame(matrix(,N,5)); names(wc) <- c('WDI.1','WDI.2','WDI.3','WDI.4','WDI.5')
for(n in 1:N){
	slice <- subset(data, Latitude>=lat.posts[n] & Latitude<=lat.posts[n+1])
	total.phases <- length(unique(slice$PhaseID))
	print(total.phases)
	for(code in 1:5){
		n.phases <- length(unique(subset(slice, WaterCode==code)$PhaseID))
		wc[n,code] <- n.phases/total.phases
		}
	}
row.names(wc) <- lat.names
#--------------------------------------------------------------------
# plot
#--------------------------------------------------------------------
cex <- 1.5
plot(wc$WDI.1,pch=16, col=cols[1],type='b', ylim=c(0,1),cex=cex,xaxt='n',xlab='Latitude (degrees N)',ylab='% of site phases')
axis(side=1, at = 1:nrow(wc), labels = row.names(wc))
lines(wc$WDI.2, pch=16, col=cols[2], type='b',cex=cex)
lines(wc$WDI.3, pch=16, col=cols[3], type='b',cex=cex)
lines(wc$WDI.4, pch=16, col=cols[4], type='b',cex=cex)
lines(wc$WDI.5, pch=16, col=cols[5], type='b',cex=cex)
#--------------------------------------------------------------------

library(sf)
require(maps)
source('functions.R')

d1 <- read.csv('../tools/animals.summary.csv')

plot(NULL, xlim=plot.xlim, ylim=plot.ylim, bty='n', xaxs='i',yaxs='i',axes=TRUE, xlab='', ylab='',main='')
polygon(x=c(plot.xlim,rev(plot.xlim)), y=rep(plot.ylim,each=2),col='steelblue',border='steelblue')
map('world', xlim=plot.xlim, ylim=plot.ylim, col='grey90', add=TRUE, fill=TRUE, border='grey90')
points(d1$Longitude+rnorm(nrow(d1),0,0.5), d1$Latitude+rnorm(nrow(d1),0,0.5),pch=16, col=d1$WDI.mean.col)

# assign data top time slices
quantile(d1$meanBP)
breaks <- seq(2000,14000,by=500)
d1$timebin <- cut(d1$meanBP, breaks, 1:(length(breaks)-1))
N <- length(breaks)-1
slice.names <- c()
for(n in 1:N)slice.names[n] <- paste0(breaks[n]/1000, '-', breaks[n+1]/1000, 'kyrBP')

# plot region
plot.xlim <- c(-20,45)
plot.ylim <- c(13,40)

par(mfrow=c(4,6),mar=c(2,2,3,1))
for(n in 1:N){
	plot(NULL, xlim=plot.xlim, ylim=plot.ylim, bty='n', xaxs='i',yaxs='i',axes=TRUE, xlab='', ylab='',main=slice.names[n])
	polygon(x=c(plot.xlim,rev(plot.xlim)), y=rep(plot.ylim,each=2),col='steelblue',border='steelblue')
	map('world', xlim=plot.xlim, ylim=plot.ylim, col='grey90', add=TRUE, fill=TRUE, border='grey90')
	sub <- subset(d1, timebin==n)
	points(sub$Longitude, sub$Latitude,pch=16, col=sub$WDI.mean.col)
	}





