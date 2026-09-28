#--------------------------------------------------------------------
# All faunal analysis and Fig 3
#--------------------------------------------------------------------
data <- read.csv(file='../tools/faunal.csv')
source('functions.R')
lat.posts <- seq(13,37,by=2)
lat.mid <- seq(14,36,by=2)
cols <- readLines('../tools/colours.txt')
#--------------------------------------------------------------------
# generate proportions of each WDI in each lat slice
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
# pre-process presence only entries for weighting the mean WDI code by NISP
# get the distribution of NISP proportions for all species that are also presenceOnly, 
# and calculate various quantiles. These can then map to rare, present, frequent, very frequent

PO.taxa <- sort(unique(data[is.na(data$PresenceOnly),]$TaxonCode))
N <- length(PO.taxa)
quantiles <- matrix(,N,4)
for(n in 1:N){
	phases <- unique(subset(data, TaxonCode==PO.taxa[n])$PhaseID)
	P <- length(phases)
	dist <- numeric(P)
	for(p in 1:P){
		sub <- subset(data, PhaseID==phases[p])
		dist[p] <- (sub$NISP/sum(sub$NISP,na.rm=TRUE))[sub$TaxonCode==PO.taxa[n]]
		}
	quantiles[n,] <- as.numeric(quantile(dist,na.rm=TRUE, prob=c(0.125,0.375,0.625,0.875)))
	}
PO.weights <- data.frame(taxa=PO.taxa, rare=quantiles[,1], present=quantiles[,2], frequent=quantiles[,3],very.frequent=quantiles[,4])
#--------------------------------------------------------------------
# generate phase summary of WDI
phases <- unique(data$PhaseID)

N <- length(phases)
names <- c('SiteID','PhaseID','StartBP','EndBP','Latitude','Longitude','WDI.mean','WDI.min','WDI.max','Richness')
summary <- as.data.frame(matrix(,N,length(names))) 
names(summary) <- names
for(n in 1:N){
	sub <- subset(data,PhaseID==phases[n])

	# mean WDI weighted by NISP proportions, including pre-processed proportions for presence data
	nisp.prop <- sub$NISP/sum(sub$NISP, na.rm=TRUE)

	# get preprocessed for taxa missing NISP
	PO.taxa <- sub$TaxonCode[is.na(sub$NISP)]
	T <- length(PO.taxa)
	if(T>0){
		for(t in 1:T){
			taxon <- PO.taxa[t]
			PO <- PO.weights[PO.weights$taxa==taxon,]
			if(nrow(PO)==0)PO[1,] <- c(taxon, 1,1,1,1)
			level <- sub$PresenceOnly[sub$TaxonCode==taxon]
			prop <- as.numeric(PO[,names(PO)==level])
			nisp.prop[sub$TaxonCode==taxon] <- prop
			}
		}
	nisp.prop <- nisp.prop/sum(nisp.prop)
			
	summary$SiteID[n] <- unique(sub$SiteID)
	summary$PhaseID[n] <- unique(sub$PhaseID)
	summary$StartBP[n] <- unique(sub$StartBP)
	summary$EndBP[n] <- unique(sub$EndBP)
	summary$Latitude[n] <- unique(sub$Latitude)
	summary$Longitude[n] <- unique(sub$Longitude)
	summary$WDI.mean[n] <- sum(sub$WaterCode * nisp.prop)
	summary$WDI.min[n] <- min(as.numeric(sub$WaterCode),na.rm=TRUE) 
	summary$WDI.max[n] <- max(as.numeric(sub$WaterCode),na.rm=TRUE) 
	summary$Richness[n] <- nrow(sub)
	}
summary$meanBP <- (summary$StartBP + summary$EndBP)/2
summary <- subset(summary, !is.na(meanBP))

# generate continuous colours for WDI.mean
summary$WDI.mean.col <- convert.stat.to.colours(summary$WDI.mean, cols, v=seq(1,5,length.out=20))

# save summary for other purposes
write.csv(summary, file='../tools/faunal phases.csv', row.names=FALSE)
#--------------------------------------------------------------------
# plot
#--------------------------------------------------------------------
pdf(file='../plots/fig 3 faunal.pdf', width=14, height=5)

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
dev.off()
#--------------------------------------------------------------------



