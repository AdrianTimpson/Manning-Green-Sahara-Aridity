#--------------------------------------------------------------------
# 3D barplots
#--------------------------------------------------------------------
data <- read.csv('../tools/animals.raw.csv')
range(data$Latitude)
length(unique(data$PhaseID))
sum(data$NISP)
#--------------------------------------------------------------------
# generate summary data
#--------------------------------------------------------------------
# choose some posts
lat.posts <- seq(15,35,by=2)

N <- length(lat.posts)-1
lat.names <- c(); for(n in 1:N)lat.names[n] <- paste(lat.posts[n],lat.posts[n+1],sep='-')
wc <- as.data.frame(matrix(,N,5)); names(wc) <- c('WDI.1','WDI.2','WDI.3','WDI.4','WDI.5')
for(n in 1:N){
	slice <- subset(data, Latitude>=lat.posts[n] & Latitude<=lat.posts[n+1])
	total.phases <- length(unique(slice$PhaseID))
	for(code in 1:5){
		n.phases <- length(unique(subset(slice, WaterCode==code)$PhaseID))
		wc[n,code] <- n.phases/total.phases
		}
	}
row.names(wc) <- lat.names
#--------------------------------------------------------------------
# plot
#--------------------------------------------------------------------
cols <- readLines('../tools/colours.txt')
colmat <- matrix(rep(1:5,each=nrow(wc)),nrow(wc),ncol(wc))

library(plot3D)
x <- 1:nrow(wc)
y <- 1:ncol(wc)
z <- as.matrix(wc)

pdf(file='../plots/barchart 3D.pdf', width=14, height=7)
par(mar=rep(0.5,4))
hist3D(	x=x,y=y,z=z*5, col=cols, colvar=colmat, border='grey', 
		space=c(0.2,0.4), shade=0.3, alpha=0.5,
		scale=FALSE, xlab='', ylab='', axes=FALSE, colkey = list(plot = FALSE),
		theta=45, phi=30, r=1)

text3D(x=1:N-0.2, y=rep(-0.15,N),  z=rep(0,N), labels=row.names(wc), add=TRUE, cex=0.7, srt=20)
text3D(x=rep(N+0.5,5), y=1:5 +0.2,  z=rep(-0.2,5), labels=names(wc), add=TRUE, cex=0.7, srt=-30)
text3D(x=rep(0,6), y=rep(0.3,6),  z=seq(0,1,by=0.2)*5, labels=seq(0,100,by=20), add=TRUE, cex=0.7, srt=0)

text3D(x=5, y=-0.5, z=0, labels='Latitude',add=TRUE, srt=-40, cex=0.8)
text3D(x=0, y=-0.1, z=0.8, labels='% of site phases',add=TRUE, srt=100, cex=0.8)
text3D(x=N+1.2, y=1.5, z=0, labels='Water Dependency Index',add=TRUE, srt=45, cex=0.8)
dev.off()
#--------------------------------------------------------------------





