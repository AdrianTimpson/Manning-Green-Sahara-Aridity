#--------------------------------------------------------------------
# biodiversity metrics
#--------------------------------------------------------------------
data <- read.csv('../tools/animals.summary.csv')
cols <- readLines('../tools/colours.txt')

#pdf(file='../plots/diversity stats.pdf',width=14,height=6)
par(mfrow=c(1,2))
plot(Richness ~ Latitude, data=data, pch=16, col=data$WDI.mean.col, cex=1.5)
plot(Shannon ~ Latitude, data=data, pch=16, col=data$WDI.mean.col, cex=1.5)
legend(x=32,y=3,bty='n', legend=paste('WDI =',1:5), pch=16, col=cols, cex=0.8)
#dev.off()
#--------------------------------------------------------------------

