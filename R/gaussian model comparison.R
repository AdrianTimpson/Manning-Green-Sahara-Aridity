#-----------------------------------------------------------------------------------
# a model of the rainfall as a north band moving down, and a south ban moving up
# should be nested for proper model comparison:
# always a mixture of two gaussians, such that both bands:
# have their own angle (from horizontal)
# have their own SD
# have their own starting latitude (at greenwich meridian)
# have their own advance/decline rate (+ve or -ve). 

# simpler models nested within:
# fix angle to zero (shifts are purely north/south)
# both bands have the same (inverse) advance/decline rate
# both bands have the same SD

# mS: 0-1 mixture proportion of the southern band (mN must be 1-mS)
# aS, aN: angle of southern band and northern band from horizontal, in degrees
# sdS, sdN: standard deviation of each band, in degrees
# yS, yN: latitude of each band in degrees, at long=0 (geenwich meridian) at start time (14000yBP)
# rS, rN: advance/decline rate of each band, degrees per year
# revD: reverse direction date
#-----------------------------------------------------------------------------------
# get hydrological data
load('hyd.cal.dates.RData')

# data should be lat, long, date(BP)
data <- hyd[,c('Latitude','Longitude','wm.date')]; names(data)[3] <- 'date'

require(DEoptimR)
x <- data$Longitude
y <- data$Latitude
#-----------------------------------------------------------------------------------
point.dist.from.origin <- function(x,y,theta,y0){
	dist <- cos(theta*pi/180)*(x*tan(theta*pi/180) + y-y0)
return(dist)}

gaussian.dist.from.origin <- function(date, revD, rate){	
	dist <- numeric(length(date))
	i <- date>=revD
	dist[i] <- (14000-date[i])*rate
	dist[!i] <- (14000-2*revD+date[!i])*rate
return(dist)}

lik <- function(x, y, date, mixS, aS, aN, sdS, sdN, yS, yN, rS, rN, revD){
	PDFO.S <- point.dist.from.origin(x=x,y=y,theta=aS,y0=yS) 
	PDFO.N <- point.dist.from.origin(x=x,y=y,theta=aN,y0=yN) 
	GDFO.S <- gaussian.dist.from.origin(date, revD, rS)
	GDFO.N <- gaussian.dist.from.origin(date, revD, rN)	
	dN <- PDFO.N - GDFO.N
	dS <- PDFO.S - GDFO.S
	res <- dnorm(dN, 0, sdN)*(1-mixS) + dnorm(dS, 0, sdS)*mixS
return(res)}

liksingle <- function(x, y, date, aS, sdS, yS, rS, revD){
	PDFO.S <- point.dist.from.origin(x=x,y=y,theta=aS,y0=yS) 
	GDFO.S <- gaussian.dist.from.origin(date, revD, rS)
	dS <- PDFO.S - GDFO.S
	res <- dnorm(dS, 0, sdS)
return(res)}
#-----------------------------------------------------------------------------------
# maximum likelihood searches
#-----------------------------------------------------------------------------------
# model 10
obj.fnc <- function(pars, x, y, date){
		mixS <- pars[1]
		aS <- pars[2]
		aN <- pars[3]
		sdS <- pars[4]
		sdN <- pars[5]
		yS <- pars[6]
		yN <- pars[7]
		rS <- pars[8]
		rN <- pars[9]
		revD <- pars[10]
	lik <- lik(x, y, date, mixS, aS, aN, sdS, sdN, yS, yN, rS, rN, revD)
return(-sum(log(lik)))}

lower <- c(0,-35,-35, 2, 2, 0,  20, -0.01, -0.01, 5000)
upper <- c(1, 35,  0, 5, 5, 20, 50,  0.01,  0.01, 12000)
m10 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, trace=TRUE, tol=1e-10, NP=300)

# model 9a
obj.fnc <- function(pars, x, y, date){
		mixS <- pars[1]
		aS <- pars[2]
		aN <- pars[3]
		sdS <- pars[4]
		sdN <- pars[5]
		yS <- pars[6]
		yN <- pars[7]
		rS <- pars[8]
		rN <- -pars[8]
		revD <- pars[9]
	lik <- lik(x, y, date, mixS, aS, aN, sdS, sdN, yS, yN, rS, rN, revD)
return(-sum(log(lik)))}

lower <- c(0,-35,-35, 2, 2,  0, 20, -0.01, 5000)
upper <- c(1, 35,  0, 5, 5, 20, 50,  0.01, 12000)
m9a <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, trace=TRUE, tol=1e-10, NP=300)

# model 9b
obj.fnc <- function(pars, x, y, date){
		mixS <- pars[1]
		aS <- pars[2]
		aN <- pars[3]
		sdS <- pars[4]
		sdN <- pars[4]
		yS <- pars[5]
		yN <- pars[6]
		rS <- pars[7]
		rN <- pars[8]
		revD <- pars[9]
	lik <- lik(x, y, date, mixS, aS, aN, sdS, sdN, yS, yN, rS, rN, revD)
return(-sum(log(lik)))}

lower <- c(0,-35,-35, 2,  0, 20, -0.01, -0.01, 5000)
upper <- c(1, 35,  0, 5, 20, 50,  0.01,  0.01, 12000)
m9b <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, trace=TRUE, tol=1e-10, NP=300)

# model 8
obj.fnc <- function(pars, x, y, date){
		mixS <- pars[1]
		aS <- pars[2]
		aN <- pars[3]
		sdS <- pars[4]
		sdN <- pars[4]
		yS <- pars[5]
		yN <- pars[6]
		rS <- pars[7]
		rN <- -pars[7]
		revD <- pars[8]
	lik <- lik(x, y, date, mixS, aS, aN, sdS, sdN, yS, yN, rS, rN, revD)
return(-sum(log(lik)))}

lower <- c(0,-35,-35, 2,  0, 20, -0.01, 5000)
upper <- c(1, 35,  0, 5, 20, 50,  0.01, 12000)
m8 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, trace=TRUE, tol=1e-10, NP=300)

# model 7
obj.fnc <- function(pars, x, y, date){
		mixS <- pars[1]
		aS <- 0
		aN <- pars[2]
		sdS <- pars[3]
		sdN <- pars[3]
		yS <- pars[4]
		yN <- pars[5]
		rS <- pars[6]
		rN <- -pars[6]
		revD <- pars[7]
	lik <- lik(x, y, date, mixS, aS, aN, sdS, sdN, yS, yN, rS, rN, revD)
return(-sum(log(lik)))}

lower <- c(0,-35, 2,  0, 20, -0.01, 5000)
upper <- c(1,  0, 5, 20, 50,  0.01, 12000)
m7 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, trace=TRUE, tol=1e-10, NP=300)

# model 6
obj.fnc <- function(pars, x, y, date){
		mixS <- pars[1]
		aS <- 0
		aN <- 0
		sdS <- pars[2]
		sdN <- pars[2]
		yS <- pars[3]
		yN <- pars[4]
		rS <- pars[5]
		rN <- -pars[5]
		revD <- pars[6]
	lik <- lik(x, y, date, mixS, aS, aN, sdS, sdN, yS, yN, rS, rN, revD)
return(-sum(log(lik)))}

lower <- c(0, 2,  0, 20, -0.01, 5000)
upper <- c(1, 5, 20, 50,  0.01, 12000)
m6 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, trace=TRUE, tol=1e-10, NP=300)

# model 5
obj.fnc <- function(pars, x, y, date){
		aS <- pars[1]
		sdS <- pars[2]
		yS <- pars[3]
		rS <- pars[4]
		revD <- pars[5]
	lik <- liksingle(x, y, date, aS, sdS, yS, rS, revD)
return(-sum(log(lik)))}

lower <- c(-35, 2, 0, -0.01, 5000)
upper <- c(35,  5, 20, 0.01, 12000)
m5 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, trace=TRUE, tol=1e-10, NP=300)

# model 4
obj.fnc <- function(pars, x, y, date){
		aS <- 0
		sdS <- pars[1]
		yS <- pars[2]
		rS <- pars[3]
		revD <- pars[4]
	lik <- liksingle(x, y, date, aS, sdS, yS, rS, revD)
return(-sum(log(lik)))}

lower <- c(2, 0, -0.01, 5000)
upper <- c(5, 20, 0.01, 12000)
m4 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, trace=TRUE, tol=1e-10, NP=300)

#-----------------------------------------------------------------------------------
# save(m4,m5,m6,m7,m8,m9a,m9b,m10, file='model results.RData')
load('model results.RData')
round(m10$par,2)
#-----------------------------------------------------------------------------------
# Model comparison
#-----------------------------------------------------------------------------------
# AIC
(2*10)-(2*-m10$value)
(2*9)-(2*-m9a$value)
(2*9)-(2*-m9b$value)
(2*8)-(2*-m8$value)
(2*7)-(2*-m7$value)
(2*6)-(2*-m6$value)
(2*5)-(2*-m5$value)
(2*4)-(2*-m4$value)

# BIC
N <- nrow(data)
(2*10)*log(N)-(2*-m10$value)
(2*9)*log(N)-(2*-m9a$value)
(2*9)*log(N)-(2*-m9b$value)
(2*8)*log(N)-(2*-m8$value)
(2*7)*log(N)-(2*-m7$value)
(2*6)*log(N)-(2*-m6$value)
(2*5)*log(N)-(2*-m5$value)
(2*4)*log(N)-(2*-m4$value)
#-----------------------------------------------------------------------------------
# plot best model at various time points
#-----------------------------------------------------------------------------------
p <- m6$par # p <- m10$par

# at 14000BP
N <- 20000
xlim <- c(-20,45)
ylim <- c(10,37)
x <- runif(N, min(xlim),max(xlim))
y <- runif(N, min(ylim),max(ylim))

breaks <- seq(14000,2000,length.out=10)
mids <- round((breaks[1:9]+breaks[2:10])/2,-2)

par(mfrow=c(3,3))
for(b in 1:9){
	date <- rep(mids[b],N)
	#liks <- lik(x, y, date, mixS=p[1], aS=0, aN=0, sdS=p[2], sdN=p[2], yS=p[3], yN=p[4], rS=p[5], rN=-p[5], revD=p[6]) #for model 6
	liks <- lik(x, y, date, mixS=p[1], aS=p[2], aN=p[3], sdS=p[4], sdN=p[5], yS=p[6], yN=p[7], rS=p[8], rN=-p[9], revD=p[10]) #for model 10
	pal <- colorRampPalette(c('red','yellow','green','blue'))
	order <- findInterval(liks, sort(liks))
	plot(x,y,col=pal(N)[order],pch=16, main=paste(mids[b],'BP'))
	maps::map('world', xlim=xlim, ylim=ylim, col=NA, add=TRUE, fill=TRUE, border='black')
	sub <- subset(data, date<=breaks[b] & date>=breaks[b+1])
	points(sub$Longitude, sub$Latitude)
	}
#-----------------------------------------------------------------------------------

N <- 100
x <- rep(0,N)
par(mfrow=c(1,2))
y <- seq(min(ylim),max(ylim),length.out=N)
dates <- round(seq(14000,p[6],length.out=5),-2) 
plot(NULL, xlim=ylim, ylim=c(0,0.12),xlab='latitude', ylab='model 6 probability density')
for(n in 1:5){
	date <- rep(dates[n],N)
	liks <- lik(x, y, date, mixS=p[1], aS=0, aN=0, sdS=p[2], sdN=p[2], yS=p[3], yN=p[4], rS=p[5], rN=-p[5], revD=p[6])
	lines(y,liks,type='l',col=n)
	}
legend(x=30,y=0.1, col=1:5, lwd=2, legend=paste(dates,'yBP'),bty='n', cex=0.7)

dates <- round(seq(p[6],2000,length.out=5),-2) 
plot(NULL, xlim=ylim, ylim=c(0,0.12),xlab='latitude', ylab='model 6 probability density')
for(n in 1:5){
	date <- rep(dates[n],N)
	liks <- lik(x, y, date, mixS=p[1], aS=0, aN=0, sdS=p[2], sdN=p[2], yS=p[3], yN=p[4], rS=p[5], rN=-p[5], revD=p[6])
	lines(y,liks,type='l',col=n)
	}
legend(x=30,y=0.1, col=1:5, lwd=2, legend=paste(dates,'yBP'),bty='n', cex=0.7)









