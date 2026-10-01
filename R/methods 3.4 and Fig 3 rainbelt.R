#-----------------------------------------------------------------------------------
# a model of the rainfall as a north band moving down, and a south band moving up
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
# only a single rainband from the south

# mS: 0-1 mixture proportion of the southern band (mN must be 1-mS)
# aS, aN: angle of southern band and northern band from vertical, in degrees
# sdS, sdN: standard deviation of each band, in degrees
# yS, yN: latitude of each band in degrees, at long=0 (geenwich meridian) at start time (14000yBP)
# rS, rN: advance/decline rate of each band, degrees per year
# revD: reverse direction date
#-----------------------------------------------------------------------------------
# get hydrological data
load('R objects/hyd.cal.dates.RData')

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

lower <- c(0,  0,-35, 1, 1, 0,  20, -0.01, -0.01, 5000)
upper <- c(1, 35,  0, 7, 7, 20, 50,  0.01,  0.01, 12000)
m10 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date,  tol=1e-10, NP=300)

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

lower <- c(0,  0,-35, 1, 1,  0, 20, -0.01, 5000)
upper <- c(1, 35,  0, 7, 7, 20, 50,  0.01, 12000)
m9a <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, tol=1e-10, NP=300)

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

lower <- c(0,  0,-35, 1,  0, 20, -0.01, -0.01, 5000)
upper <- c(1, 35,  0, 7, 20, 50,  0.01,  0.01, 12000)
m9b <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, tol=1e-10, NP=300)

# model 9c
obj.fnc <- function(pars, x, y, date){
		mixS <- pars[1]
		aS <- 0
		aN <- pars[2]
		sdS <- pars[3]
		sdN <- pars[4]
		yS <- pars[5]
		yN <- pars[6]
		rS <- pars[7]
		rN <- pars[8]
		revD <- pars[9]
	lik <- lik(x, y, date, mixS, aS, aN, sdS, sdN, yS, yN, rS, rN, revD)
return(-sum(log(lik)))}

lower <- c(0,-35, 1, 1, 0,  20, -0.01, -0.01, 5000)
upper <- c(1,  0, 7, 7, 20, 50,  0.01,  0.01, 12000)
m9c <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date,  tol=1e-10, NP=300)

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

lower <- c(0,  0,-35, 1,  0, 20, -0.01, 5000)
upper <- c(1, 35,  0, 7, 20, 50,  0.01, 12000)
m8 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, NP=300)

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

lower <- c(0,-35, 1,  0, 20, -0.01, 5000)
upper <- c(1,  0, 7, 20, 50,  0.01, 12000)
m7 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, NP=300)

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

lower <- c(0, 1,  0, 20, -0.01, 5000)
upper <- c(1, 7, 20, 50,  0.01, 12000)
m6 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, NP=300)

# model 5
obj.fnc <- function(pars, x, y, date){
		aS <- pars[1]
		sdS <- pars[2]
		yS <- pars[3]
		rS <- pars[4]
		revD <- pars[5]
	lik <- liksingle(x, y, date, aS, sdS, yS, rS, revD)
return(-sum(log(lik)))}

lower <- c(0, 1, 0, -0.01, 5000)
upper <- c(35, 7, 20, 0.01, 12000)
m5 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, NP=300)

# model 4
obj.fnc <- function(pars, x, y, date){
		aS <- 0
		sdS <- pars[1]
		yS <- pars[2]
		rS <- pars[3]
		revD <- pars[4]
	lik <- liksingle(x, y, date, aS, sdS, yS, rS, revD)
return(-sum(log(lik)))}

lower <- c(1, 0, -0.01, 5000)
upper <- c(7, 20, 0.01, 12000)
m4 <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, x=x, y=y, date=data$date, NP=300)

#-----------------------------------------------------------------------------------
#save(m4,m5,m6,m7,m8,m9a,m9b,m9c,m10, file='R objects/model comparison results.RData')
load('R objects/model comparison results.RData')
#-----------------------------------------------------------------------------------
# Model comparison
#-----------------------------------------------------------------------------------
# BIC
N <- nrow(data)
(2*10)*log(N)-(2*-m10$value)
(2*9)*log(N)-(2*-m9a$value)
(2*9)*log(N)-(2*-m9b$value)
(2*9)*log(N)-(2*-m9c$value)
(2*8)*log(N)-(2*-m8$value)
(2*7)*log(N)-(2*-m7$value)
(2*6)*log(N)-(2*-m6$value)
(2*5)*log(N)-(2*-m5$value)
(2*4)*log(N)-(2*-m4$value)
#-----------------------------------------------------------------------------------
# plot best model at various time points
#-----------------------------------------------------------------------------------
p <- m9c$par
N <- 50
xlim <- c(-20,45)
ylim <- c(13,40)

x <- seq(min(xlim), max(xlim),by=0.2)
y <- seq(min(ylim), max(ylim),by=0.2)
cols <- rev(readLines('../tools/colours.txt'))
pal <- colorRampPalette(cols)
dates <- c(14000,round(p[9]))

#pdf(file='../plots/fig 3 rainbelts.pdf',width=14,height=5)

par(mfrow=c(1,2))
for(n in 1:length(dates)){
	date <- dates[n]
	liks <- lik(rep(x,each=length(y)), rep(y,times=length(x)), date, mixS=p[1], aS=0, aN=p[2], sdS=p[3], sdN=p[4], yS=p[5], yN=p[6], rS=p[7], rN=p[8], revD=p[9]) #for model 9c
	lik.mat <- t(matrix(liks,length(y),length(x)))
	image(x,y,sqrt(lik.mat),col=scales::alpha(pal(50),alpha=0), main=paste(date,'BP'),add=FALSE, xlab='Longitude', ylab='Latitude')	
	polygon(x=c(xlim,rev(xlim)), y=rep(ylim,each=2),col='steelblue',border='steelblue')
	maps::map('world', xlim=xlim, ylim=ylim, col='grey90', add=TRUE, fill=TRUE, border='black')
	image(x,y,sqrt(lik.mat),col=scales::alpha(pal(50),alpha=0.8), main=paste(date,'BP'),add=TRUE)	
	}

# fiddly legend
leg.col <- pal(20)
y.scale <- 1.5
Y <- 32
X <- c(38,42)
for (i in 1:(length(leg.col)-1)) {
	y <- Y - (i-1)/y.scale
	rect(X[1],y,X[2],y-1/y.scale, col=leg.col[i], border=NA)
	}

rect(X[1],y-1/y.scale,X[2],Y, col=NA, border='grey')
text(x=X[1], y=Y, 'low density rainfall', cex=0.7, pos=2, offset=0)
text(x=X[1], y=y, 'high density rainfall', cex=0.7, pos=2, offset=0)

#dev.off()
#-----------------------------------------------------------------------------------


