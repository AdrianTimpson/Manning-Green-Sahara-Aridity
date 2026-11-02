#------------------------------------------------------------
# Not for public
# fit gaussian date model to each phase with 14C
# Will not work without access to database!! 
#------------------------------------------------------------
library(BIADconnect)
library(ADMUR)

for(cycle in 1:100){
print(paste('cycle',cycle,'of 100'))

# get required tables from database
conn  <-  init.conn(dbname ='greensahara')
sit <- query.database("SELECT * FROM greensahara.CommonSites")
pha <- query.database("SELECT * FROM greensahara.CommonPhases")
c14 <- query.database("SELECT * FROM greensahara.DatingSamples")
pha <- merge(pha,sit,by='SiteID', all.y=FALSE)
c14 <- subset(c14, !is.na(PhaseID))
c14 <- subset(c14, Method=='C14')
#------------------------------------------------------------
N <- 50
mu.bw <- 100
sigma.bw <- 25
plot <- FALSE
res <- 250

# prioritise phases with no modelled chronology yet
priority <- rep(1,nrow(pha))
priority[is.na(pha$GMM)] <- 2
#--------------------------------------------------------------------------------------

for(n in 1:N){

	# pick a phase at random
	i <- sample(1:nrow(pha),size=1,prob=priority)
	phase <- pha[i,]

	# get other phases with same culture and have a model estimate
	cultures <- phase[,c('Culture')]
	cultures <- cultures[!is.na(cultures)]
	near.phases <- subset(pha, Culture%in%cultures & PhaseID!=phase$PhaseID)
	near.phases <- subset(near.phases, !is.na(GMM))
	if(nrow(near.phases)>0){
		near.phases$dist <- slc(x=phase$Longitude, y=phase$Latitude, ax=near.phases$Longitude, ay=near.phases$Latitude, input='deg') * 6378.1
		}

	# if  no phases available, get phases with same period
	if(nrow(near.phases)==0){
		near.phases <- subset(pha, Period%in%phase$Period & PhaseID!=phase$PhaseID)
		near.phases <- subset(near.phases, !is.na(GMM))
		if(nrow(near.phases)>0){
			near.phases$dist <- slc(x=phase$Longitude, y=phase$Latitude, ax=near.phases$Longitude, ay=near.phases$Latitude, input='deg') * 6378.1
			}
		}

	# weighting by distance from target site using a Gaussian, requires a parameter
	if(nrow(near.phases)!=0){
		weights <- dnorm(near.phases$dist, 0 , 100)
		local.mu <- near.phases$GMM
		local.sigma <- near.phases$GMS
		i <- !is.na(local.mu) & weights!=0
		local.mu <- local.mu[i]
		local.sigma <- local.sigma[i]
		weights <- weights[i]
		weights <- weights/sum(weights)
		NL <- length(weights)
		}
	if(nrow(near.phases)==0)NL <- 0

	# get phase c14 dates
	d <- subset(c14, PhaseID==phase$PhaseID)
	data <- data.frame(age=d$C14Age, sd=d$C14SD)

	# remove outlier C14 dates retaining dates 3x the width of the 50% quantile
	if(nrow(data)>1){
		q <- quantile(data$age,prob=c(0.25,0.5,0.75))
		h <- (q[3]-q[1])*1.5
		r <- q[2] + c(-h,h)
		i <- data$age>r[1] & data$age<r[2]
		data <- data[i,]
		}		

	# Update chronology according to following rules. Note mu and sigma are independent in the prior. 

	if(nrow(data)==0 & NL==0){
		# no useful information. Ensure value is null
		mu <- 'NULL '
		sigma <- 'NULL '
		}

	if(nrow(data)==0 & NL>0){
		# no target 14C, so just use the weighted mean estimates from the local phases
		mu <- round(sum(local.mu * weights))
		sigma <- round(sum(local.sigma * weights))	
		}

	if(nrow(data)>0 & NL==0){
		# no local phases available

		# a  wide domain for mu and sigma	
		mu.range <- estimateDataDomain(data, calcurve=intcal20) + c(-500,500)
		sigma.range <- c(30,500)

		# uniform prior across domain
		prior.matrix <- matrix(1/(res^2),res,res)
		row.names(prior.matrix) <- seq(min(mu.range),max(mu.range),length.out=res)
		colnames(prior.matrix) <- seq(min(sigma.range),max(sigma.range),length.out=res)
		mod.gaussian <- phaseModel(data, calcurve=intcal20, prior.matrix=prior.matrix, model='norm', plot = plot)
		mu <- mod.gaussian$mu
		sigma <- mod.gaussian$sigma
		}

	if(nrow(data)>0 & NL>0){

		# a sensibly wide domain for mu and sigma	
		m1 <- estimateDataDomain(data, calcurve=intcal20) + c(-500,500)
		m2 <- range(local.mu) + c(-500,500)
		mu.range <- c(min(m1[1],m2[1]),max(m1[2],m2[2]))
		sigma.range <- c(30, max(local.sigma)+100)

		# prior joint density using kernels
		d.mu <- density(local.mu,from=mu.range[1],to=mu.range[2],n=res, weights=weights, bw=mu.bw)
		d.sigma <- density(local.sigma,from=sigma.range[1],to=sigma.range[2],n=res, weights=weights, bw=sigma.bw)
		prior.matrix <- matrix(d.mu$y,res,res) * t(matrix(d.sigma$y,res,res))
		prior.matrix <- prior.matrix/sum(prior.matrix)
		row.names(prior.matrix) <- d.mu$x
		colnames(prior.matrix) <- d.sigma$x

		# update prior with data
		mod.gaussian <- phaseModel(data, calcurve=intcal20, prior.matrix=prior.matrix, model='norm', plot = plot)
		mu <- mod.gaussian$mu
		sigma <- mod.gaussian$sigma
		}

	sql.command <- paste("UPDATE `greensahara`.`CommonPhases` SET `GMM`=",mu,", `GMS`=",sigma," WHERE `PhaseID`='",phase$PhaseID,"';",sep='')
	if(!is.nan(mu) & !is.nan(sigma))query.database(conn=conn,sql.command=sql.command)
	}
#--------------------------------------------------------------------------------------
disconnect()
#--------------------------------------------------------------------------------------
}
