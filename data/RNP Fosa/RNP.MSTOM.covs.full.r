###################################################
#This script fits the Ranomafana National Park fosa occupancy camera trap data - Data from Gerber et al. 2012
#Gerber, B. D., Karpanty, S. M., & Randrianantenaina, J. (2012). The impact of forest logging and fragmentation on carnivore species composition, density and occupancy in Madagascar's rainforests. Oryx, 46(3), 414-422.
#
#Specifically, it uses the camera surveys at Sahamalaotra and Valohoaka-Vatoharanana within the park.
#The states are assigned as,

# 4 States:
# 1: No use
# 2: Day only
# 3: Night Only
# 4: Night and Day


#Day was assigned as the period between sunrise and sunset
#Nigt was assigned as the period between sunset and sunrise

rm(list=ls())
setwd("C:/Users/Brian Gerber/Google Drive/MSTOM/multi.state.temporal.activity/data/RNP Fosa")
logit=function(x){log(x/(1-x))}
expit=function(x){exp(x)/(exp(x)+1)}
library(rjags)
library(runjags)
library(coda)
source("multi.state.likelihood.r")

#load the prepared data file
load("RNP.data")

#assign data to objects
y=RNP.data[[1]] #detection history
covs=RNP.data[[2]]

#FYI-The Locals covariate- all detections of people are during the day. There are no nighttime detections.
#Locals covariate is probably not a good measure of overall human/forest disturbance.

#Look at covariate correlations
cor(covs$DistMatrix,covs$Locals)
cor(covs$DistMatrix,covs$DistTown)
cor(covs$DistTown,covs$Locals)

hist(covs$DistMatrix)
hist(covs$Locals)
hist(covs$DistTown)

#Hypothesis:
#Each covariate, representing an effect of edge/people will influence fosa occurence where
#the fosa will have lower occurence near edge/people during the day and this will increase
#the further a site is from edge/people.

#First consider the FUll MSTOM with each covariate separately. Assume that the major source in
#detection is based on the state (2,3,4).


#################################################################
#################################################################
#MODEL 1- DistTown covariate - FULL model
cov=covs$DistTown
hist(cov)
cov1.unscaled=as.numeric(cov)
cov1.scaled=as.numeric(scale(cov,center = TRUE,scale = TRUE))

K=6 # number of alpha parameters

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],
                   cov1=cov1.scaled,K=K)

# Initial values- max state in hierachy as starting value
zst=rep(4,dim(y)[1])
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha", "pNight", "pDay","pND","prob") #"p"

# MCMC settings
ni <- 20000  ;       nt <- 2;    nb <- 4000;    nc <- 3;   adapt=4000

#Fit the model to do adapt phase
#this is the full model
model.jags <- jags.model(file="jags.MSTOM.full.covs", 
                                  data = data.input,
                                  inits=inits,
                                  n.chains = nc,
                                  n.adapt=adapt)


#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M1.full <- coda.samples(model.jags, variable.names=params, 
                                      n.iter=ni, 
                                      thin=nt,
                                      progress.bar="text")
save(M1.full,file="M1.full")

#############################################################
#############################################################
#plot(M1.full,ask=TRUE)

#Inspect diagnostics
gelman.diag(M1.full,multivariate = FALSE)

# Summarize posteriors
print(M1.full, dig = 2)


#combine the chains
fit <- combine.mcmc(M1.full)

str(fit)
attributes(fit)

dim(fit)
colnames(fit)

#To do model selection with CPO, we need to get site-level occurence and detection
#probabilities.

pDay=fit[,which(grepl("pDay",colnames(fit)))]
pNight=fit[,which(grepl("pNight",colnames(fit)))]

pND.0=fit[,which(grepl("pND",colnames(fit)))[1]]
pND.D=fit[,which(grepl("pND",colnames(fit)))[2]]
pND.N=fit[,which(grepl("pND",colnames(fit)))[3]]
pND.ND=fit[,which(grepl("pND",colnames(fit)))[4]]

det.matrix=cbind(pDay,
                 pNight,
                 pND.0,
                 pND.D,
                 pND.N,
                 pND.ND)

#the likelihood is calcaulted for each site k
n.mcmc=length(pDay)

lik.save=matrix(NA, nrow=dim(y)[1],ncol=n.mcmc)

#loop through each site k
for(k in 1:dim(lik.save)[1]){
  
  site.num=k
  index.occ=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("prob[",site.num,",",sep=""), x, fixed = TRUE)}))
  
  #For site q in order of states by 1,2,3,4

  psi.matrix=fit[,index.occ]

  #calcualte site level likelihood
  lik.save[k,]=multi.state.likelihood(psi.matrix,det.matrix,y[k,])
}


#CPO calculation from site-level liklihood
CPO.sites=n.mcmc/(apply(1/lik.save,1,sum))
length(CPO.sites)
M1.full.CPO=(-1)*sum(log(CPO.sites))
M1.full.CPO

############################################################################
############################################################################
############################################################################
#MODEL 2- DistMatrix covariate - FULL model
cov=covs$DistMatrix
hist(cov)
cov1.unscaled=as.numeric(cov)
cov1.scaled=as.numeric(scale(cov,center = TRUE,scale = TRUE))

K=6 # number of alpha parameters

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],
                   cov1=cov1.scaled,K=K)

# Initial values- max state in hierachy as starting value
zst=rep(4,dim(y)[1])
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha", "pNight", "pDay","pND","prob")  #include p?

# MCMC settings
ni <- 20000  ;       nt <- 2;    nb <- 4000;    nc <- 3;   adapt=4000

#Fit the model to do adapt phase
#this is the full model
model.jags <- jags.model(file="jags.MSTOM.full.covs", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)


#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M2.full <- coda.samples(model.jags, variable.names=params, 
                        n.iter=ni, 
                        thin=nt,
                        progress.bar="text")
save(M2.full,file="M2.full")

#############################################################
#############################################################
#plot(M2.full,ask=TRUE)

#Inspect Diagnostics
gelman.diag(M2.full,multivariate = FALSE)

# Summarize posteriors
print(M2.full, dig = 2)


#combine the chains
fit <- combine.mcmc(M2.full)

#To do model selection with CPO, we need to get site-level occurence and detection
#probabilities.

pDay=fit[,which(grepl("pDay",colnames(fit)))]
pNight=fit[,which(grepl("pNight",colnames(fit)))]

pND.0=fit[,which(grepl("pND",colnames(fit)))[1]]
pND.D=fit[,which(grepl("pND",colnames(fit)))[2]]
pND.N=fit[,which(grepl("pND",colnames(fit)))[3]]
pND.ND=fit[,which(grepl("pND",colnames(fit)))[4]]

det.matrix=cbind(pDay,
                 pNight,
                 pND.0,
                 pND.D,
                 pND.N,
                 pND.ND)

n.mcmc=length(pDay)

#the likelihood is calcaulted for each site k
lik.save=matrix(NA, nrow=dim(y)[1],ncol=n.mcmc)

#loop through each site k
for(k in 1:dim(lik.save)[1]){
  
  site.num=k
  index.occ=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("prob[",site.num,",",sep=""), x, fixed = TRUE)}))
  
  psi.matrix=fit[,index.occ]

  #calcualte site levle likelihood
  lik.save[k,]=multi.state.likelihood(psi.matrix,det.matrix,y[k,])
}


#CPO calculation from site-level liklihood
CPO.sites=n.mcmc/(apply(1/lik.save,1,sum))
length(CPO.sites)
M2.full.CPO=(-1)*sum(log(CPO.sites))
M2.full.CPO


###################################################################
###################################################################
###################################################################
#Consider the same covariate models as above, but with the Reduced model

#MODEL 3- DistTown covariate - Reduced model
cov=covs$DistTown
hist(cov)
cov1.unscaled=as.numeric(cov)
cov1.scaled=as.numeric(scale(cov,center = TRUE,scale = TRUE))

K=4 # number of alpha parameters

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],
                   cov1=cov1.scaled,K=K)

# Initial values- max state in hierachy as starting value
zst=rep(4,dim(y)[1])
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha", "pNight", "pDay","prob") #"p"

# MCMC settings
ni <- 10000  ;       nt <- 2;    nb <- 2000;    nc <- 3;   adapt=1000

#Fit the model to do adapt phase
#this is the full model
model.jags <- jags.model(file="jags.MSTOM.reduced.covs", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)


#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M3.red <- coda.samples(model.jags, variable.names=params, 
                        n.iter=ni, 
                        thin=nt,
                        progress.bar="text")
save(M3.red,file="M3.red")

#############################################################
#plot(M3.red,ask=TRUE)

#Inspect Diagnostics
gelman.diag(M3.red,multivariate = FALSE)

# Summarize posteriors
print(M3.red, dig = 2)


#combine the chains
fit <- combine.mcmc(M3.red)

#To do model selection with CPO, we need to get site-level occurence and detection
#probabilities.

pDay=fit[,which(grepl("pDay",colnames(fit)))]
pNight=fit[,which(grepl("pNight",colnames(fit)))]

n.mcmc=length(pDay)

pND.0 <- (1-pNight)*(1-pDay)
pND.D <- pDay*(1-pNight)
pND.N <- pNight*(1-pDay)
pND.ND <- pDay*pNight 


det.matrix=cbind(pDay,
                 pNight,
                 pND.0,
                 pND.D,
                 pND.N,
                 pND.ND)

n.mcmc=length(pDay)

#the likelihood is calcaulted for each site k
lik.save=matrix(NA, nrow=dim(y)[1],ncol=n.mcmc)

#loop through each site k
for(k in 1:dim(lik.save)[1]){
  
  site.num=k
  index.occ=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("prob[",site.num,",",sep=""), x, fixed = TRUE)}))
  
  psi.matrix=fit[,index.occ]

  #calcualte site level likelihood
  lik.save[k,]=multi.state.likelihood(psi.matrix,det.matrix,y[k,])
}


#CPO calculation from site-level liklihood
CPO.sites=n.mcmc/(apply(1/lik.save,1,sum))
length(CPO.sites)
M3.red.CPO=(-1)*sum(log(CPO.sites))
M3.red.CPO

#################################################
#################################################
#MODEL 4- DistMatrix covariate - Reduced model
cov=covs$DistMatrix
hist(cov)
cov1.unscaled=as.numeric(cov)
cov1.scaled=as.numeric(scale(cov,center = TRUE,scale = TRUE))

K=4 # number of alpha parameters

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],
                   cov1=cov1.scaled,K=K)

# Initial values- max state in hierachy as starting value
zst=rep(4,dim(y)[1])
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha", "pNight", "pDay","prob") #"p"

# MCMC settings
ni <- 10000  ;       nt <- 2;    nb <- 2000;    nc <- 3;   adapt=1000

#Fit the model to do adapt phase
#this is the full model
model.jags <- jags.model(file="jags.MSTOM.reduced.covs", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)


#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M4.red <- coda.samples(model.jags, variable.names=params, 
                       n.iter=ni, 
                       thin=nt,
                       progress.bar="text")
save(M4.red,file="M4.red")

#############################################################
#############################################################
#plot(M4.red,ask=TRUE)

#Inspect Diagnostics
gelman.diag(M4.red,multivariate = FALSE)

# Summarize posteriors
print(M4.red, dig = 2)


#combine the chains
fit <- combine.mcmc(M4.red)

#To do model selection with CPO, we need to get site-level occurence and detection
#probabilities.

pDay=fit[,which(grepl("pDay",colnames(fit)))]
pNight=fit[,which(grepl("pNight",colnames(fit)))]

n.mcmc=length(pDay)

pND.0 <- (1-pNight)*(1-pDay)
pND.D <- pDay*(1-pNight)
pND.N <- pNight*(1-pDay)
pND.ND <- pDay*pNight 


det.matrix=cbind(pDay,
                 pNight,
                 pND.0,
                 pND.D,
                 pND.N,
                 pND.ND)


#the likelihood is calcaulted for each site k
lik.save=matrix(NA, nrow=dim(y)[1],ncol=n.mcmc)

#loop through each site k
for(k in 1:dim(lik.save)[1]){
  
  site.num=k
  index.occ=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("prob[",site.num,",",sep=""), x, fixed = TRUE)}))
  
  psi.matrix=fit[,index.occ]
  
  #calcualte site level likelihood
  lik.save[k,]=multi.state.likelihood(psi.matrix,det.matrix,y[k,])
}


#CPO calculation from site-level liklihood
CPO.sites=n.mcmc/(apply(1/lik.save,1,sum))
length(CPO.sites)
M4.red.CPO=(-1)*sum(log(CPO.sites))
M4.red.CPO

######################################################
######################################################
#Need to fit the Null model with covariates

#MODEL 5- Distown covariate - Null model
cov=covs$DistTown
hist(cov)
cov1.unscaled=as.numeric(cov)
cov1.scaled=as.numeric(scale(cov,center = TRUE,scale = TRUE))

K=2 # number of alpha parameters

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],
                   cov1=cov1.scaled,K=K)

# Initial values- max state in hierachy as starting value
zst=rep(4,dim(y)[1])
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha", "pdet", "prob") #"p"

# MCMC settings
ni <- 10000  ;       nt <- 2;    nb <- 2000;    nc <- 3;   adapt=1000

#Fit the model to do adapt phase
#this is the full model
model.jags <- jags.model(file="jags.MSTOM.null.covs", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)


#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M5.null <- coda.samples(model.jags, variable.names=params, 
                       n.iter=ni, 
                       thin=nt,
                       progress.bar="text")
save(M5.null,file="M5.null")

#############################################################
#############################################################
#plot(M5.null)

#Inspect Diagnostics
gelman.diag(M5.null,multivariate = FALSE)

# Summarize posteriors
print(M5.null, dig = 2)


#combine the chains
fit <- combine.mcmc(M5.null)

str(fit)
attributes(fit)

dim(fit)
colnames(fit)

#To do model selection with CPO, we need to get site-level occurence and detection
#probabilities.

pDay=fit[,which(grepl("pdet",colnames(fit)))]
pNight=pDay

n.mcmc=length(pDay)

pND.0 <- (1-pNight)*(1-pDay)
pND.D <- pDay*(1-pNight)
pND.N <- pNight*(1-pDay)
pND.ND <- pDay*pNight 


det.matrix=cbind(pDay,
                 pNight,
                 pND.0,
                 pND.D,
                 pND.N,
                 pND.ND)

n.mcmc=length(pDay)

#the likelihood is calcaulted for each site k
lik.save=matrix(NA, nrow=dim(y)[1],ncol=n.mcmc)

#loop through each site k
for(k in 1:dim(lik.save)[1]){
  
  site.num=k
  index.occ=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("prob[",site.num,",",sep=""), x, fixed = TRUE)}))
  
  psi.matrix=fit[,index.occ]
  
  #calcualte site level likelihood
  lik.save[k,]=multi.state.likelihood(psi.matrix,det.matrix,y[k,])
}


#CPO calculation from site-level liklihood
CPO.sites=n.mcmc/(apply(1/lik.save,1,sum))
length(CPO.sites)
M5.null.CPO=(-1)*sum(log(CPO.sites))
M5.null.CPO

#####################################
#MODEL 6- DistMatrix covariate - Null model
cov=covs$DistMatrix
hist(cov)
cov1.unscaled=as.numeric(cov)
cov1.scaled=as.numeric(scale(cov,center = TRUE,scale = TRUE))

K=2 # number of alpha parameters

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],
                   cov1=cov1.scaled,K=K)

# Initial values- max state in hierachy as starting value
zst=rep(4,dim(y)[1])
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha", "pdet", "prob") #"p"

# MCMC settings
ni <- 10000  ;       nt <- 2;    nb <- 2000;    nc <- 3;   adapt=1000

#Fit the model to do adapt phase
#this is the full model
model.jags <- jags.model(file="jags.MSTOM.null.covs", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)


#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M6.null <- coda.samples(model.jags, variable.names=params, 
                        n.iter=ni, 
                        thin=nt,
                        progress.bar="text")
save(M6.null,file="M6.null")

#############################################################
#############################################################
#plot(M6.null)
#gelman.diag(M6.null,multivariate = FALSE)

# Summarize posteriors
print(M6.null, dig = 2)


#combine the chains
fit <- combine.mcmc(M6.null)

str(fit)
attributes(fit)

dim(fit)
colnames(fit)

#To do model selection with CPO, we need to get site-level occurence and detection
#probabilities.

pDay=fit[,which(grepl("pdet",colnames(fit)))]
pNight=pDay

n.mcmc=length(pDay)

pND.0 <- (1-pNight)*(1-pDay)
pND.D <- pDay*(1-pNight)
pND.N <- pNight*(1-pDay)
pND.ND <- pDay*pNight 


det.matrix=cbind(pDay,
                 pNight,
                 pND.0,
                 pND.D,
                 pND.N,
                 pND.ND)

n.mcmc=length(pDay)

#the likelihood is calcaulted for each site k
lik.save=matrix(NA, nrow=dim(y)[1],ncol=n.mcmc)

#loop through each site k
for(k in 1:dim(lik.save)[1]){
  
  site.num=k
  index.occ=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("prob[",site.num,",",sep=""), x, fixed = TRUE)}))
  
  psi.matrix=fit[,index.occ]
  
  #calcualte site level likelihood
  lik.save[k,]=multi.state.likelihood(psi.matrix,det.matrix,y[k,])
}


#CPO calculation from site-level liklihood
CPO.sites=n.mcmc/(apply(1/lik.save,1,sum))
length(CPO.sites)
M6.null.CPO=(-1)*sum(log(CPO.sites))
M6.null.CPO

##########################################
##########################################
#Fit the null model with no covariates

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2])

# Initial values- max state in hierachy as starting value
zst=rep(4,dim(y)[1])
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha", "pdet", "prob") 

# MCMC settings
ni <- 10000  ;       nt <- 2;    nb <- 2000;    nc <- 3;   adapt=1000

#Fit the model to do adapt phase
#this is the full model
model.jags <- jags.model(file="jags.MSTOM.null", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)


#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M7.null <- coda.samples(model.jags, variable.names=params, 
                        n.iter=ni, 
                        thin=nt,
                        progress.bar="text")
save(M7.null,file="M7.null")

#############################################################
#############################################################
#plot(M7.null)

#Inspect R-hat diagnostics
gelman.diag(M7.null,multivariate = FALSE)

#combine the chains
fit <- combine.mcmc(M7.null)

# Summarize posteriors
print(fit, dig = 2)


str(fit)
attributes(fit)

dim(fit)
colnames(fit)

#To do model selection with CPO, we need to get site-level occurence and detection
#probabilities.

pDay=fit[,which(grepl("pdet",colnames(fit)))]
pNight=pDay

n.mcmc=length(pDay)

pND.0 <- (1-pNight)*(1-pDay)
pND.D <- pDay*(1-pNight)
pND.N <- pNight*(1-pDay)
pND.ND <- pDay*pNight 


det.matrix=cbind(pDay,
                 pNight,
                 pND.0,
                 pND.D,
                 pND.N,
                 pND.ND)

n.mcmc=length(pDay)

#the likelihood is calcaulted for each site k
lik.save=matrix(NA, nrow=dim(y)[1],ncol=n.mcmc)

#loop through each site k
for(k in 1:dim(lik.save)[1]){
  
  site.num=k
  index.occ=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("prob[",site.num,",",sep=""), x, fixed = TRUE)}))
  
  psi.matrix=fit[,index.occ]
  
  #calcualte site level likelihood
  lik.save[k,]=multi.state.likelihood(psi.matrix,det.matrix,y[k,])
}


#CPO calculation from site-level liklihood
CPO.sites=n.mcmc/(apply(1/lik.save,1,sum))
length(CPO.sites)
M7.null.CPO=(-1)*sum(log(CPO.sites))
M7.null.CPO


########################################
########################################
#Model Comparison
CPOs=c(M1.full.CPO,M2.full.CPO,M3.red.CPO,M4.red.CPO,M5.null.CPO,M6.null.CPO,M7.null.CPO)
names(CPOs)=c("M1","M2","M3","M4","M5","M6","M7")
CPOs[order(CPOs)]
