### Setup ############################# 
#This script makes use of fosa detection non-detection data from Makira National Park 
#in northern Madagascar. These data come from Z.J. Farris. Data were collected
#systematically at 7 different sites. 
#
#This script will fit the data prepared in the file - makira.data.script.r

#We will treat sites using a random effect approach, where coefieints come from a higher-level
#distribution. This will account for site-level variation.

# 4 States:
# 1: No use
# 2: Day only
# 3: Night Only
# 4: Night and Day

#Day was assigned as the period between sunrise and sunset
#Nigt was assigned as the period between sunset and sunrise

rm(list=ls())
setwd("C:/Users/Brian Gerber/Google Drive/MSTOM Project/multi.state.temporal.activity")
setwd("C:/Users/bgerber/Google Drive/MSTOM Project/multi.state.temporal.activity")
logit=function(x){log(x/(1-x))}
expit=function(x){exp(x)/(exp(x)+1)}
library(rjags)
library(runjags)
library(coda)
source("Makira Fosa2/multi.state.likelihood.r")
source("Makira Fosa2/CPO.function.RE.r") 


### Data ############################# 
#load the prepared data file
load("Makira Fosa2/Makira.data2")

#assign data to objects
y=Makira.data2[[2]] #detection history, sites x occs x survey area

#assign covariate data to object
cov=Makira.data2[[3]]

#How many camera sites are sampled for each survey area?
N=Makira.data2[[1]]

#How many sites surveyd?
T=length(N)

#How many occasions sampled (this is the same for all survey areas)
K=dim(y)[2]

### Set Modeling ############################# 

# Initial values- max state in hierarchy as starting value
zst=matrix(rep(rep(4,dim(y)[1]),dim(y)[3]),ncol=dim(y)[3])
#Make sure that camera sites that are not surveyed (and not sampled via MCMC) have 
#no initial values.

zst[(N[2]+1):N[1],2]=NA
zst[(N[3]+1):N[1],3]=NA
zst[(N[4]+1):N[1],4]=NA
zst[(N[5]+1):N[1],5]=NA
zst[(N[6]+1):N[1],6]=NA
zst[(N[7]+1):N[1],7]=NA

inits <- function(){list(z = zst)}

# MCMC settings
ni <- 20000  ;       nt <- 2;    nb <- 4000;    nc <- 3;   adapt=4000


### Fit Model1 - Full model - No Covariates ############################# 

#jags data input
data.input <- list(y = y, N = N, K = K, T = T)

# Parameters monitored
params <- c("alpha1", "alpha2", "alpha3",
            "mu.alpha1","mu.alpha2","mu.alpha3",
            "tau.alpha1","tau.alpha2","tau.alpha3",
            "beta1","beta2","beta3","beta4","beta5",
            "mu.beta1","mu.beta2","mu.beta3","mu.beta4","mu.beta5",
            "tau.beta1","tau.beta2","tau.beta3","tau.beta4","tau.beta5",
            "psiDay","psiNight","psiND","pDay","pNight","pND.ND",
            "pND.N","pND.D","pND.0","PSI","PSIM")

#Fit the model to do adapt phase
model.jags <- jags.model(file="JAGS/jags.multistate.occ.full.alt.RE.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M1.full.no.covs <- coda.samples(model.jags, variable.names=params, 
                                n.iter=ni, 
                                thin=nt,
                                progress.bar="text")
save(M1.full.no.covs,file="Makira Fosa2/M1.full.no.covs.out")

#load("Makira Fosa2/M1.full.no.covs.out")

#plot(M1.full.no.covs,ask=TRUE)

#Inspect diagnostics
gelman.diag(M1.full.no.covs,multivariate = FALSE)

#combine the chains
fit <- combine.mcmc(M1.full.no.covs)

M1.full.no.covs.CPO=CPO.function.RE(fit,y)
CPO.out=t(matrix(c("M1.full.no.covs",M1.full.no.covs.CPO)))
write.table(CPO.out,file="Makira Fosa2/CPO.out.Makira.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


### Model1 -Reduced model - No Covariates ############################# 
#jags data input
data.input <- list(y = y, N = N, K = K, T = T)

params <- c("alpha1", "alpha2", 
            "mu.alpha1","mu.alpha2",
            "tau.alpha1","tau.alpha2",
            "beta1","beta2",
            "mu.beta1","mu.beta2",
            "tau.beta1","tau.beta2",
            "psiDay","psiNight","psiND","pDay","pNight","pND.ND",
            "pND.N","pND.D","pND.0","PSI","PSIM")
#Fit the model to do adapt phase
model.jags <- jags.model(file="JAGS/jags.multistate.occ.reduced.alt.RE.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M1.red.no.covs <- coda.samples(model.jags, variable.names=params, 
                                n.iter=ni, 
                                thin=nt,
                                progress.bar="text")
save(M1.red.no.covs,file="Makira Fosa2/M1.red.no.covs.out")

#load("Makira Fosa2/M1.red.no.covs.out")

#plot(M1.red.no.covs,ask=TRUE)

#Inspect diagnostics
gelman.diag(M1.red.no.covs,multivariate = FALSE)

#combine the chains
fit <- combine.mcmc(M1.red.no.covs)

M1.red.no.covs.CPO=CPO.function.RE(fit,y)
CPO.out=t(matrix(c("M1.red.no.covs",M1.red.no.covs.CPO)))
write.table(CPO.out,file="Makira Fosa2/CPO.out.Makira.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Model1 -Null model - No Covariates ############################# 

#jags data input
data.input <- list(y = y, N = N, K = K, T = T)

params <- c("alpha", "mu.alpha", "tau.alpha",
            "beta","mu.beta","tau.beta",
            "psiDay","psiNight","psiND","pDay","pNight","pND.ND",
            "pND.N","pND.D","pND.0","PSI","PSIM")
#Fit the model to do adapt phase
model.jags <- jags.model(file="JAGS/jags.multistate.occ.null.alt.RE.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M1.null.no.covs <- coda.samples(model.jags, variable.names=params, 
                               n.iter=ni, 
                               thin=nt,
                               progress.bar="text")
save(M1.null.no.covs,file="Makira Fosa2/M1.null.no.covs.out")

#load("Makira Fosa2/M1.null.no.covs.out")

#plot(M1.null.no.covs,ask=TRUE)

#Inspect diagnostics
gelman.diag(M1.null.no.covs,multivariate = FALSE)

#combine the chains
fit <- combine.mcmc(M1.null.no.covs)

M1.null.no.covs.CPO=CPO.function.RE(fit,y)
CPO.out=t(matrix(c("M1.null.no.covs",M1.null.no.covs.CPO)))
write.table(CPO.out,file="Makira Fosa2/CPO.out.Makira.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Fit Model2 - Full model - Human TS ############################# 

#jags data input
data.input <- list(y = y, N = N, K = K, T = T)

#Need to define cov1, cov2, and cov3, that respectively influence day, night, and ND occurence
head(cov)

cov1=cov[,4,]
cov2=cov[,5,]
cov3=cov[,4,]+cov[,5,]

#jags data input
data.input <- list(y = y, K = dim(y)[2], N=N,T=T,
                   cov1=cov1,cov2=cov2,cov3=cov3)


# Parameters monitored
params <- c("alpha1", "alpha2", "alpha3",
            "alpha1.1","alpha2.1","alpha3.1",
            "mu.alpha.day","mu.alpha.night","mu.alpha.ND",
            "tau.alpha.day","tau.alpha.night","tau.alpha.ND",
            "tau.beta.day","tau.beta.night","tau.beta.ND.D","tau.beta.ND.N","tau.beta.ND.ND",
            "mu.beta.day","mu.beta.night","mu.beta.ND.D","mu.beta.ND.N","mu.beta.ND.ND",
            "psiDay","psiNight","psiND","pDay","pNight","pND.ND",
            "pND.N","pND.D","pND.0","PSI","PSIM")

#Fit the model to do adapt phase
model.jags <- jags.model(file="JAGS/jags.multistate.occ.full.site.covs.RE.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M2.full.covs <- coda.samples(model.jags, variable.names=params, 
                                n.iter=ni, 
                                thin=nt,
                                progress.bar="text")
save(M2.full.covs,file="Makira Fosa2/M2.full.covs.out")

#load("Makira Fosa2/M2.full.covs.out")

#plot(M2.full.covs,ask=TRUE)

#Inspect diagnostics
gelman.diag(M2.full.covs,multivariate = FALSE)

#combine the chains
fit <- combine.mcmc(M2.full.covs)

M2.full.covs.CPO=CPO.function.RE(fit,y)
CPO.out=t(matrix(c("M2.full.covs",M2.full.covs.CPO)))
write.table(CPO.out,file="Makira Fosa2/CPO.out.Makira.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Fit Model2 - Reduced model - Human TS ############################# 


#Need to define cov1, cov2, and cov3, that respectively influence day, night, and ND occurence
head(cov)

cov1=cov[,4,]
cov2=cov[,5,]

#jags data input
data.input <- list(y = y, K = dim(y)[2], N=N,T=T,
                   cov1=cov1,cov2=cov2)

# Parameters monitored
params <- c("alpha1", "alpha2", "alpha1.1","alpha2.1",
            "mu.alpha.day","mu.alpha.night",
            "tau.alpha.day","tau.alpha.night",
            "tau.beta.day","tau.beta.night",
            "mu.beta.day","mu.beta.night",
            "psiDay","psiNight","psiND","pDay","pNight","pND.ND",
            "pND.N","pND.D","pND.0","PSI","PSIM")

#Fit the model to do adapt phase
model.jags <- jags.model(file="JAGS/jags.multistate.occ.reduced.site.covs.RE.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M2.red.covs <- coda.samples(model.jags, variable.names=params, 
                             n.iter=ni, 
                             thin=nt,
                             progress.bar="text")
save(M2.red.covs,file="Makira Fosa2/M2.red.covs.out")

#load("Makira Fosa2/M2.red.covs.out")

#plot(M2.red.covs,ask=TRUE)

#Inspect diagnostics
gelman.diag(M2.red.covs,multivariate = FALSE)

#combine the chains
fit <- combine.mcmc(M2.red.covs)

M2.red.covs.CPO=CPO.function.RE(fit,y)
CPO.out=t(matrix(c("M2.red.covs",M2.red.covs.CPO)))
write.table(CPO.out,file="Makira Fosa2/CPO.out.Makira.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


### Fit Model2 - Null model - Human TS ############################# 

#Need to define cov1, cov2, and cov3, that respectively influence day, night, and ND occurence
head(cov)
cov=cov[,3,]

#jags data input
data.input <- list(y = y, K = dim(y)[2], N=N,T=T,cov=cov)

# Parameters monitored
params <- c("alpha1", "alpha1.1",
            "mu.alpha",
            "tau.alpha",
            "tau.beta",
            "mu.beta",
            "psiDay","psiNight","psiND","pDay","pNight","pND.ND",
            "pND.N","pND.D","pND.0","PSI","PSIM")

#Fit the model to do adapt phase
model.jags <- jags.model(file="JAGS/jags.multistate.occ.null.site.covs.RE.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M2.null.covs <- coda.samples(model.jags, variable.names=params, 
                            n.iter=ni, 
                            thin=nt,
                            progress.bar="text")
save(the,file="Makira Fosa2/M2.null.covs.out")

#load("Makira Fosa2/M2.null.covs.out")

#plot(M2.null.covs,ask=TRUE)

#Inspect diagnostics
gelman.diag(M2.null.covs,multivariate = FALSE)

#combine the chains
fit <- combine.mcmc(M2.null.covs)

M2.null.covs.CPO=CPO.function.RE(fit,y)
CPO.out=t(matrix(c("M2.null.covs",M2.null.covs.CPO)))
write.table(CPO.out,file="Makira Fosa2/CPO.out.Makira.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


