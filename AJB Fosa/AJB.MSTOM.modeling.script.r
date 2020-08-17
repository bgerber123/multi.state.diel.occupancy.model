### Setup ############################# 
#This script makes use of fosa detection non-detection data from the site AJB
#in northern Madagascar. These data come from Z.J. Farris.
#
#This script will fit the data prepared in the file - AJB.data.script.r

# 4 States:
# 1: No use
# 2: Day only
# 3: Night Only
# 4: Night and Day

#Day was assigned as the period between sunrise and sunset
#Nigt was assigned as the period between sunset and sunrise

rm(list=ls())
logit=function(x){log(x/(1-x))}
expit=function(x){exp(x)/(exp(x)+1)}
library(rjags)
library(runjags)
library(coda)
source("AJB Fosa/multi.state.likelihood.r")
source("AJB Fosa/CPO.function.r")
#Random effects across multiple surveys can be evaluated by CPO with this function:
source("AJB Fosa/CPO.function.RE.r") 


### Data ############################# 
#load the prepared data file
load("AJB Fosa/AJB.data")

#assign data to objects
y=as.matrix(AJB.data$data) #detection history

#Drop sites with no detection data
index.site.drop=which(apply(y,1,FUN=function(x){all(is.na(x))}))
y=y[-index.site.drop,]

#assign covariate data to object
cov=AJB.data$cov
cov=cov[-index.site.drop,]
head(cov)

#Effect Coding Matrix for Survey
survey.cov=AJB.data$conmat

#Look at correlation among human trap success covariates
cor(cov,use="complete.obs")

#compare day and night trap success
plot(cov$day.TS,cov$night.TS)

hist(cov$day.TS,breaks=20)
hist(cov$night.TS,breaks=20)

#There is very little variability in human night TS. Only 5 sites have TS >0
#for all surveys. No reason to include this covariate.

### Set Modeling ############################# 

# Initial values- max state in hierarchy as starting value
zst=rep(4,dim(y)[1])
inits <- function(){list(z = zst)}

# MCMC settings
ni <- 20000  ;       nt <- 2;    nb <- 4000;    nc <- 3;   adapt=4000


### Fit Model1 - Full model - No Covariates ############################# 

#jags data input
data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2])

# Parameters monitored
params <- c("PSI", "pNight", "pDay","pND")

#Fit the model to do adapt phase
model.jags <- jags.model(file="JAGS/jags.multistate.occ.full.R", 
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
save(M1.full.no.covs,file="AJB Fosa/M1.full.no.covs.out")

#load("AJB Fosa/M1.full.no.covs")

#plot(M1.full.no.covs,ask=TRUE)

#Inspect diagnostics
gelman.diag(M1.full.no.covs,multivariate = FALSE)

#combine the chains
fit <- combine.mcmc(M1.full.no.covs)

M1.full.no.covs.CPO=CPO.function(fit,y,"full")
CPO.out=t(matrix(c("M1.full.no.covs",M1.full.no.covs.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Model2 -Reduced model - No Covariates ############################# 

data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2])
params <- c("pNight", "pDay","PSI","psiNight","psiDay")  

model.jags <- jags.model(file="JAGS/jags.multistate.occ.reduced.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)


update(model.jags, n.iter=nb)

M2.reduced.no.covs <- coda.samples(model.jags, variable.names=params, 
                                   n.iter=ni, 
                                   thin=nt,
                                   progress.bar="text")
save(M2.reduced.no.covs,file="AJB Fosa/M2.reduced.no.covs.out")

#load("AJB Fosa/M2.reduced.no.covs")

#plot(M2.reduced.no.covs,ask=TRUE)

gelman.diag(M2.reduced.no.covs,multivariate = FALSE)

fit <- combine.mcmc(M2.reduced.no.covs)

M2.red.no.covs.CPO=CPO.function(fit,y,"reduced")
CPO.out=t(matrix(c("M2.red.no.covs",M2.red.no.covs.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Model3 -Null model - No Covariates ############################# 

data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2])
params <- c("alpha","beta", "p.overall", "PSI") 

model.jags <- jags.model(file="JAGS/jags.multistate.occ.null.alt.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M3.null.no.covs <- coda.samples(model.jags, variable.names=params, 
                                n.iter=ni, 
                                thin=nt,
                                progress.bar="text")
save(M3.null.no.covs,file="AJB Fosa/M3.null.no.covs.out")

#load("AJB Fosa/M3.null.no.covs.out")

gelman.diag(M3.null.no.covs,multivariate = FALSE)

fit <- combine.mcmc(M3.null.no.covs)

M3.null.no.covs.CPO=CPO.function(fit,y,"null")
CPO.out=t(matrix(c("M3.null.no.covs",M3.null.no.covs.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Model1 -Full model - Survey Covariate ############################# 
#Fit the full model with survey/year covariate (effect coding)

#use the survey.cov variable
dim(survey.cov)
head(survey.cov)
dim(survey.cov)

#The X design matrix is the same for all states, but the effects are different
Xday=Xnight=Xnd=survey.cov
Q.day=dim(Xday)[2] # number of alpha parameters for day state
Q.night=dim(Xnight)[2] # number of alpha parameters for night state
Q.nd=dim(Xnd)[2] # number of alpha parameters for night and day state

data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2],Q.day=Q.day,
                   Q.night=Q.night,Q.nd=Q.nd,Xday=Xday,Xnight=Xnight,Xnd=Xnd)

params <- c("alpha.day","alpha.night","alpha.nd", "pNight", "pDay","pND","PSI")

model.jags <- jags.model(file="JAGS/jags.multistate.occ.full.site.covs.by.state.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M1.full.covs.1 <- coda.samples(model.jags, variable.names=params, 
                                n.iter=ni, 
                                thin=nt,
                                progress.bar="text")
save(M1.full.covs.1,file="AJB Fosa/M1.full.covs.1.out")

#plot(M1.full.covs.1,ask=TRUE)

gelman.diag(M1.full.covs.1,multivariate = FALSE)

fit <- combine.mcmc(M1.full.covs.1)

M1.full.covs.1.CPO=CPO.function(fit,y,"full")
CPO.out=t(matrix(c("M1.full.covs.1",M1.full.covs.1.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


### Model2 -Reduced model - Survey Covariate ############################# 
#Fit the reduced model with survey/year covariate (effect coding)

#The X design matrix is the same for both states
Xday=Xnight=survey.cov

Q.day=ncol(Xday) # number of alpha parameters
Q.night=ncol(Xnight) # number of alpha parameters

data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2],Q.day=Q.day,
                   Q.night=Q.night,Xday=Xday,Xnight=Xnight)

params <- c("alpha.day","alpha.night", "pNight", "pDay","PSI","psiDay","psiNight")

model.jags <- jags.model(file="JAGS/jags.multistate.occ.reduced.site.covs.by.state.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M2.red.covs.1 <- coda.samples(model.jags, variable.names=params, 
                               n.iter=ni, 
                               thin=nt,
                               progress.bar="text")
save(M2.red.covs.1,file="AJB Fosa/M2.red.covs.1.out")

#plot(M2.red.covs.1,ask=TRUE)

gelman.diag(M2.red.covs.1,multivariate = FALSE)

fit <- combine.mcmc(M2.red.covs.1)

M2.red.covs.1.CPO=CPO.function(fit,y,"reduced")
CPO.out=t(matrix(c("M2.red.covs.1",M2.red.covs.1.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


### Model3 -Null model - Survey Covariate ############################# 
X=survey.cov
Q=6 # number of alpha parameters

data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2],Q = Q,X=X)

# Parameters monitored
params <- c("alpha","beta","p.overall","psi.overall", "PSI")

model.jags <- jags.model(file="JAGS/jags.multistate.occ.null.site.covs.by.state.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M3.null.covs.1 <- coda.samples(model.jags, variable.names=params, 
                              n.iter=ni, 
                              thin=nt,
                              progress.bar="text")
save(M3.null.covs.1,file="AJB Fosa/M3.null.covs.1.out")

#plot(M3.null.covs.1,ask=TRUE)

gelman.diag(M3.null.covs.1,multivariate = FALSE)

fit <- combine.mcmc(M3.null.covs.1)

M3.null.covs.1.CPO=CPO.function(fit,y,"null")
CPO.out=t(matrix(c("M3.null.covs.1",M3.null.covs.1.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Model1 -Full model - Human TS Night/Day Covariate ############################# 
#Since Night TS has no variabilty. Consider Day TS to effect probability of states 2 and 4 only.
Xday=model.matrix(~day.TS,data=cov)
Xnight=model.matrix(~1,data=cov)
Xnd=model.matrix(~day.TS,data=cov)

data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2],Q.day=ncol(Xday),
                  Q.night=ncol(Xnight),Q.nd=ncol(Xnd),Xday=Xday,Xnight=Xnight,
                  Xnd=Xnd)

params <- c("alpha.day","alpha.night","alpha.nd", "pNight", "pDay","pND","PSI")
model.jags <- jags.model(file="JAGS/jags.multistate.occ.full.site.covs.by.state.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M1.full.covs.2 <- coda.samples(model.jags, variable.names=params, 
                               n.iter=ni, 
                               thin=nt,
                               progress.bar="text")
save(M1.full.covs.2,file="AJB Fosa/M1.full.covs.2.out")

#plot(M1.full.covs.2,ask=TRUE)

gelman.diag(M1.full.covs.2,multivariate = FALSE)

fit <- combine.mcmc(M1.full.covs.2)

M1.full.covs.2.CPO=CPO.function(fit,y,"full")
CPO.out=t(matrix(c("M1.full.covs.2",M1.full.covs.2.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


### Model2 -Reduced model - Human TS Night/Day Covariate ############################# 
#Since Night TS has no variabilty. Consider Day TS to effect probability of states 2 and 4 only.
Xday=model.matrix(~day.TS,data=cov)
Xnight=model.matrix(~1,data=cov)

data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2],Q.day=ncol(Xday),
                   Q.night=ncol(Xnight),Xday=Xday,Xnight=Xnight)

params <- c("alpha.day","alpha.night", "pNight", "pDay","PSI")
model.jags <- jags.model(file="JAGS/jags.multistate.occ.reduced.site.covs.by.state.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M2.red.covs.2 <- coda.samples(model.jags, variable.names=params, 
                               n.iter=ni, 
                               thin=nt,
                               progress.bar="text")
save(M2.red.covs.2,file="AJB Fosa/M2.red.covs.2.out")

#plot(M2.red.covs.2,ask=TRUE)

gelman.diag(M2.red.covs.2,multivariate = FALSE)

fit <- combine.mcmc(M2.red.covs.2)

M2.red.covs.2.CPO=CPO.function(fit,y,"reduced")
CPO.out=t(matrix(c("M2.red.covs.2",M2.red.covs.2.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)



### Model3 -Null model - Human TS Night/Day Covariate ############################# 
#Since we treat each of the occupied states as the same, we should use mean trap success
#covariate on each state
X=model.matrix(~TS,data=cov)


data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2], Q=ncol(X),X=X)

params <- c("alpha","beta","p.overall","psi.overall","psi","PSI")
model.jags <- jags.model(file="JAGS/jags.multistate.occ.null.site.covs.by.state.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M3.null.covs.2 <- coda.samples(model.jags, variable.names=params, 
                              n.iter=ni, 
                              thin=nt,
                              progress.bar="text")
save(M3.null.covs.2,file="AJB Fosa/M3.null.covs.2.out")

#plot(M3.null.covs.2,ask=TRUE)

gelman.diag(M3.null.covs.2,multivariate = FALSE)

fit <- combine.mcmc(M3.null.covs.2)

M3.null.covs.2.CPO=CPO.function(fit,y,"null")
CPO.out=t(matrix(c("M3.null.covs.2",M3.null.covs.2.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


### Model1 -hybrid - no covariates - occupancy is null and detection is full############################# 
data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2])

params <- c("alpha","pND","psi.overall","pNight","pDay","PSI")
model.jags <- jags.model(file="JAGS/jags.multistate.occ.null.det.null.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M1.hybrid.no.covs <- coda.samples(model.jags, variable.names=params, 
                              n.iter=ni, 
                              thin=nt,
                              progress.bar="text")
save(M1.hybrid.no.covs,file="AJB Fosa/M1.hybrid.no.covs.out")

#plot(M1.hybrid.no.covs,ask=TRUE)

gelman.diag(M1.hybrid.no.covs,multivariate = FALSE)

fit <- combine.mcmc(M1.hybrid.no.covs)

M1.hybrid.no.covs.CPO=CPO.function(fit,y,"full") #CPO function will only need "full" for the detection parameters
CPO.out=t(matrix(c("M1.hybrid.no.covs",M1.hybrid.no.covs.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Model2-hybrid - reduced occupancy and detection is full############################# 
data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2])

params <- c("psiDay","pND","psiNight","pNight","pDay","PSI")
model.jags <- jags.model(file="JAGS/jags.multistate.occ.red.det.null.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M2.hybrid.no.covs <- coda.samples(model.jags, variable.names=params, 
                              n.iter=ni, 
                              thin=nt,
                              progress.bar="text")
save(M2.hybrid.no.covs,file="AJB Fosa/M2.hybrid.no.covs.out")

#plot(M2.hybrid.no.covs,ask=TRUE)

gelman.diag(M2.hybrid.no.covs,multivariate = FALSE)

fit <- combine.mcmc(M2.hybrid.no.covs)

M2.hybrid.no.covs.CPO=CPO.function(fit,y,"full") #CPO function will only need "full" for the detection parameters
CPO.out=t(matrix(c("M2.hybrid.no.covs",M2.hybrid.no.covs.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


##########################################
### Fit Model1-Site-Level-Random Effect - Full model - No Covariates ############################# 

#Bring in the array data
y=AJB.data$fosa.array

# Initial values- max state in hierarchy as starting value
zst=matrix(rep(rep(4,dim(y)[1]),dim(y)[3]),ncol=dim(y)[3])
inits <- function(){list(z = zst)}

# MCMC settings
ni <- 40000  ;       nt <- 2;    nb <- 3000;    nc <- 3;   adapt=2000



dim(y)
#jags data input
data.input <- list(y = y, N = dim(y)[1], K = dim(y)[2],Ns=dim(y)[3])

# Parameters monitored
params <- c("alpha1", "alpha2", "alpha3",
            "mu.alpha1","mu.alpha2","mu.alpha3",
            "tau.alpha1","tau.alpha2","tau.alpha3",
            "beta1","beta2","beta3","beta4","beta5",
            "mu.beta1","mu.beta2","mu.beta3","mu.beta4","mu.beta5",
            "tau.beta1","tau.beta2","tau.beta3","tau.beta4","tau.beta5",
            "psiDay","psiNight","psiND","pDay","pNight","pND.ND",
            "pND.N","pND.D","pND.0","PSIM","PSI")

#Fit the model to do adapt phase
model.jags <- jags.model(file="JAGS/jags.multistate.occ.full.alt.RE.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

#burn in period
update(model.jags, n.iter=nb)

#Fit the model and get posterior samples
M1.full.no.covs.RE <- coda.samples(model.jags, variable.names=params, 
                                n.iter=ni, 
                                thin=nt,
                                progress.bar="text")
save(M1.full.no.covs.RE,file="AJB Fosa/M1.full.no.covs.RE.out")
#load("AJB Fosa/M1.full.no.covs.RE.out")

gelman.diag(M1.full.no.covs.RE,multivariate = FALSE)

fit <- combine.mcmc(M1.full.no.covs.RE)

M1.full.no.covs.RE.CPO=CPO.function.RE(fit,y)
CPO.out=t(matrix(c("M1.full.no.covs.RE",M1.full.no.covs.RE.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


