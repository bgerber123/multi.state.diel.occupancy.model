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

#TS vs mean.TS is irrelevant
plot(cov$TS,cov$mean.TS)

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

# number of alpha parameters
K=3 

#jags data input
data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],K=K)

# Parameters monitored
params <- c("alpha", "pNight", "pDay","pND","prob")

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
save(M1.full.no.covs,file="M1.full.no.covs")

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

K=2 
data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],K=K)
params <- c("alpha", "pNight", "pDay","prob")  

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
save(M2.reduced.no.covs,file="AJB Fosa/M2.reduced.no.covs")

#load("AJB Fosa/M2.reduced.no.covs")

#plot(M2.reduced.no.covs,ask=TRUE)

gelman.diag(M2.reduced.no.covs,multivariate = FALSE)

fit <- combine.mcmc(M2.reduced.no.covs)

M2.red.no.covs.CPO=CPO.function(fit,y,"reduced")
CPO.out=t(matrix(c("M2.red.no.covs.CPO",M2.red.no.covs.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Model3 -Null model - No Covariates ############################# 


data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2])
params <- c("alpha", "pdet", "prob") 

model.jags <- jags.model(file="JAGS/jags.multistate.occ.null.R", 
                         data = data.input,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

update(model.jags, n.iter=nb)

M3.null.no.covs <- coda.samples(model.jags, variable.names=params, 
                                n.iter=ni, 
                                thin=nt,
                                progress.bar="text")
save(M3.null.no.covs,file="AJB Fosa/M3.null.no.covs")

#load("AJB Fosa/M3.null.no.covs")

gelman.diag(M3.null.no.covs,multivariate = FALSE)

fit <- combine.mcmc(M3.null.no.covs)

M3.null.no.covs.CPO=CPO.function(fit,y,"null")
CPO.out=t(matrix(c("M3.null.no.covs.CPO",M3.null.no.covs.CPO)))
write.table(CPO.out,file="CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Model1 -Full model - Survey Covariate ############################# 
#Fit the full model with survey/year covariate (effect coding)

#use the survey.cov variable
dim(survey.cov)
head(survey.cov)
dim(survey.cov)

#The X design matrix is the same for all states, but the effects are different
Xday=Xnight=Xnd=survey.cov
K.day=dim(Xday)[2] # number of alpha parameters for day state
K.night=dim(Xnight)[2] # number of alpha parameters for night state
K.nd=dim(Xnd)[2] # number of alpha parameters for night and day state

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],K.day=K.day,
                   K.night=K.night,K.nd=K.nd,Xday=Xday,Xnight=Xnight,Xnd=Xnd)

params <- c("alpha.day","alpha.night","alpha.nd", "pNight", "pDay","pND","prob")

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
save(M1.full.covs.1,file="AJB Fosa/M1.full.covs.1")

#plot(M1.full.covs.1,ask=TRUE)

gelman.diag(M1.full.covs.1,multivariate = FALSE)

fit <- combine.mcmc(M1.full.covs.1)

M1.full.covs.1.CPO=CPO.function(fit,y,"full")
CPO.out=t(matrix(c("M1.full.covs.1.CPO",M1.full.covs.1.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


### Model2 -Reduced model - Survey Covariate ############################# 
#Fit the reduced model with survey/year covariate (effect coding)

#The X design matrix is the same for both states
Xday=Xnight=survey.cov

K.day=ncol(Xday) # number of alpha parameters
K.night=ncol(Xnight) # number of alpha parameters

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],K.day=K.day,
                   K.night=K.night,Xday=Xday,Xnight=Xnight)

params <- c("alpha.day","alpha.night", "pNight", "pDay","prob")

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
save(M2.red.covs.1,file="AJB Fosa/M2.red.covs.1")

#plot(M2.red.covs.1,ask=TRUE)

gelman.diag(M2.red.covs.1,multivariate = FALSE)

fit <- combine.mcmc(M2.red.covs.1)

M2.red.covs.1.CPO=CPO.function(fit,y,"reduced")
CPO.out=t(matrix(c("M2.red.covs.1.CPO",M2.red.covs.1.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


### Model3 -Null model - Survey Covariate ############################# 
X=survey.cov
K=6 # number of alpha parameters

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],K=K,X=X)

# Parameters monitored
params <- c("alpha","pdet", "prob")

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
save(M3.null.covs.1,file="AJB Fosa/M3.null.covs.1")

#plot(M3.null.covs.1,ask=TRUE)

gelman.diag(M3.null.covs.1,multivariate = FALSE)

fit <- combine.mcmc(M3.null.covs.1)

M3.null.covs.1.CPO=CPO.function(fit,y,"null")
CPO.out=t(matrix(c("M3.null.covs.1.CPO",M3.null.covs.1.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)

### Model1 -Full model - Human TS Night/Day Covariate ############################# 
#Since Night TS has no variabilty. Consider Day TS to effect probability of states 2 and 4 only.
Xday=model.matrix(~day.TS,data=cov)
Xnight=model.matrix(~1,data=cov)
Xnd=model.matrix(~day.TS,data=cov)

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],K.day=ncol(Xday),
                  K.night=ncol(Xnight),K.nd=ncol(Xnd),Xday=Xday,Xnight=Xnight,
                  Xnd=Xnd)

params <- c("alpha.day","alpha.night","alpha.nd", "pNight", "pDay","pND","prob")
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
save(M1.full.covs.2,file="AJB Fosa/M1.full.covs.2")

#plot(M1.full.covs.2,ask=TRUE)

gelman.diag(M1.full.covs.2,multivariate = FALSE)

fit <- combine.mcmc(M1.full.covs.2)

M1.full.covs.2.CPO=CPO.function(fit,y,"full")
CPO.out=t(matrix(c("M1.full.covs.2.CPO",M1.full.covs.2.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)


### Model2 -Reduced model - Human TS Night/Day Covariate ############################# 
#Since Night TS has no variabilty. Consider Day TS to effect probability of states 2 and 4 only.
Xday=model.matrix(~day.TS,data=cov)
Xnight=model.matrix(~1,data=cov)

data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],K.day=ncol(Xday),
                   K.night=ncol(Xnight),Xday=Xday,Xnight=Xnight)

params <- c("alpha.day","alpha.night", "pNight", "pDay","prob")
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
save(M2.red.covs.2,file="AJB Fosa/M2.red.covs.2")

#plot(M2.red.covs.2,ask=TRUE)

gelman.diag(M2.red.covs.2,multivariate = FALSE)

fit <- combine.mcmc(M2.red.covs.2)

M2.red.covs.2.CPO=CPO.function(fit,y,"reduced")
CPO.out=t(matrix(c("M2.red.covs.2.CPO",M2.red.covs.2.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)



### Model3 -Null model - Human TS Night/Day Covariate ############################# 
#Since we treat each of the occupied states as the same, we should use mean trap success
#covariate on each state
X=model.matrix(~TS,data=cov)


data.input <- list(y = y, R = dim(y)[1], T = dim(y)[2],K=ncol(X),X=X)

params <- c("alpha","pdet","prob")
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
save(M3.null.covs.2,file="AJB Fosa/M3.null.covs.2")

#plot(M3.null.covs.2,ask=TRUE)

gelman.diag(M3.null.covs.2,multivariate = FALSE)

fit <- combine.mcmc(M3.null.covs.2)

M3.null.covs.2.CPO=CPO.function(fit,y,"null")
CPO.out=t(matrix(c("M3.null.covs.2.CPO",M3.null.covs.2.CPO)))
write.table(CPO.out,file="AJB Fosa/CPO.out.AJB.csv",append=TRUE,col.names = FALSE,sep=",",row.names = FALSE)
