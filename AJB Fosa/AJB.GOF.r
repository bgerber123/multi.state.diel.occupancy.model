#Bayesian p-value for goodness of fit script.

rm(list=ls())
library(runjags)
#library(parallel)
#library(foreach)

source("simulation study/det.matrix.func.r")
source("AJB Fosa/multi.state.likelihood.r")
source("AJB Fosa/GOF.r")

#load most supported model and data
#assign data to objects
load("AJB Fosa/AJB.data")

#assign data to objects
y=as.matrix(AJB.data$data) #detection history

#Drop sites with no detection data
index.site.drop=which(apply(y,1,FUN=function(x){all(is.na(x))}))
y=y[-index.site.drop,]

load("AJB Fosa/M1.full.no.covs.out")
fit <- combine.mcmc(M1.full.no.covs)

model.type="full"

#The GOF function is not very efficient due to a for loop across sites and across mcmc samples.
#This could be sped up significantly.
dev=GOF(fit,y,model.type)
save(dev,file="AJB Fosa/dev.RNP.out")

#plot the Predicted and observed deviances
hist(dev$Deviance.Predicted,breaks=50,col=2)
hist(dev$Deviance.Observed,breaks=50,add=TRUE)
  
#Calculate Bayesian P-vale. Small (<0.1) or large values (>0.9) indicate
#a lack of model fit to the data

length(which(dev$Deviance.Observed>dev$Deviance.Predicted))/length(dev$Deviance.Predicted)
