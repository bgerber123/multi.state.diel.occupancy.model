#Bayesian p-value for goodness of fit script.

rm(list=ls())
library(runjags)
#library(parallel)
#library(foreach)

source("./Simulation Files/det.matrix.func.r")
source("./RNP Fosa/multi.state.likelihood.r")
source("./RNP Fosa/GOF.r")

#load most supported model and data
load("./RNP Fosa/RNP.data")

#assign data to objects
y=RNP.data[[1]] #detection history

load("./RNP Fosa/M1.full.out")
fit=combine.mcmc(M1.full)


model.type="full"

#The GOF function takes some time. It loops over each site and does foreach for all mcmc samples
#This could be sped up.
dev=GOF(fit,y,model.type)
save(dev,file="./RNP Fosa/dev.RNP.out")

#plot the Predicted and observed deviances
hist(dev$Deviance.Predicted,breaks=50,col=2)
hist(dev$Deviance.Observed,breaks=50,add=TRUE)
  
#Calculate Bayesian P-vale. Small (<0.1) or large values (>0.9) indicate
#a lack of model fit to the data

length(which(dev$Deviance.Observed>dev$Deviance.Predicted))/length(dev$Deviance.Predicted)
