#######################################
# Brian D. Gerber 
# 06/29/2020
#
########################################
#We need to compare fitted models (Null, reduced, Full).

#We will do this using CPO.

#This script assumes that simulated data have been fit and saved as model outputs using the 
#other scripts

#States
#1: No use
#2: Day only
#3: Night Only
#4: Night and Day
################################################
rm(list=ls())
setwd("C:/Users/Brian Gerber/Google Drive/MSTOM/multi.state.temporal.activity/data simulation")
setwd("C:/Users/bgerber/Google Drive/MSTOM/multi.state.temporal.activity/data simulation")

#The likelihood funciton is defined in separate script. This is used to calculate CPO
source("multi.state.likelihood.r")


#All three models are fit to each simulated data sets where the underlying true process
#can be different and is either based on the null, reduced, or full models

#This script will compare 3 models fit from a single simulated data set.

#Which simulated data and model fitting iteration to load?
q=1

#We want to output 3 CPO values for each model selection technique
CPO.save=matrix(NA, nrow=1,ncol=3)
colnames(CPO.save)=c("Null","Reduced","Full")

#load model fits and dataset
load(paste("data.null.3.models",q,sep=""))
load(paste("data.reduced.3.models",q,sep=""))
load(paste("data.full.3.models",q,sep=""))
obs.matrix=save.model$obs.matrix #observed data
  
  ##############################################  
  ##############################################  
#Loop over the results for each model fit to get each CPO value

  for (m in 1:3){
    if(m==1){psi.matrix=save.model$models$Null$omega.samples
    det.matrix=save.model$models$Null$det.samples}
    
    if(m==2){psi.matrix=save.model$models$Reduced$omega.samples
    det.matrix=save.model$models$Reduced$det.samples}
    
    if(m==3){psi.matrix=save.model$models$Full$omega.samples
    det.matrix=save.model$models$Full$det.samples}
    
    n.mcmc=dim(det.matrix)[1]
    
    #the likelihood is calcaulted for each site i
    lik.save=matrix(NA, nrow=dim(obs.matrix)[1],ncol=n.mcmc)
    for(k in 1:dim(lik.save)[1]){
      lik.save[k,]=multi.state.likelihood(psi.matrix,det.matrix,obs.matrix[k,])
    }
    #this is the likelihood of each site (row) by n.mcmc (columns)
    #lik.save
    
    #CPO
    CPO.sites=n.mcmc/(apply(1/lik.save,1,sum))
    length(CPO.sites)
    CPO=(-1)*sum(log(CPO.sites))
    CPO.save[q,m]=CPO
}    
    
##################################
#get CPO values. Smaller values indicate a more appropriate model.
CPO.save


