#######################################
# Brian D. Gerber 
# 06/29/2020
#
########################################

#Assuming data has been simulated and an object has been output that stored the simulated data
# See scripts sim.data.MSTOM.full.model.r
#             sim.data.MSTOM.reduced.model.r)
#             sim.data.MSTOM.null.model.r).

#States
#1: No use
#2: Day only
#3: Night Only
#4: Night and Day


#Below, we will fit three multi-state temporal occurence models, in the following order to the loaded 
#simulated data: 
#1) Null
#2) Reduced
#3) Full

################################################
#Setup workspace
  rm(list=ls())
  setwd("C:/Users/Brian Gerber/Google Drive/MSTOM/multi.state.temporal.activity/data simulation")
  setwd("C:/Users/bgerber/Google Drive/MSTOM/multi.state.temporal.activity/data simulation")
  
  #Requires the program JAGS to be installed
  library(jagsUI)
  library(rjags)

#Load the simulated data - choose one
  load("sim.data.multistate.null") 
#  load("sim.data.multistate.reduced") 
#  load("sim.data.multistate.full") 
  
  
  # MCMC settings
  ni <- 10000  ;       nt <- 1;        nb <- 1000;  nc <- 1;  adapt <- 1000


#For each qth iteration, fit the three different models
  n.sim=sim.data$n.sim

########################################
########################################
#Start for loop
for(q in 1:n.sim){

#the first saved list object will be the data
  save.model=vector("list",2)

#save the observed data
  obs.matrix=sim.data$obs.matrix[[q]] #observed data
  save.model$obs.matrix=obs.matrix

#The second saved list item will be a list for each of the fitted models and outputs: null, reduced, full 
  save.model$models=vector("list",3)
  names(save.model$models)=c("Null","Reduced","Full")


#Bundle data for jags
  data.list <- list(y = obs.matrix, R = dim(obs.matrix)[1], T = dim(obs.matrix)[2])

#Initial values
  zst=rep(4,dim(data.list$y)[1])
  inits <- function(){list(z = zst)}

################################################################
#Fit the Null Model 
  params <- c("psi","pdet","alpha")
#Prepare the model and data
  model.null <- jags.model(file="jags.multistate.occ.null", 
                         data = data.list,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

#Do a burn in period
  update(model.null,n.iter=nb,progress.bar="none")

#Fit the model  
  model.null.fit <- jags.samples(model.null, variable.names=params, 
                               n.iter=ni, 
                               thin=nt,
                               progress.bar="none")

#save the model  
  save.model$models$Null$model=model.null.fit

#save the parameters
  alpha.samples=model.null.fit$alpha[,,1]
  psi.samples=model.null.fit$psi[,,1]
  pdet.samples=model.null.fit$pdet[,,1]

#We need to derive state occupancy and detection probs

psiDay=psi.samples
psiNight=psi.samples
psiND=psi.samples
psi0=1-psiDay-psiNight-psiND

#psi.overall=psiDay+psiNight+psiND
#hist(psi.overall)

omega.samples=cbind(psi0,psiDay,psiNight,psiND)
#table(rowSums(omega.samples))

#save the omega samples
save.model$models$Null$omega.samples=omega.samples

psiDay.marginal.samples=psiDay+psiND
psiNight.marginal.samples=psiNight+psiND

#Derive the state-dependent detection probabilities
pDay.samples=pdet.samples
pNight.samples=pdet.samples

pND.ND=pdet.samples*pdet.samples
pND.D=(1-pdet.samples)*pdet.samples
pND.N=(1-pdet.samples)*pdet.samples
pND.0=(1-pdet.samples)*(1-pdet.samples)

pND.samples=cbind(pND.0,pND.D,pND.N,pND.ND)

det.samples=cbind(pDay.samples,pNight.samples,pND.samples)
save.model$models$Null$det.samples=det.samples


#alpha=alpha.samples,
save.model$models$Null$samples=data.frame(pDet=pdet.samples,
                                          psi=psi.samples,pDay=pdet.samples,pNight=pdet.samples,
                                          pND=pND.samples,psiDay.marginal=psiDay.marginal.samples,
                                          psiNight.marginal=psiNight.marginal.samples)


save(save.model,file=paste("data.null.3.models",q,sep=""))


######################################################
######################################################
#Fit the reduced model  

# Parameters monitored
# Note that psi and p parameters are marginal probabilities
# We need to derive the state-occurence probabilities
params <- c("psiDay","psiNight", "pNight", "pDay") 

model.reduced <- jags.model(file="jags.multistate.occ.reduced",
                            data = data.list,
                            inits=inits,
                            n.chains = nc,
                            n.adapt=adapt)

#this is a burn in period
update(model.reduced,n.iter=nb)

#now fit the model
model.reduced.fit <- jags.samples(model.reduced, variable.names=params, 
                                  n.iter=ni, 
                                  thin=nt,
                                  progress.bar="none")

#Save the fitted model
save.model$models$Reduced$model=model.reduced.fit

#Get posterior distribution samples
psiNight.marginal.samples=model.reduced.fit$psiNight[,,1]
psiDay.marginal.samples=model.reduced.fit$psiDay[,,1]
pNight.samples=model.reduced.fit$pNight[,,1]
pDay.samples=model.reduced.fit$pDay[,,1]


save.model$models$Reduced$samples=data.frame(psiNight.marginal=psiNight.marginal.samples,
                                             psiDay.marginal=psiDay.marginal.samples,
                                             pNight.marginal=pNight.samples,
                                             pDay.marginal=pDay.samples)

#Get state dependent occupancy estimates
psiND.samples=psiNight.marginal.samples*psiDay.marginal.samples  
psiNight.samples=psiNight.marginal.samples*(1-psiDay.marginal.samples)
psiDay.samples=psiDay.marginal.samples*(1-psiNight.marginal.samples)
psi0.samples=(1-psiNight.marginal.samples)*(1-psiDay.marginal.samples)

omega.samples=cbind(psi0.samples,psiDay.samples,psiNight.samples,psiND.samples)
#table(rowSums(omega.samples))
save.model$models$Reduced$omega.samples=omega.samples

#Get state-dependent detection estimates
pDay.samples=pDay.samples
pNight.samples=pNight.samples
pND.ND=pDay.samples*pNight.samples
pND.D=pDay.samples*(1-pNight.samples)
pND.N=pNight.samples*(1-pDay.samples)
pND.0=(1-pDay.samples)*(1-pNight.samples)

pND.samples=cbind(pND.0,pND.D,pND.N,pND.ND)

det.samples=cbind(pDay.samples,pNight.samples,pND.samples)
save.model$models$Reduced$det.samples=det.samples

save(save.model,file=paste("data.reduced.3.models",q,sep=""))
################################################################
################################################################
################################################################
#Fit the Full Model 

# Parameters to monitor
params <- c("psi", "pNight", "pDay","pND") 


#Prepare the model and data
model.full <- jags.model(file="jags.multistate.occ.full", 
                         data = data.list,
                         inits=inits,
                         n.chains = nc,
                         n.adapt=adapt)

#Do a burn in period
update(model.full,n.iter=nb,progress.bar="none")

#fit the model
model.full.fit <- jags.samples(model.full, variable.names=params, 
                               n.iter=ni, 
                               thin=nt,
                               progress.bar="none")

#save the model
save.model$models$Full$model=model.full.fit

#occupancy parameters- these are alrady state-specific
psi.samples=t(model.full.fit$psi[1:4,,1])

save.model$models$Full$omega.samples=psi.samples

#Derive marginal occurecnce probs
psiDay.marginal.samples=psi.samples[,2]+psi.samples[,4]
psiNight.marginal.samples=psi.samples[,3]+psi.samples[,4]


#Detection parameters
pNight.samples=model.full.fit$pNight[,,1]
pDay.samples=model.full.fit$pDay[,,1]
pND.samples=t(model.full.fit$pND[1:4,,1])

det.samples=cbind(pDay.samples,pNight.samples,pND.samples)
save.model$models$Full$det.samples=det.samples



save.model$models$Full$samples=data.frame(psi=psi.samples,pDay=pDay.samples,pNight=pNight.samples,
                                          pND=pND.samples,psiDay.marginal=psiDay.marginal.samples,
                                          psiNight.marginal=psiNight.marginal.samples)



save(save.model,file=paste("/data.reduced.3.models",q,sep=""))



print(paste(q," - ",round(q/n.sim,digits=4))  )
}#END MODEL FITTING LOOP
########################################
########################################
########################################