#######################################
# Brian D. Gerber 
# 06/29/2020
#
# This script is intended to simulate multi-state occuancy data where the states are defined
# based on temporal use of a site.
#
# 4 States
# 1: No use
# 2: Day only
# 3: Night Only
# 4: Night and Day

#  ****This script simulates data from the Full model without covariates****
#
# The Full model - state detection and occurence can be uniquly different without shared relationships.

########################################
#######################################
#Setup the workspace
rm(list=ls())
setwd("C:/Users/Brian Gerber/Google Drive/MSTOM/multi.state.temporal.activity/data simulation")
setwd("C:/Users/bgerber/Google Drive/MSTOM/multi.state.temporal.activity/data simulation")

#Which file name should we save the simulated data
#The saved simulations will be saved within the working directory
filename="sim.data.multistate.full"

#How many datasets to simualte
n.sim=100

#How many sites and occasions
n.sites=200 #number of sites
n.occs=5    #number of surveys per site

#################################################
#Define the model type
model.type="Full"

#If TIF == 1 then the probability of night (state 3) and day (state 2) occurence is the product 
#of the maringal occurences of day and night states.If TIF does not equal one, it changes how 
#state 4 probability is different from the product

TIF.occ=1.1

psiNight.Marginal=0.8  #prob combination of states 3 and 4
psiDay.Marginal=0.6   #prob combination of states 2 and 4

#Define State Occurence Probabilitities
psiND=psiNight.Marginal*psiDay.Marginal*TIF.occ
psiNight=psiNight.Marginal-psiND
psiDay=psiDay.Marginal-psiND
psi0=1-psiND-psiNight-psiDay

#Alternatively  
#Define State Occurence Probabilitities  
#  psiNight=0.2  #prob combination of states 3 and 4
#  psiDay=0.1    #prob combination of states 2 and 4
#  psiND= 0.4
#  psi0=1-psiND-psiNight-psiDay


#probability of site i for each state - in order  - states 1,2,3,4
#These are mututally exclsive probabilitities and must sum to 1
omega=c(psi0,psiDay,psiNight,psiND) 
omega
sum(omega)

#################################################
#Define the detection matrix

pDay=0.4
pNight=0.6

#These do not have to be equal to the above pDay and pNight
pND.D.marginal=pDay
pND.N.marginal=pNight
TIF.det=1.3

#Define State 4 Detection Probabilitities 
pND.ND=pND.D.marginal*pND.N.marginal*TIF.det  #Detection of state 4 in state 4
pND.D=pND.D.marginal-pND.ND                   #Detection of state 2 in state 4
pND.N=pND.N.marginal-pND.ND                   #Detection of state 3 in state 4
pND.0=1-pND.ND- pND.D-pND.N                   #Detection of state 1 in state 4

#This function creates the detection matrix
source("det.matrix.func.R")
det.matrix=det.matrix.func(pNight,pDay,pND.N,pND.D,pND.0)
det.matrix

#Check to make sure the sums of the rows (true states) equal one and no negative values
rowSums(det.matrix)

length(which(det.matrix<0))

#################################################
#Setup storage of simulation information
names.list=c("Model","n.sim","n.sites","n.occs","omega","marginal.probs","TIF.occ","TIF.det","det.matrix",
             "state.occ.matrix","state.occurence","obs.matrix")

sim.data=vector("list",length(names.list))
names(sim.data)=names.list

#Save defined inputs of simulation
sim.data$Full.Model=model.type
sim.data$n.sim=n.sim
sim.data$n.sites=n.sites
sim.data$n.occs=n.occs
sim.data$omega=omega
sim.data$TIF.occ=TIF.occ
sim.data$TIF.det=TIF.det
sim.data$marginal.probs=data.frame(psiNight.Marginal=psiNight.Marginal,psiDay.Marginal=psiDay.Marginal)
sim.data$det.matrix=det.matrix

#True occurence and observed data change for each iteration of n.sim
sim.data$state.occ.matrix=vector("list",n.sim)
sim.data$state.occurence=vector("list",n.sim)
sim.data$obs.matrix=vector("list",n.sim)

#start data simulation loop
for(z in 1:n.sim){
  
  #Simualte true occurence states using the omega probs
  state.occurence.matrix=rmultinom(n.sites,1,prob=omega)
  rownames(state.occurence.matrix)=rownames(det.matrix)
  
  #reducd matrix to simple vector format
  state.occurence=apply(state.occurence.matrix,2,FUN=function(x){which(x==1)})
  
  #save true occurences
  sim.data$state.occ.matrix[[z]]=state.occurence.matrix
  sim.data$state.occurence[[z]]=state.occurence
  
  #Simulate the detection process
  obs.matrix=matrix(0,nrow=n.sites,ncol=n.occs)
  for(i in 1:n.sites){
    for(j in 1:n.occs){
      temp=rmultinom(1,1,state.occurence.matrix[,i]%*%det.matrix)
      obs.matrix[i,j]= apply(temp,2,FUN=function(x){which(x==1)})
    }
  }
  
  #Here is the data we observe, n.sites by n.occs
  
  sim.data$obs.matrix[[z]]=obs.matrix
  print(z)
} ##END SIMULATION LOOP
########################################
########################################

#save sim.data to the correct directory
save(sim.data,file=paste(filename))
