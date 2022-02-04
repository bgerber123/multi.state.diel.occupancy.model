#######################################
# Brian D. Gerber 
# 06/29/2020
#
# This script is intended to simulate multi-state occupancy data where the states are defined
# based on temporal use of a site.
#
# 4 States
# 1: No use
# 2: Day only
# 3: Night Only
# 4: Night and Day

#  ****This script simulates data from the Reduced model without covariates****

# The Reduced model - assume that the probability of detection and occurrence in state 4 is a combination
#                  of probabilities associated with states 2 and 3.
########################################
########################################
#Setup the workspace
rm(list=ls())

#Which file name should we save the simulated data
#The saved simulations will be saved within the working directory
filename="./Simulation Files/sim.reduced.data"

#How many datasets to simualte
n.sim=100

#How many sites and occasions
n.sites=200 #number of sites
n.occs=5    #number of surveys per site

#################################################

#Define the model type
model.type="Reduced"

#For the reduced model, the probabilty of Night&Day (state 4) are simply the product of the 
#marginal probabilities of day and night occurence 

psiNight.Marginal=0.6  #prob combination of states 3 and 4

psiDay.Marginal=0.3    #prob combination of states 2 and 4

#Define State Occurence Probabilitities
psiND=psiNight.Marginal*psiDay.Marginal

#These can be derived this way or the way below
psiNight=psiNight.Marginal*(1-psiDay.Marginal)
psiDay=psiDay.Marginal*(1-psiNight.Marginal)
psi0=(1-psiDay.Marginal)*(1-psiNight.Marginal)

#psiNight=psiNight.Marginal-psiND
#psiDay=psiDay.Marginal-psiND
#psi0=1-psiNight-psiDay-psiND


#Omega are the probabilities of site i for each state - in order  - states 1,2,3,4
#These are mututally exclsive probabilitities and must sum to 1
omega=c(psi0,psiDay,psiNight,psiND) 
omega

#make sure probabilities sum to one
sum(omega)

#################################################
#Define the detection matrix

#If pDay and pNight are the same then there is no heterogeneity in site level of detection
#due to temporal use. They do not have to be the same though.
pDay=0.6
pNight=0.4

#Same goes for these.    
pND.D.marginal=pDay
pND.N.marginal=pNight


#Define State 4 Detection Probabilities 
pND.0=(1-pND.D.marginal)*(1-pND.N.marginal)   #Detection of state 1 in state 4
pND.D=pND.D.marginal*(1-pND.N.marginal)       #Detection of state 2 in state 4
pND.N=pND.N.marginal*(1-pND.D.marginal)       #Detection of state 3 in state 4
pND.ND=pND.D.marginal*pND.N.marginal          #Detection of state 4 in state 4

#This function creates the detection matrix
source("./Simulation Files/det.matrix.func.r")
det.matrix=det.matrix.func(pNight,pDay,pND.N,pND.D,pND.ND,pND.0)
det.matrix

#Check to make sure the sums of the rows (true states) equal one and no negative values
rowSums(det.matrix)

length(which(det.matrix<0))

#################################################
#Setup storage of simulation information
names.list=c("Model","n.sim","n.sites","n.occs","omega","marginal.probs","det.matrix",
             "state.occ.matrix","state.occurence","obs.matrix")

sim.data=vector("list",length(names.list))
names(sim.data)=names.list

#Save defined inputs of simulation
sim.data$Full.Model=model.type
sim.data$n.sim=n.sim
sim.data$n.sites=n.sites
sim.data$n.occs=n.occs
sim.data$omega=omega
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
