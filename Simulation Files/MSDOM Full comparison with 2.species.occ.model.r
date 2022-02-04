#Comparison of the original unconditional two species occupancy model (MacKenzie et al. 2004)
#with that of the MSDOM Full model. These models are equivalent for 4 states (2-species).
#However, the Bayesian MSTOM does not have the same estimation optimization issues, with
#or without covariates.

#MacKenzie, D. I., Bailey, L. L., & Nichols, J. D. (2004). Investigating species co‐occurrence patterns when species are detected imperfectly. Journal of Animal Ecology, 73(3), 546-555.

#The 2-species occupancy model is fit in MARK. See,
#https://sites.warnercnr.colostate.edu/gwhite/occupancy-estimation-two-species/

#This model can only be fit in MARK. It can not be fit in RMark, see
#http://www.phidot.org/forum/viewtopic.php?f=21&t=3995&p=13211&hilit=mlogit#p13211

#This model is also not setup to be fit in RPresence.

#Setup the workspace and load the first simulated data set
rm(list=ls())
load("./Simulation Files/sim.full.data")
length(sim.data$obs.matrix)
q=1
y=sim.data$obs.matrix[[1]]

#Need to turn this matrix into a 2-species occupancy
#detection matrix that is suitable for MARK
n.occs=dim(y)[2]
n.sites=dim(y)[1]
night.mat=matrix(0,n.sites,n.occs)
day.mat=matrix(0,n.sites,n.occs)
joint=matrix(0,n.sites,n.occs*2) #This is for MARK
for(i in 1:n.sites){
  for(j in 1:n.occs){
    if(y[i,j]==1){night.mat[i,j]=0;day.mat[i,j]=0}
    if(y[i,j]==2){night.mat[i,j]=1;day.mat[i,j]=0}
    if(y[i,j]==3){night.mat[i,j]=0;day.mat[i,j]=1}
    if(y[i,j]==4){night.mat[i,j]=1;day.mat[i,j]=1}
  }
  joint[i,]=c(rbind(night.mat[i,], day.mat[i,]))  
}

night.mat
day.mat
joint

#Create the CH inp file for MARK
mark_input=as.data.frame(apply(joint,1,paste,collapse=""))
mark_input <- data.frame(lapply(mark_input, as.character), stringsAsFactors=FALSE)

#Output the inp file for MARK
mark.out=cbind(mark_input,rep("1;",dim(mark_input)[1]))
colnames(mark.out)= NULL
mark.out
write.table(mark.out,file="./Simulation Files/MARK.2.species/markout1.inp",sep=" ",quote = FALSE,row.names = FALSE)


#Go to MARK and fit 2-species occupancy model; see the folder MARK.2.species
#Note, an mlogit(1) link funciton is needed for rAB,rAb, and raB parameters.
#Also, a design matrix that constrains beta parameters is required for the optimization
#to correctly return the MLE's. An identity matrix will not work.

#Note, the MSTOM Reduced model could be fit if MARK allows parameters to be functions
#of other parameters. I do not know how to do this.

###############################################################
#Read in the saved real parameters estimated from MARK. 
#A constant model was fit (no site or temporal variation,
#thus correctly specifying how the data were simulated).
MARK.est=read.table("./Simulation Files/MARK.2.species/real.parms.txt")
colnames(MARK.est)=c("Parm","MLE","SE","LCL","UCL")

#Derive state-specific probabilities, 
#Note that A=Day, B= Night. Also MARK estimates the marginal probabilityies
#of Day and Night occurence. One needs to derive state-specific occupancy
#pobabilities by subtracting the probabilty of night&Day
PSI.AB=MARK.est[1,2]
PSI.A=MARK.est[2,2]-PSI.AB
PSI.B=MARK.est[3,2]-PSI.AB

#Get the detection probabilities
pA=MARK.est[4,2]
pB=MARK.est[5,2]
rAB=MARK.est[6,2]
rAb=MARK.est[7,2]
raB=MARK.est[8,2]
##########################################
##########################################

#Load Bayesian MSTOM FUll model fit
load("./Simulation Files/fit.simdata.Full.1.out")

#Assign posterior samples to objects
omega.samples=save.model$models$Full$omega.samples
det.samples=save.model$models$Full$det.samples


#Plot and compare posterior distribution (occupancy) with the MLE's from MARK
hist(omega.samples[,2],main="Day Occurence",xlab="Proabability",freq = FALSE)
abline(v=PSI.A,col=2,lty=3,lwd=3)

hist(omega.samples[,3],main="Night Occurence",xlab="Proabability",freq = FALSE)
abline(v=PSI.B,col=2,lty=3,lwd=3)

hist(omega.samples[,4],main="Day & Night Occurence",xlab="Proabability",freq = FALSE)
abline(v=PSI.AB,col=2,lty=3,lwd=3)

#Plot and compare posterior distribution (detection) with the MLE's from MARK
head(det.samples)

hist(det.samples[,1],main="Day Detection",xlab="Proabability",freq = FALSE)
abline(v=pA,col=2,lty=3,lwd=3)

hist(det.samples[,2],main="Night Detection",xlab="Proabability",freq = FALSE)
abline(v=pB,col=2,lty=3,lwd=3)

hist(det.samples[,6],main="ND.ND Detection",xlab="Proabability",freq = FALSE)
abline(v=rAB,col=2,lty=3,lwd=3)

hist(det.samples[,4],main="ND.D Detection",xlab="Proabability",freq = FALSE)
abline(v=rAb,col=2,lty=3,lwd=3)

hist(det.samples[,5],main="ND.N Detection",xlab="Proabability",freq = FALSE)
abline(v=raB,col=2,lty=3,lwd=3)
