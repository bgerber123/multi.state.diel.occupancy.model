#Comparison of the conditional two species occupancy (Richmond et al. 2010) model with that of the MSDOM Full model. These models are equivalent for 4 states (2-species).

#Richmond, O. M., Hines, J. E., & Beissinger, S. R. (2010). Two-species occupancy models:
#   a new parameterization applied to co-occurrence of secretive rails. Ecological Applications,
#   20(7), 2036-2046.

#This codes was graciously prepared by Jim Hines.

#The 2-species occupancy model is fit in MARK. See,
#https://sites.warnercnr.colostate.edu/gwhite/occupancy-estimation-two-species/

#Setup the workspace and load the first simulated data set
rm(list=ls())
setwd("./multi.state.temporal.activity/Simulation Files")
load("sim.full.data")
#setwd('Simulation Files/MARK.2.species/')
length(sim.data$obs.matrix)

S1='proc title mstom test data;
proc chmatrix occasions=5 groups=1 etype=2SpecConOccup mixtures=2 Nodes=101 hist=300;
glabel(1)=Group 1; time interval 1 1 1 1;'

S2='proc estimate link=Sin varest=2ndPart ;
model={psiA,psiBA,psiBa,pA,pB,rA,rBA,rBa};
group=1 PsiA rows=1 cols=1 Square Constant=1;
group=1 PsiBA rows=1 cols=1 Square Constant=2;
group=1 PsiBa rows=1 cols=1 Square Constant=3;
group=1 pA rows=1 cols=5 Square Constant=4;
group=1 pB rows=1 cols=5 Square Constant=5;
group=1 rA rows=1 cols=5 Square Constant=6;
group=1 rBA rows=1 cols=5 Square Constant=7;
group=1 rBa rows=1 cols=5 Square Constant=8;
design matrix constraints=8 covariates=8 identity;
blabel(1)=PsiA;
blabel(2)=PsiBA;
blabel(3)=PsiBa;
blabel(4)=pA;
blabel(5)=pB;
blabel(6)=rA;
blabel(7)=rBA;
blabel(8)=rBa;
rlabel(1)=PsiA;
rlabel(2)=PsiBA;
rlabel(3)=PsiBa;
rlabel(4)=pA;
rlabel(5)=pB;
rlabel(6)=rA;
rlabel(7)=rBA;
rlabel(8)=rBa;
dlabel(1)=SIF;
dlabel(2)=PsiB;
dlabel(3)=PsiAB;
proc stop;
'

simrslt=NULL
#                                psiAB
TIF.occ=1.1           # TIF = -----------   psiAB = psiA*psiB*TIF
#                               psiA*psiB

psiNight.Marginal=0.8  #prob combination of states 3 and 4  (psiB)
psiDay.Marginal=0.6   #prob combination of states 2 and 4    (psiA)
#Define State Occurence Probabilitities
psiND=psiNight.Marginal*psiDay.Marginal*TIF.occ       #   (psiAB)
psiNight=psiNight.Marginal-psiND                      #   (psiBa)
psiDay=psiDay.Marginal-psiND                          #   (psiAb)
psi0=1-psiND-psiNight-psiDay                          #   1-psiAB-psiBa-psiAb

pDay=0.4
pNight=0.6

#These do not have to be equal to the above pDay and pNight
pND.D.marginal=pDay                                  #  (rA)
pND.N.marginal=pNight                                #  (rB)
TIF.det=1.3

#Define State 4 Detection Probabilities
pND.ND=pND.D.marginal*pND.N.marginal*TIF.det  #Detection of state 4 in state 4
pND.D=pND.D.marginal-pND.ND                   #Detection of state 2 in state 4
pND.N=pND.N.marginal-pND.ND                   #Detection of state 3 in state 4
pND.0=1-pND.ND- pND.D-pND.N                   #Detection of state 1 in state 4

#      psiAB,  psiA            psiB             pA    pB   rAB     rA             rB
inpt=c(psiND,psiDay.Marginal,psiNight.Marginal,pDay,pNight,pND.ND,pND.D.marginal,pND.N.marginal)

#loop over simulated data
for (q in 1:length(sim.data$obs.matrix)) {

   y=sim.data$obs.matrix[[q]]
   joint=c('00','10','01','11')[y] # convert 4 state codes into 2-digit codes for MARK
   joint=matrix(joint,ncol=ncol(y)) # convert resulting vector back into a matrix

   #  collapse matrix into histories and add frequency, '1;'
   mark.inp=paste(apply(joint,1,paste,collapse=''),'1;')

   #Go to MARK and fit 2-species occupancy model; see the folder MARK.2.species
   #Note, an mlogit(1) link funciton is needed for rAB,rAb, and raB parameters.
   #Also, a design matrix that constrains beta parameters is required for the optimization
   #to correctly return the MLE's. An identity matrix will not work.

   fname=paste0("markout",q,".inp")
   write(c(S1,mark.inp,S2),file=fname)
   oname=paste0("markout",q,".out")
   #Note, the MSTOM Reduced model could be fit if MARK allows parameters to be functions
   #of other parameters. I do not know how to do this.
   system(paste0('c:/progra~2/mark/mark64 i=',fname,' l=',oname,' lines=0'),
          show.output.on.console=T)
   ###############################################################
   #Read in the saved real parameters estimated from MARK.
   #A constant model was fit (no site or temporal variation,
   #thus correctly specifying how the data were simulated).
   i=grep('Real Function',readLines(oname))
   #MARK.est=read.table("simulation study/MARK.2.species/real.parms.txt")
   MARK.est=read.table(oname,skip=i+3,nrows=8)
   colnames(MARK.est)=c("Parm","MLE","SE","LCL","UCL")

   #Derive state-specific probabilities,
   #Note that A=Day, B= Night. Also MARK estimates the marginal probabilityies
   #of Day and Night occurence. One needs to derive state-specific occupancy
   #pobabilities by subtracting the probabilty of night&Day

   # get ests from cond. parm.       #  get ests from uncond. param.
   PSI.A = MARK.est[1,2]             # PSI.AB=MARK.est[1,2]
   PSI.BA = MARK.est[2,2]            # PSI.A=MARK.est[2,2]-PSI.AB
   PSI.Ba = MARK.est[3,2]            # PSI.B=MARK.est[3,2]-PSI.AB

   #Get the detection probabilities
   pA=MARK.est[4,2]
   pB=MARK.est[5,2]
   rA = MARK.est[6,2]                # rAB=MARK.est[6,2]
   rBA = MARK.est[7,2]               # rAb=MARK.est[7,2]
   rBa = MARK.est[8,2]               # raB=MARK.est[8,2]

   simrslt=rbind(simrslt,c(PSI.A*PSI.BA,                 # psiAB
                           PSI.A,                        # psiA
                           PSI.A*PSI.BA+(1-PSI.A)*PSI.Ba,# psiB
                           pA,pB,                        #   pA, pB
                           rA*rBA,                       #   rAB
                           rA,                           #   rA
                           rA*rBA+(1-rA)*rBa))           #   rB
   ##########################################
   ##########################################
   if (FALSE) {
      #Load Bayesian MSTOM FUll model fit
      load("Simulation Files/fit.simdata.Full.1.out")

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
   }
}

v=cbind(inpt,colMeans(simrslt))
rownames(v)=c('psiAB','psiA','psiB','pA','pB','rAB','rA','rB')
#compare known parameters (1st columns) with means of estimated parameters
print(v)
