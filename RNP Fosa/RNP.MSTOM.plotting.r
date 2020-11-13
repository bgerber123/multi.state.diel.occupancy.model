#M1 is the most supported model
#Load the model and plots results

library(rjags)
library(runjags)
library(coda)
library(bayesplot)
library(ggplot2)


#for first fit
load("RNP Fosa/M1.fit")

#rm(list=ls())
#load("RNP Fosa/M3.red.out")
#fit <- combine.mcmc(M3.red)
#load the prepared data file
load("RNP Fosa/RNP2.data")
covs=RNP2.data[[2]]
cov=covs$DistTown
hist(cov)
cov1.unscaled=as.numeric(cov)
cov1.scaled=as.numeric(scale(cov,center = TRUE,scale = TRUE))


#Plot detection parameters
pDay=fit[,which(grepl("pDay",colnames(fit)))]
pNight=fit[,which(grepl("pNight",colnames(fit)))]

pND.0=fit[,which(grepl("pND",colnames(fit)))[1]]
pND.D=fit[,which(grepl("pND",colnames(fit)))[2]]
pND.N=fit[,which(grepl("pND",colnames(fit)))[3]]
pND.ND=fit[,which(grepl("pND",colnames(fit)))[4]]

det.matrix=cbind(pDay,
                 pNight,
                 pND.0,
                 pND.D,
                 pND.N,
                 pND.ND)

par(mfrow=c(1,2))
plot(density(pDay),xlim=c(0,0.2),lwd=3)
lines(density(pNight),xlim=c(0,0.2),lwd=3,col=2)
plot(density(pND.D),xlim=c(0,0.2),lwd=3)
lines(density(pND.N),xlim=c(0,0.2),lwd=3,col=2)
lines(density(pND.ND),xlim=c(0,0.2),lwd=3,col=3)


png(file="RNP Fosa/RNP.fosa.det.parms.png",res=200,units = "in",height=8,width=8)
color_scheme_set("red")
mcmc_intervals(det.matrix, pars = colnames(det.matrix)[-3],xlab="Detection Parameters")+ labs(x = "Detection Probability",y="State")+
     theme(text = element_text(size=20))
dev.off()

colnames(fit)
#The first 6 columns are the alpha parameters
fit.matrix=as.matrix(fit)

png(file="RNP Fosa/RNP.fosa.parms.png",res=200,units = "in",height=8,width=8)
color_scheme_set("red")
mcmc_areas(fit.matrix,
           pars = c("alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]","alpha[5]","alpha[6]"),
           prob = 0.5,
           prob_outer=0.99) + geom_vline(xintercept=0, linetype="solid", 
                                         color = "black", size=1)+
  theme(text = element_text(size=20))
dev.off()

length(which(fit.matrix[,2]>0))/dim(fit.matrix)[1]
length(which(fit.matrix[,4]<0))/dim(fit.matrix)[1]
length(which(fit.matrix[,6]>0))/dim(fit.matrix)[1]


#################################################
#Need to predict site occupancy using newcov

n.mcmc=dim(fit.matrix)[1]
n.sites=dim(y)[1]



alpha1=fit[,which(grepl("alpha[1]",colnames(fit), fixed = TRUE))]
alpha2=fit[,which(grepl("alpha[2]",colnames(fit), fixed = TRUE))]
alpha3=fit[,which(grepl("alpha[3]",colnames(fit), fixed = TRUE))]
alpha4=fit[,which(grepl("alpha[4]",colnames(fit), fixed = TRUE))]
alpha5=fit[,which(grepl("alpha[5]",colnames(fit), fixed = TRUE))]
alpha6=fit[,which(grepl("alpha[6]",colnames(fit), fixed = TRUE))]


#cov1.scaled2=cov1.scaled[order(cov1.scaled)]
#cov1.unscaled2=cov1.unscaled[order(cov1.unscaled)]

x.pred=seq(500,4000,length.out = 100)
x.pred.scaled=as.numeric(scale(x.pred,center = TRUE,scale = TRUE))
psi.state.mcmc=array(NA, dim=c(length(x.pred),n.mcmc,4))
dim(psi.state.mcmc)
#Get psi estimates for each state by site (for all n.mcmc)
for (s in 1:length(x.pred)){
  phi1 <- 1
  phi2 <- exp(alpha1+alpha2*x.pred.scaled[s])
  phi3 <- exp(alpha3+alpha4*x.pred.scaled[s])
  phi4 <- exp(alpha1+alpha2*x.pred.scaled[s]+
                alpha3+alpha4*x.pred.scaled[s]+
                alpha5+alpha6*x.pred.scaled[s])
  constant=phi1+phi2+phi3+phi4
  psi.state.mcmc[s,,1]     <- phi1/constant
  psi.state.mcmc[s,,2]     <- phi2/constant
  psi.state.mcmc[s,,3]     <- phi3/constant
  psi.state.mcmc[s,,4]     <- phi4/constant
  
}

#site 1 and state 1
hist(psi.state.mcmc[1,,3])

#xpred, mcmc, 4 states
dim(psi.state.mcmc)

#We want to plot cov1.unscaled (x axis) by prob for each state

y <- seq(0, 1, length=1000)
#Store matrix values and loop through and calcaulte the density
#of each value of probability at each level of weight
z1=z2=z3=z4 <- matrix(nrow=length(x.pred), ncol=length(y))

i=1
y=density(psi.state.mcmc[i,,4],n = 1000, from=0, to=1)$x

for(i in 1:length(x.pred)){ z1[i,] <- density(psi.state.mcmc[i,,1],n = 1000, from=0, to=1,adjust=1.25)$y}
for(i in 1:length(x.pred)){ z2[i,] <- density(psi.state.mcmc[i,,2],n = 1000, from=0, to=1,adjust=1.25)$y}
for(i in 1:length(x.pred)){ z3[i,] <- density(psi.state.mcmc[i,,3],n = 1000, from=0, to=1,adjust=1.25)$y}
for(i in 1:length(x.pred)){ z4[i,] <- density(psi.state.mcmc[i,,4],n = 1000, from=0, to=1,adjust=1.25)$y}

z1.center=apply(psi.state.mcmc[,,1],1,median)
z2.center=apply(psi.state.mcmc[,,2],1,median)
z3.center=apply(psi.state.mcmc[,,3],1,median)
z4.center=apply(psi.state.mcmc[,,4],1,median)

#Plot the posterior distribution where the shading represents higher/lower
#probability within a value of weight
library(denstrip)
nlevel=200
xlim1=c(min(x.pred),max(x.pred))

png(file="RNP Fosa/RNP.fosa.State.Prob.png",res=200,units = "in",height=10,width=10)
par(mfrow=c(2,2))
plot(0, type="n", ylim=c(0, 1),xlim=c(xlim1[1],xlim1[2]),ylab="Probability of Occurence",xlab="Distance to Nearest Village (m)",main="State 1 (No Use)")
densregion(x.pred,y,z1,nlevels = nlevel,colmax = "black",pointwise = TRUE)
lines(x.pred,z1.center,col="white",lwd=2)

plot(0, type="n", ylim=c(0, 1),xlim=c(xlim1[1],xlim1[2]),ylab="Probability of Occurence",xlab="Distance to Nearest Village (m)",main="State 2 (Day Use)")
densregion(x.pred,y,z2,nlevels = nlevel,colmax = "darkgreen",pointwise = TRUE)
lines(x.pred,z2.center,col="white",lwd=2)


plot(0, type="n", ylim=c(0, 1),xlim=c(xlim1[1],xlim1[2]),ylab="Probability of Occurence",xlab="Distance to Nearest Village (m)",main="State 3 (Night Use)")
densregion(x.pred,y,z3,nlevels = nlevel,colmax = "darkblue",pointwise = TRUE)
lines(x.pred,z3.center,col="white",lwd=2)


plot(0, type="n", ylim=c(0, 1),xlim=c(xlim1[1],xlim1[2]),ylab="Probability of Occurence",xlab="Distance to Nearest Village (m)",main="State 4 (Day & Night Use)")
densregion(x.pred,y,z4,nlevels = nlevel,colmax = "darkred",pointwise = TRUE)
lines(x.pred,z4.center,col="white",lwd=2)

dev.off()


