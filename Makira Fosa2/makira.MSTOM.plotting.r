#M1 is the most supported model
#Load the model and plots results

library(rjags)
library(runjags)
library(coda)
library(bayesplot)
library(ggplot2)


rm(list=ls())
load("Makira Fosa2/M1.full.no.covs.out")
fit <- combine.mcmc(M1.full.no.covs)
#load the prepared data file
load("Makira Fosa2/Makira.data2")
y=Makira.data2[[2]] #detection history, sites x occs x survey area

#assign covariate data to object
cov=Makira.data2[[3]]

#These are the mean (across survey areas- i.e., mean of random effect) for each logit scaled parameter for states 2, 3, and 4
#alpha 1 applies to state 2 and so on
mu.alpha.1=fit[,which(grepl("mu.alpha1",colnames(fit)))]
mu.alpha.2=fit[,which(grepl("mu.alpha2",colnames(fit)))]
mu.alpha.3=fit[,which(grepl("mu.alpha3",colnames(fit)))]

#Since there are no covariates, these are simply the intercepts on the logit scale
hist(mu.alpha.1)
hist(mu.alpha.2)

#The alpha 3 is how much different state 4 is from simply the combinations of states 2 and 3
#The distribution being a bit negative suggests that we are observing fosa a bit less than expected
#based on simply the combinations of states 2 and 3
hist(mu.alpha.3)

#Derive the mean-level (across all survey areas) state specific occupancy probability
denominator=1+exp(mu.alpha.1)+exp(mu.alpha.2)+exp(mu.alpha.1+mu.alpha.2+mu.alpha.3)
psi.day=exp(mu.alpha.1)/denominator
psi.night=exp(mu.alpha.2)/denominator
psi.ND=exp(mu.alpha.1+mu.alpha.2+mu.alpha.3)/denominator

occ.matrix=cbind(psi.day,
                 psi.night,
                 psi.ND)


#There is not much support for variation in occurence
png(file="Makira Fosa2/Makira.fosa.occ.parms.png",res=200,units = "in",height=8,width=8)
color_scheme_set("red")
mcmc_intervals(occ.matrix, pars = colnames(occ.matrix))+ 
               labs(x = "Occupancy Probability",y="State")+
               theme(text = element_text(size=20))
dev.off()

#########################################################
#########################################################
#Extract and Plot detection parameters
#These are the mean level (across all survey areas) detection parameters on the logit scale.
mu.beta1=fit[,which(grepl("mu.beta1",colnames(fit)))]
mu.beta2=fit[,which(grepl("mu.beta2",colnames(fit)))]
mu.beta3=fit[,which(grepl("mu.beta3",colnames(fit)))]
mu.beta4=fit[,which(grepl("mu.beta4",colnames(fit)))]
mu.beta5=fit[,which(grepl("mu.beta5",colnames(fit)))]

#Derive detection probability (mean across all survey areas)

pDay=exp(mu.beta1)/(1+exp(mu.beta1))
pNight=exp(mu.beta2)/(1+exp(mu.beta2))
denominator1=1+exp(mu.beta3) + exp(mu.beta4) + exp(mu.beta5) +exp(mu.beta3 + mu.beta4 + mu.beta5)
pND.0 =1/denominator1
pND.D =exp(mu.beta3)/denominator1
pND.N = exp(mu.beta4)/denominator1
pND.ND =exp(mu.beta3 + mu.beta4 + mu.beta5)/denominator1

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
lines(density(pND.N),lwd=3,col=2)
lines(density(pND.ND),lwd=3,col=3)


png(file="Makira Fosa2/Makira.fosa.det.parms.png",res=200,units = "in",height=8,width=8)
color_scheme_set("red")
mcmc_intervals(det.matrix, pars = colnames(det.matrix)[-3])+ 
  labs(x = "Detection Probability",y="State")+
     theme(text = element_text(size=20))
dev.off()

