#Plotting the most supported model
library(runjags)
library(bayesplot)
library(ggplot2)
rm(list=ls())
load("AJB Fosa/M1.full.fit.plotting")

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


png(file="AJB Fosa/AJB.fosa.det.parms.png",res=200,units = "in",height=8,width=8)
color_scheme_set("red")
mcmc_intervals(det.matrix, pars = colnames(det.matrix)[-3])+ labs(x = "Detection Probability",y="State")+
  theme(text = element_text(size=20))
dev.off()


head(fit)

prob.state1=fit[,match("prob[1,1]",colnames(fit))]
prob.state2=fit[,match("prob[1,2]",colnames(fit))]
prob.state3=fit[,match("prob[1,3]",colnames(fit))]
prob.state4=fit[,match("prob[1,4]",colnames(fit))]

prob.matrix=cbind(prob.state1,prob.state2,prob.state3,prob.state4)
hist(prob.state4)

png(file="AJB Fosa/AJB.fosa.state.occ.png",res=200,units = "in",height=8,width=8)
color_scheme_set("red")
mcmc_intervals(prob.matrix)+ scale_y_discrete(labels = c("State 1 (No Use)","State 2 (Day)","State 3 (Night)","State 4 (Night&Day)"))+
  labs(x = "Occupancy Probability",y="State")+
  theme(text = element_text(size=20))
dev.off()

