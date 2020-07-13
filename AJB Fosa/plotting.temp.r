#plotting for investigation of models
library(runjags)
library(bayesplot)
library(ggplot2)
load("AJB Fosa/M1.full.covs.1")
fit <- combine.mcmc(M1.full.covs.1)
head(fit)
head(survey.cov)
survey.cov[,6]

#row 1 is effect diff from grand mean of survey 2008
#row 25 is effect diff from grand mean of survey 2010
#row 49 is effect diff from grand mean of survey 2011
#row 73 is effect diff from grand mean of survey 2012
#row 97 is effect diff from grand mean of survey 2013
#effect diff from grand of survey 2015 is grand mean (intercept) minus all coefs


colnames(fit)[1:20]


hist(fit[,7],breaks=20)

hist(plogis(fit[,7]))


###########z##########################
load("AJB Fosa/M1.full.no.covs")
fit <- combine.mcmc(M1.full.no.covs)

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


png(file="AJB Fosa/B.fosa.det.parms.png",res=200,units = "in",height=8,width=8)
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
#################################################

load("AJB Fosa/M1.full.covs.2")
fit <- combine.mcmc(M1.full.covs.2)
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


fit.matrix=as.matrix(fit)

colnames(fit)[1:10]

length(which(fit.matrix[,2]<0))/dim(fit.matrix)[1]
length(which(fit.matrix[,4]<0))/dim(fit.matrix)[1]


