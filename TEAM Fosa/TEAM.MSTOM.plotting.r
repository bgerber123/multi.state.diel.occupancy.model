#M1 is the most supported model
#Load the model and plots results

library(rjags)
library(runjags)
library(coda)
library(bayesplot)
library(ggplot2)


load("TEAM Fosa/M1.full.no.covs.out")

fit=combine.mcmc(M1.full.no.covs)

#load the prepared data file
#load("RNP Fosa/RNP.data")
#covs=RNP.data[[2]]
#cov=covs$DistTown
#hist(cov)
#cov1.unscaled=as.numeric(cov)
#cov1.scaled=as.numeric(scale(cov,center = TRUE,scale = TRUE))


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


png(file="TEAM Fosa/RNP.fosa.det.parms.png",res=200,units = "in",height=8,width=8)
color_scheme_set("red")
mcmc_intervals(det.matrix, pars = colnames(det.matrix)[-3],xlab="Detection Parameters")+ labs(x = "Detection Probability",y="State")+
     theme(text = element_text(size=20))
dev.off()

colnames(fit)
#The first 6 columns are the alpha parameters
fit.matrix=as.matrix(fit)

    index.occ1=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("PSI[1,1]"), x, fixed = TRUE)}))
    index.occ2=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("PSI[1,2]"), x, fixed = TRUE)}))
    index.occ3=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("PSI[1,3]"), x, fixed = TRUE)}))
    index.occ4=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("PSI[1,4]"), x, fixed = TRUE)}))

fit.psi=fit.matrix[,c(index.occ1,index.occ2,index.occ3,index.occ4)]
    
png(file="TEAM Fosa/RNP.fosa.parms.png",res=200,units = "in",height=8,width=8)
color_scheme_set("red")
mcmc_areas(fit.psi,
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
n.sites=53



alpha1=fit[,1]
alpha2=fit[,2]
alpha3=fit[,3]
alpha4=fit[,4]


