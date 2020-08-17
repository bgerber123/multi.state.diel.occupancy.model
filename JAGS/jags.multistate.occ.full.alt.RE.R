    model{ 
    # Priors
      #Priors for mean population level intercepts and effects
      mu.alpha1~dnorm(0,0.5)
      mu.alpha2~dnorm(0,0.5)
      mu.alpha3~dnorm(0,0.5)
      mu.beta1~dnorm(0,0.5)
      mu.beta2~dnorm(0,0.5)
      mu.beta3~dnorm(0,0.5)
      mu.beta4~dnorm(0,0.5)
      mu.beta5~dnorm(0,0.5)
      #Priors for the variance of the population level effects
      s2.alpha1~dgamma(1,1)
      s2.alpha2~dgamma(1,1)
      s2.alpha3~dgamma(1,1)
      s2.beta1~dgamma(1,1)
      s2.beta2~dgamma(1,1)
      s2.beta3~dgamma(1,1)
      s2.beta4~dgamma(1,1)
      s2.beta5~dgamma(1,1)
      #inverse of variance
      tau.alpha1<-1/s2.alpha1
      tau.alpha2<-1/s2.alpha2
      tau.alpha3<-1/s2.alpha3
      tau.beta1<-1/s2.beta1
      tau.beta2<-1/s2.beta2
      tau.beta3<-1/s2.beta3
      tau.beta4<-1/s2.beta4
      tau.beta5<-1/s2.beta5
    # Define state vector
    for (u in 1:Ns){
      #These the population-level random effects for the intercepts and effects for
      #the logit-linear occupancy model (b0, b1) and detection model (g0 and g1)
      alpha1[u]~dnorm(mu.alpha1, tau.alpha1)
      alpha2[u]~dnorm(mu.alpha2, tau.alpha2)
      alpha3[u]~dnorm(mu.alpha3, tau.alpha3)
    for (i in 1:N){
      phi[i,1,u] <- 1
      phi[i,2,u] <- exp(alpha1[u])
      phi[i,3,u] <- exp(alpha2[u])
      phi[i,4,u] <- exp(alpha1[u]+alpha2[u]+alpha3[u])
      PSI[i,1,u] <- phi[i,1,u]/sum(phi[i,,u])
      PSI[i,2,u] <- phi[i,2,u]/sum(phi[i,,u])
      PSI[i,3,u] <- phi[i,3,u]/sum(phi[i,,u])
      PSI[i,4,u] <- phi[i,4,u]/sum(phi[i,,u])
      }
    }
    # Define observation matrix
    # Order of indices: true state, time, observed state, survey
    for (u in 1:Ns){
      beta1[u]~dnorm(mu.beta1, tau.beta1)
      beta2[u]~dnorm(mu.beta2, tau.beta2)
      beta3[u]~dnorm(mu.beta3, tau.beta3)
      beta4[u]~dnorm(mu.beta4, tau.beta4)
      beta5[u]~dnorm(mu.beta5, tau.beta5)
    for (j in 1:K){
    q[1,j,1,u] <- 1
    q[1,j,2,u] <- 0
    q[1,j,3,u] <- 0
    q[1,j,4,u] <- 0
    q[2,j,1,u] <- 1
    q[2,j,2,u] <- exp(beta1[u])
    q[2,j,3,u] <- 0
    q[2,j,4,u] <- 0
    q[3,j,1,u] <- 1
    q[3,j,2,u] <- 0
    q[3,j,3,u] <- exp(beta2[u])
    q[3,j,4,u] <- 0
    q[4,j,1,u] <- 1
    q[4,j,2,u] <- exp(beta3[u])
    q[4,j,3,u] <- exp(beta4[u])
    q[4,j,4,u] <- exp(beta3[u] + beta4[u] + beta5[u])
    }
      }
    # State-space likelihood
    # State equation: model of true states (z)
    for (u in 1:Ns){
    for (i in 1:N){
     z[i,u] ~ dcat(phi[i,,u])
    }
    }
    # Observation equation
    for (u in 1:Ns){
    for (i in 1:N){
       for (j in 1:K){ 
        y[i,j,u] ~ dcat(q[z[i,u],j,,u])
       } #j
    } #i
    }#u

#Derive outputs      
for (u in 1:Ns){      
      psiDay[u] <- PSI[1,2,u]  
      psiNight[u] <- PSI[1,3,u]
      psiND[u] <- PSI[1,4,u]
      #Derive detection parameters for occasion 1
      pDay[u] <- q[2,1,2,u]/sum(q[2,1,,u])
      pNight[u] <- q[3,1,3,u]/sum(q[3,1,,u])
      pND.ND[u] <- q[4,1,4,u]/sum(q[4,1,,u])
      pND.N[u] <- q[4,1,3,u]/sum(q[4,1,,u])
      pND.D[u] <- q[4,1,2,u]/sum(q[4,1,,u])
      pND.0[u] <- q[4,1,1,u]/sum(q[4,1,,u])      
}
#Mean population-level occurence
for(i in 1:N){
  phiM[i,1] <- 1
  phiM[i,2] <- exp(mu.alpha1)
  phiM[i,3] <- exp(mu.alpha2)
  phiM[i,4] <- exp(mu.alpha1+mu.alpha2+mu.alpha3)
  PSIM[i,1] <- phiM[i,1]/sum(phiM[i,])
  PSIM[i,2] <- phiM[i,2]/sum(phiM[i,])
  PSIM[i,3] <- phiM[i,3]/sum(phiM[i,])
  PSIM[i,4] <- phiM[i,4]/sum(phiM[i,])
}
      
}