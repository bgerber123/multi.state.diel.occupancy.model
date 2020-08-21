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
    for (t in 1:T){
      #These the population-level random effects for the intercepts 
      alpha1[t]~dnorm(mu.alpha1, tau.alpha1)
      alpha2[t]~dnorm(mu.alpha2, tau.alpha2)
      alpha3[t]~dnorm(mu.alpha3, tau.alpha3)
    for (i in 1:N[t]){
      phi[i,1,t] <- 1
      phi[i,2,t] <- exp(alpha1[t])
      phi[i,3,t] <- exp(alpha2[t])
      phi[i,4,t] <- exp(alpha1[t]+alpha2[t]+alpha3[t])
      PSI[i,1,t] <- phi[i,1,t]/sum(phi[i,,t])
      PSI[i,2,t] <- phi[i,2,t]/sum(phi[i,,t])
      PSI[i,3,t] <- phi[i,3,t]/sum(phi[i,,t])
      PSI[i,4,t] <- phi[i,4,t]/sum(phi[i,,t])
      #Mean population-level occurence
      phiM[i,1,t] <- 1
      phiM[i,2,t] <- exp(mu.alpha1)
      phiM[i,3,t] <- exp(mu.alpha2)
      phiM[i,4,t] <- exp(mu.alpha1+mu.alpha2+mu.alpha3)
      PSIM[i,1,t] <- phiM[i,1,t]/sum(phiM[i,,t])
      PSIM[i,2,t] <- phiM[i,2,t]/sum(phiM[i,,t])
      PSIM[i,3,t] <- phiM[i,3,t]/sum(phiM[i,,t])
      PSIM[i,4,t] <- phiM[i,4,t]/sum(phiM[i,,t])
      }
    }
    # Define observation matrix
    # Order of indices: true state, time, observed state, survey
    for (t in 1:T){
      beta1[t]~dnorm(mu.beta1, tau.beta1)
      beta2[t]~dnorm(mu.beta2, tau.beta2)
      beta3[t]~dnorm(mu.beta3, tau.beta3)
      beta4[t]~dnorm(mu.beta4, tau.beta4)
      beta5[t]~dnorm(mu.beta5, tau.beta5)
    for (j in 1:K){
    q[1,j,1,t] <- 1
    q[1,j,2,t] <- 0
    q[1,j,3,t] <- 0
    q[1,j,4,t] <- 0
    q[2,j,1,t] <- 1
    q[2,j,2,t] <- exp(beta1[t])
    q[2,j,3,t] <- 0
    q[2,j,4,t] <- 0
    q[3,j,1,t] <- 1
    q[3,j,2,t] <- 0
    q[3,j,3,t] <- exp(beta2[t])
    q[3,j,4,t] <- 0
    q[4,j,1,t] <- 1
    q[4,j,2,t] <- exp(beta3[t])
    q[4,j,3,t] <- exp(beta4[t])
    q[4,j,4,t] <- exp(beta3[t] + beta4[t] + beta5[t])
    }
      }
    # State-space likelihood
    # State equation: model of true states (z)
    for (t in 1:T){
    for (i in 1:N[t]){
     z[i,t] ~ dcat(phi[i,,t])
    }
    }
    # Observation equation
    for (t in 1:T){
    for (i in 1:N[t]){
       for (j in 1:K){ 
        y[i,j,t] ~ dcat(q[z[i,t],j,,t])
       } #j
    } #i
    }#t

#Derive outputs      
for (t in 1:T){      
      psiDay[t] <- PSI[1,2,t]  
      psiNight[t] <- PSI[1,3,t]
      psiND[t] <- PSI[1,4,t]
      #Derive detection parameters for occasion 1
      pDay[t] <- q[2,1,2,t]/sum(q[2,1,,t])
      pNight[t] <- q[3,1,3,t]/sum(q[3,1,,t])
      pND.ND[t] <- q[4,1,4,t]/sum(q[4,1,,t])
      pND.N[t] <- q[4,1,3,t]/sum(q[4,1,,t])
      pND.D[t] <- q[4,1,2,t]/sum(q[4,1,,t])
      pND.0[t] <- q[4,1,1,t]/sum(q[4,1,,t])      
}      
}