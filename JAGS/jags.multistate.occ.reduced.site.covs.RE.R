    model  { 
    # Priors for slopes
      alpha1.1~dlogis(0,1)
      alpha2.1~dlogis(0,1)
      #Priors for mean population level intercepts and effects
      mu.alpha.day~dnorm(0,0.5)
      mu.alpha.night~dnorm(0,0.5)
      mu.beta.day~dnorm(0,0.5)
      mu.beta.night~dnorm(0,0.5)
      #Priors for the variance of the population level effects
      s2.alpha.day~dgamma(1,1)
      s2.alpha.night~dgamma(1,1)
      s2.beta.day~dgamma(1,1)
      s2.beta.night~dgamma(1,1)
      #inverse of variance
      tau.alpha.day<-1/s2.alpha.day
      tau.alpha.night<-1/s2.alpha.night
      tau.beta.day<-1/s2.beta.day
      tau.beta.night<-1/s2.beta.night
  for (t in 1:T){
      #These the population-level random effects 
        alpha1[t]~dnorm(mu.alpha.day, tau.alpha.day)
        alpha2[t]~dnorm(mu.alpha.day, tau.alpha.night)
    # Define state vector
    for (i in 1:N[t]){
    phi[i,1,t] <- 1
    phi[i,2,t] <- exp(alpha1[t]+alpha1.1*cov1[i,t])
    phi[i,3,t] <- exp(alpha2[t]+alpha2.1*cov2[i,t])
    phi[i,4,t] <- exp(alpha1[t]+alpha1.1*cov1[i,t]+alpha2[t]+alpha2.1*cov2[i,t])
    PSI[i,1,t] <- phi[i,1,t]/sum(phi[i,,t])
    PSI[i,2,t] <- phi[i,2,t]/sum(phi[i,,t])
    PSI[i,3,t] <- phi[i,3,t]/sum(phi[i,,t])
    PSI[i,4,t] <- phi[i,4,t]/sum(phi[i,,t])
    }
}
  # Define observation matrix
    for (t in 1:T){
        beta1[t]~dnorm(mu.beta.day, tau.beta.day)
        beta2[t]~dnorm(mu.beta.night, tau.beta.night)
    # Order of indices: true state, time, observed state, survey
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
    q[4,j,2,t] <- exp(beta1[t])
    q[4,j,3,t] <- exp(beta2[t])
    q[4,j,4,t] <- exp(beta1[t]+beta2[t])
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
      
      #Mean population-level occurence
      # for(i in 1:N[t]){
      #   phiM[i,1] <- 1
      #   phiM[i,2] <- exp(mu.alpha1)
      #   phiM[i,3] <- exp(mu.alpha2)
      #   phiM[i,4] <- exp(mu.alpha1+mu.alpha2)
      #   PSIM[i,1] <- phiM[i,1]/sum(phiM[i,])
      #   PSIM[i,2] <- phiM[i,2]/sum(phiM[i,])
      #   PSIM[i,3] <- phiM[i,3]/sum(phiM[i,])
      #   PSIM[i,4] <- phiM[i,4]/sum(phiM[i,])
      # }
    }
}