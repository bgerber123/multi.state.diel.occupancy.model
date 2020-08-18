    model  { 
    # Priors
      alpha1.1~dlogis(0,1)
      #Priors for mean population level intercepts and effects
      mu.alpha~dnorm(0,0.5)
      mu.beta~dnorm(0,0.5)
      #Priors for the variance of the population level effects
      s2.alpha~dgamma(1,1)
      s2.beta~dgamma(1,1)
      #inverse of variance
      tau.alpha<-1/s2.alpha
      tau.beta<-1/s2.beta
    for (t in 1:T){
        #These the population-level random effects 
        alpha[t]~dnorm(mu.alpha, tau.alpha)
        # Define state vector
    for (i in 1:N[t]){
      phi[i,1,t] <- 3
      phi[i,2,t] <- exp(alpha[t]+alpha1.1*cov[i,t])
      phi[i,3,t] <- exp(alpha[t]+alpha1.1*cov[i,t])
      phi[i,4,t] <- exp(alpha[t]+alpha1.1*cov[i,t])
      #Save the site-level occupancy probabilties      
      PSI[i,1,t] <- phi[i,1,t]/sum(phi[i,,t])
      PSI[i,2,t] <- phi[i,2,t]/sum(phi[i,,t])
      PSI[i,3,t] <- phi[i,3,t]/sum(phi[i,,t])
      PSI[i,4,t] <- phi[i,4,t]/sum(phi[i,,t])
    }
  }
    # Define observation matrix
      for (t in 1:T){
        beta[t]~dnorm(mu.beta, tau.beta)
        # Order of indices: true state, occasion, observed state, survey
    for (j in 1:K){
    q[1,j,1,t] <- 1
    q[1,j,2,t] <- 0
    q[1,j,3,t] <- 0
    q[1,j,4,t] <- 0
    q[2,j,1,t] <- 1
    q[2,j,2,t] <- exp(beta[t])
    q[2,j,3,t] <- 0
    q[2,j,4,t] <- 0
    q[3,j,1,t] <- 1
    q[3,j,2,t] <- 0
    q[3,j,3,t] <- exp(beta[t])
    q[3,j,4,t] <- 0
    q[4,j,1,t] <- 3
    q[4,j,2,t] <- exp(beta[t])
    q[4,j,3,t] <- exp(beta[t])
    q[4,j,4,t] <- exp(beta[t])
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
   #The prob of occupancy derived for site 1 in state 2.   
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
      
      # #Mean population-level occurence
      # for(i in 1:N[t]){
      #   phiM[i,1] <- 1
      #   phiM[i,2] <- exp(mu.alpha)
      #   phiM[i,3] <- exp(mu.alpha)
      #   phiM[i,4] <- exp(mu.alpha)
      #   PSIM[i,1] <- phiM[i,1]/sum(phiM[i,])
      #   PSIM[i,2] <- phiM[i,2]/sum(phiM[i,])
      #   PSIM[i,3] <- phiM[i,3]/sum(phiM[i,])
      #   PSIM[i,4] <- phiM[i,4]/sum(phiM[i,])
      # }
      }
    }