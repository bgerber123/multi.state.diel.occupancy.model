
    model  { 
    
    # Priors
    pNight ~ dunif(0, 1)  
    pDay ~ dunif(0, 1)


    for (i in 1:4) {
      beta[i] ~ dgamma(1, 1)   # Induce Dirichlet prior
      pND[i] <- beta[i]/sum(beta[])

      alpha[i] ~ dgamma(1, 1)   # Induce Dirichlet prior
      psi[i] <- alpha[i]/sum(alpha[])
    }

    
    
    # Define state vector
    for (s in 1:R){
      prob[s,1] <- psi[1]
      prob[s,2] <- psi[2]
      prob[s,3] <- psi[3]
      prob[s,4] <- psi[4]
    }
    
    # Define observation matrix
    # Order of indices: true state, time, observed state
    for (t in 1:T){
    p[1,t,1] <- 1
    p[1,t,2] <- 0
    p[1,t,3] <- 0
    p[1,t,4] <- 0
    p[2,t,1] <- 1-pDay
    p[2,t,2] <- pDay
    p[2,t,3] <- 0
    p[2,t,4] <- 0
    p[3,t,1] <- 1-pNight
    p[3,t,2] <- 0
    p[3,t,3] <- pNight
    p[3,t,4] <- 0
    p[4,t,1] <- pND[1]
    p[4,t,2] <- pND[2]
    p[4,t,3] <- pND[3]
    p[4,t,4] <- pND[4]
    }
    
    # State-space likelihood
    # State equation: model of true states (z)
    for (s in 1:R){
     z[s] ~ dcat(prob[s,])
    }
    
    # Observation equation
    for (s in 1:R){
       for (t in 1:T){ 
        y[s,t] ~ dcat(p[z[s],t,])
       } #t
    } #s
    

    }
    
