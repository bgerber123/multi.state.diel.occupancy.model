
    model  { 
    
    # Priors for day/night detection
    pNight ~ dunif(0, 1)
    pDay ~ dunif(0, 1)

    # Priors for day&night detection. Dirchlet prior forces pND for each 
    #state i to sum to 1.
    for (i in 1:4) {
      beta[i] ~ dgamma(1, 1)   # Induce Dirichlet prior
      pND[i] <- beta[i]/sum(beta[])
    }

    #Define priors for alpha parameters on logit-scale
    for (i in 1:K.day) {
    alpha.day[i] ~ dlogis(0,1)
    }

    for (i in 1:K.night) {
    alpha.night[i] ~ dlogis(0,1)
    }

    for (i in 1:K.nd) {
    alpha.nd[i] ~ dlogis(0,1)
    }

   
    
    # Define state vector for each s site
    for (s in 1:R){
      log(phi[s,1]) <- 1
      log(phi[s,2]) <- inprod(Xday[s,],alpha.day)
      log(phi[s,3]) <- inprod(Xnight[s,],alpha.night)
      log(phi[s,4]) <- inprod(Xnd[s,],alpha.nd)

      prob[s,1]     <- phi[s,1]/sum(phi[s,])
      prob[s,2]     <- phi[s,2]/sum(phi[s,])
      prob[s,3]     <- phi[s,3]/sum(phi[s,])
      prob[s,4]     <- phi[s,4]/sum(phi[s,])

    }
    
    # Define observation matrix
    # Order of indices: true state, time, observed state
    for (t in 1:T){
    p[1,t,1] <- 1
    p[1,t,2] <- 0
    p[1,t,3] <- 0
    p[1,t,4] <- 0
    p[2,t,1] <- 1-pNight
    p[2,t,2] <- pNight
    p[2,t,3] <- 0
    p[2,t,4] <- 0
    p[3,t,1] <- 1-pDay
    p[3,t,2] <- 0
    p[3,t,3] <- pDay
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
    