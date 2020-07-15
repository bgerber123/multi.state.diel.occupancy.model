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
    for (i in 1:K) {
    alpha[i] ~ dlogis(0,1)
    }
    # Define state vector for each s site
    for (s in 1:R){
      phi[s,1] <- 1
      phi[s,2] <- exp(alpha[1]+alpha[2]*cov1[s])
      phi[s,3] <- exp(alpha[3]+alpha[4]*cov1[s])
      phi[s,4] <- exp(alpha[5]+alpha[6]*cov1[s])
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