    model  { 
    # Priors for day/night detection
    pNight ~ dbeta(1, 1)
    pDay ~ dbeta(1, 1)
    # Priors for day&night detection. 
    #state i to sum to 1.
    for (m in 1:4) {
      beta[m] ~ dgamma(1, 1)   # Induce Dirichlet prior
      pND[m] <- beta[m]/sum(beta[])
    }
    #Define priors for alpha parameters on logit-scale
    for (q in 1:Q.day) {
    alpha.day[q] ~ dlogis(0,1)
    }
    for (q in 1:Q.night) {
    alpha.night[q] ~ dlogis(0,1)
    }
    for (q in 1:Q.nd) {
    alpha.nd[q] ~ dlogis(0,1)
    }
    # Define state vector for each s site
    for (i in 1:N){
      phi[i,1] <- 1
      phi[i,2] <- exp(inprod(Xday[i,],alpha.day))
      phi[i,3] <- exp(inprod(Xnight[i,],alpha.night))
      phi[i,4] <- exp(inprod(Xday[i,],alpha.day)+inprod(Xnight[i,],alpha.night)+inprod(Xnd[i,],alpha.nd))
      PSI[i,1] <- phi[i,1]/sum(phi[i,])
      PSI[i,2] <- phi[i,2]/sum(phi[i,])
      PSI[i,3] <- phi[i,3]/sum(phi[i,])
      PSI[i,4] <- phi[i,4]/sum(phi[i,])
    }
    # Define observation matrix
    # Order of indices: true state, time, observed state
    for (j in 1:K){
    p[1,j,1] <- 1
    p[1,j,2] <- 0
    p[1,j,3] <- 0
    p[1,j,4] <- 0
    p[2,j,1] <- 1-pNight
    p[2,j,2] <- pNight
    p[2,j,3] <- 0
    p[2,j,4] <- 0
    p[3,j,1] <- 1-pDay
    p[3,j,2] <- 0
    p[3,j,3] <- pDay
    p[3,j,4] <- 0
    p[4,j,1] <- pND[1]
    p[4,j,2] <- pND[2]
    p[4,j,3] <- pND[3]
    p[4,j,4] <- pND[4]
    }
    # State-space likelihood
    # State equation: model of true states (z)
    for (i in 1:N){
     z[i] ~ dcat(PSI[i,])
    }
    # Observation equation
    for (i in 1:N){
       for (j in 1:K){ 
        y[i,j] ~ dcat(p[z[i],j,])
       } #j
    } #i
}