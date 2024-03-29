    model  { 
    # Priors
    pNight ~ dbeta(1, 1)  
    pDay ~ dbeta(1, 1)
    for (m in 1:4) {
      beta[m] ~ dgamma(1, 1)   # Induce Dirichlet prior
      pND[m] <- beta[m]/sum(beta[])

      alpha[m] ~ dgamma(1, 1)   # Induce Dirichlet prior
      psi[m] <- alpha[m]/sum(alpha[])
    }
    # Define state vector
    for (i in 1:N){
      PSI[i,1] <- psi[1]
      PSI[i,2] <- psi[2]
      PSI[i,3] <- psi[3]
      PSI[i,4] <- psi[4]
    }
    # Define observation matrix
    # Order of indices: true state, survey occ, observed state
    for (j in 1:K){
    p[1,j,1] <- 1
    p[1,j,2] <- 0
    p[1,j,3] <- 0
    p[1,j,4] <- 0
    p[2,j,1] <- 1-pDay
    p[2,j,2] <- pDay
    p[2,j,3] <- 0
    p[2,j,4] <- 0
    p[3,j,1] <- 1-pNight
    p[3,j,2] <- 0
    p[3,j,3] <- pNight
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