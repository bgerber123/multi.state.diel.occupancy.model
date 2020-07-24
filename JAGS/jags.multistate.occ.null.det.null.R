    model  { 
    # Priors
    alpha ~ dlogis(0,1)
    pNight ~ dbeta(1, 1)  
    pDay ~ dbeta(1, 1)
    for (m in 1:4) {
      beta[m] ~ dgamma(1, 1)   # Induce Dirichlet prior
      pND[m] <- beta[m]/sum(beta[])
    }
    # Define state vector
    for (i in 1:N){
      phi[i,1] <- 3
      phi[i,2] <- exp(alpha)
      phi[i,3] <- exp(alpha)
      phi[i,4] <- exp(alpha)
      #Save the site-level occupancy probabilties      
      PSI[i,1] <- phi[i,1]/sum(phi[i,])
      PSI[i,2] <- phi[i,2]/sum(phi[i,])
      PSI[i,3] <- phi[i,3]/sum(phi[i,])
      PSI[i,4] <- phi[i,4]/sum(phi[i,])
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
   #The prob of occupancy derived for site 1 in state 2.   
   psi<-phi[1,2]/sum(phi[1,])
   #Overall occurence is the sum of states 2,3, and 4. Here, fixed for site 1
   psi.overall <- PSI[1,2]+PSI[1,3]+PSI[1,4]
}