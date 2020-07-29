    model  { 
    # Priors
    beta1 ~ dlogis(0,1)
    beta2 ~ dlogis(0,1)
    beta3 ~ dlogis(0,1)
    alpha1 ~ dlogis(0,1)
    alpha2 ~ dlogis(0,1)
    alpha3 ~ dlogis(0,1)
    # Define state vector
    for (i in 1:N){
      phi[i,1] <- 1
      phi[i,2] <- exp(alpha1)
      phi[i,3] <- exp(alpha2)
      phi[i,4] <- exp(alpha1+alpha2+alpha3)
      PSI[i,1] <- phi[i,1]/phi[i,]
      PSI[i,2] <- phi[i,2]/phi[i,]
      PSI[i,3] <- phi[i,3]/phi[i,]
      PSI[i,4] <- phi[i,4]/phi[i,]
    }
    # Define observation matrix
    # Order of indices: true state, time, observed state
    for (j in 1:K){
    q[1,j,1] <- 1
    q[1,j,2] <- 0
    q[1,j,3] <- 0
    q[1,j,4] <- 0
    q[2,j,1] <- 1
    q[2,j,2] <- exp(beta1)
    q[2,j,3] <- 0
    q[2,j,4] <- 0
    q[3,j,1] <- 1
    q[3,j,2] <- 0
    q[3,j,3] <- exp(beta2)
    q[3,j,4] <- 0
    q[4,j,1] <- 1
    q[4,j,2] <- exp(beta1)
    q[4,j,3] <- exp(beta2)
    q[4,j,4] <- exp(beta1 + beta2 + beta3)
    }
    # State-space likelihood
    # State equation: model of true states (z)
    for (i in 1:N){
     z[i] ~ dcat(phi[i,])
    }
    # Observation equation
    for (i in 1:N){
       for (j in 1:K){ 
        y[i,j] ~ dcat(q[z[i],j,])
       } #j
    } #i
    #State Occupancy Parametsr for site 1
    psiDay <- PSI[1,2]
    psiNight <- PSI[1,3]
    psiND <- PSI[1,4]
    #Derive detection parameters for occasion 1
    pDay <- q[2,1,2]/sum(q[2,1,])
    pNight <- q[3,1,3]/sum(q[3,1,])
    pND.ND <- q[4,1,4]/sum(q[4,1,])
}