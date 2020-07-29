    model  { 
    # Priors
    psi.overall ~ dbeta(1,1)
    p.overall ~ dbeta(1,1)
 #Define state vector for each i site
    for (i in 1:N){
      PSI[i,1] <- 1-psi.overall
      PSI[i,2] <- psi.overall/3
      PSI[i,3] <- psi.overall/3
      PSI[i,4] <- psi.overall/3 
    }
    # Define observation matrix
    # Order of indices: true state, time, observed state
    for (j in 1:K){
    p[1,j,1] <- 1-p.overall
    p[1,j,2] <- 0
    p[1,j,3] <- 0
    p[1,j,4] <- 0
    p[2,j,1] <- 1-p.overall
    p[2,j,2] <- p.overall
    p[2,j,3] <- 0
    p[2,j,4] <- 0
    p[3,j,1] <- 1-p.overall
    p[3,j,2] <- 0
    p[3,j,3] <- p.overall
    p[3,j,4] <- 0
    p[4,j,2] <- p.overall/3
    p[4,j,3] <- p.overall/3
    p[4,j,4] <- p.overall/3
    p[4,j,1] <- 1-p.overall
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
   psi<-psi.overall/3 
    }