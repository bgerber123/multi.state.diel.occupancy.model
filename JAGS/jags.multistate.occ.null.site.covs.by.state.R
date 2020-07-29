    model  { 
    # Priors
    for(q in 1:Q){
      alpha[q]~dlogis(0,1)
    }
    beta ~ dlogis(0,1)
 #Define state vector for each i site
    for (i in 1:N){
      phi[i,1] <- 3
      phi[i,2] <- exp(inprod(X[i,],alpha))
      phi[i,3] <- exp(inprod(X[i,],alpha))
      phi[i,4] <- exp(inprod(X[i,],alpha))
      #Save the site-level occupancy probabilties      
      PSI[i,1]  <- phi[i,1]/sum(phi[i,])
      PSI[i,2]  <- phi[i,2]/sum(phi[i,])
      PSI[i,3]  <- phi[i,3]/sum(phi[i,])
      PSI[i,4]  <- phi[i,4]/sum(phi[i,])
    }
    # Define observation matrix
    # Order of indices: true state, time, observed state
    for (j in 1:K){
    q[1,j,1] <- 1
    q[1,j,2] <- 0
    q[1,j,3] <- 0
    q[1,j,4] <- 0
    q[2,j,1] <- 1
    q[2,j,2] <- exp(beta)
    q[2,j,3] <- 0
    q[2,j,4] <- 0
    q[3,j,1] <- 1
    q[3,j,2] <- 0
    q[3,j,3] <- exp(beta)
    q[3,j,4] <- 0
    q[4,j,1] <- 3
    q[4,j,2] <- exp(beta)
    q[4,j,3] <- exp(beta)
    q[4,j,4] <- exp(beta)
    }
    # State-space likelihood
    # State equation: model of true states (z)
    for (i in 1:N){
    z[i] ~ dcat(PSI[i,])
    }
    # Observation equation
    for (i in 1:N){
    for (j in 1:K){ 
    y[i,j] ~ dcat(q[z[i],j,])
    } #j
    } #i
    #The prob of occupancy derived for site 1 in state 2.   
    psi<-phi[1,2]/sum(phi[1,])
    #Overall occurence is the sum of states 2,3, and 4. Here, fixed for site 1
    psi.overall <- PSI[1,2]+PSI[1,3]+PSI[1,4]
   #The probability of detection is equal to the sum of detection in state 2,3, and 4.
   #These estimates come from the beta parameter. Here, I have fixed this to be occasion 1 (t=1).
  p.overall <- q[2,1,2]/sum(q[2,1,])
}
    
