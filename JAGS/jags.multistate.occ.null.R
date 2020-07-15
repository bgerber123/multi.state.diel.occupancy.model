    model  { 
    # Priors
    alpha ~ dlogis(0,1)
    beta ~ dlogis(0,1)
 #Define state vector for each s site
    for (s in 1:R){
      phi[s,1] <- 3
      phi[s,2] <- exp(alpha)
      phi[s,3] <- exp(alpha)
      phi[s,4] <- exp(alpha)
    }
    # Define observation matrix
    # Order of indices: true state, time, observed state
    for (t in 1:T){
    p[1,t,1] <- 1
    p[1,t,2] <- 0
    p[1,t,3] <- 0
    p[1,t,4] <- 0
    p[2,t,1] <- 1
    p[2,t,2] <- exp(beta)
    p[2,t,3] <- 0
    p[2,t,4] <- 0
    p[3,t,1] <- 1
    p[3,t,2] <- 0
    p[3,t,3] <- exp(beta)
    p[3,t,4] <- 0
    p[4,t,1] <- 3
    p[4,t,2] <- exp(beta)
    p[4,t,3] <- exp(beta)
    p[4,t,4] <- exp(beta)
    }
    # State-space likelihood
    # State equation: model of true states (z)
    for (s in 1:R){
    z[s] ~ dcat(phi[s,])
    }
    # Observation equation
    for (s in 1:R){
    for (t in 1:T){ 
    y[s,t] ~ dcat(p[z[s],t,])
    } #t
    } #s
   #The prob of occupancy derived for site 1 in state 2.   
   psi<-phi[1,2]/sum(phi[1,])
   #Overall occurence is the sum of states 2,3, and 4. Here, fixed for site 1
   psi.overall <- phi[1,2]/sum(phi[1,])+phi[1,3]/sum(phi[1,])+phi[1,4]/sum(phi[1,])
   #The probability of detection is equal to the sum of detection in state 2,3, and 4.
   #These estimates come from the beta parameter. Here, I have fixed this to be occasion 1 (t=1).
   pdet <- p[4,1,2]/sum(p[4,1,])+p[4,1,3]/sum(p[4,1,])+p[4,1,4]/sum(p[4,1,])
    }
    
