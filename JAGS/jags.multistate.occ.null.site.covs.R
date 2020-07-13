    model  { 
    
    # Priors
    for(i in 1:K){
      alpha[i]~dlogis(0,1)
    }
    beta ~ dlogis(0,1)


 #Define state vector for each s site
    for (s in 1:R){
      log(phi[s,1]) <- log(3)
      log(phi[s,2]) <- alpha[1]+alpha[2]*cov1[s]
      log(phi[s,3]) <- alpha[1]+alpha[2]*cov1[s]
      log(phi[s,4]) <- alpha[1]+alpha[2]*cov1[s]

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
    p[2,t,1] <- 1-pdet
    p[2,t,2] <- pdet
    p[2,t,3] <- 0
    p[2,t,4] <- 0
    p[3,t,1] <- 1-pdet
    p[3,t,2] <- 0
    p[3,t,3] <- pdet
    p[3,t,4] <- 0

    log(q[4,t,1]) <- log(3)
    log(q[4,t,2]) <- beta
    log(q[4,t,3]) <- beta
    log(q[4,t,4]) <- beta

    p[4,t,1] <- q[4,t,1]/sum(q[4,t,])
    p[4,t,2] <- q[4,t,2]/sum(q[4,t,])
    p[4,t,3] <- q[4,t,3]/sum(q[4,t,])
    p[4,t,4] <- q[4,t,4]/sum(q[4,t,])
 
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

   #The prob of occupancy derived for site 1 in state 2.   
   psi<-phi[1,2]/sum(phi[1,])

   #Overall occurence is the sum of states 2,3, and 4. Here, fixed for site 1
   psi.overall <- phi[1,2]/sum(phi[1,])+phi[1,3]/sum(phi[1,])+phi[1,4]/sum(phi[1,])

   #The probability of detection is equal to the sum of detection in state 2,3, and 4.
   #These estimates come from the beta parameter. Here, I have fixed this to be occasion 1 (t=1).
   pdet <- p[4,1,2]+p[4,1,3]+p[4,1,4]


    }
    
