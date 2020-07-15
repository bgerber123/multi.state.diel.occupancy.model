    model  { 
    # Priors
    pNight ~ dunif(0, 1)
    pDay ~ dunif(0, 1)
    psiNight ~ dunif(0, 1)
    psiDay ~ dunif(0, 1)
    # Define state vector
    for (s in 1:R){
    prob[s,1] <-  (1-psiDay)*(1-psiNight)                  #1-(psiDay-psiDay*psiNight)-(psiNight-psiDay*psiNight)-psiDay*psiNight 
    prob[s,2] <-  psiDay*(1-psiNight)  
    prob[s,3] <-  psiNight*(1-psiDay)
    prob[s,4] <-  psiDay*psiNight
    }
    # Define observation matrix
    # Order of indices: true state, time, observed state
    for (t in 1:T){
    p[1,t,1] <- 1
    p[1,t,2] <- 0
    p[1,t,3] <- 0
    p[1,t,4] <- 0
    p[2,t,1] <- 1-pDay
    p[2,t,2] <- pDay
    p[2,t,3] <- 0
    p[2,t,4] <- 0
    p[3,t,1] <- 1-pNight
    p[3,t,2] <- 0
    p[3,t,3] <- pNight
    p[3,t,4] <- 0
    p[4,t,1] <- (1-pDay)*(1-pNight)
    p[4,t,2] <- pDay*(1-pNight)
    p[4,t,3] <- pNight*(1-pDay)
    p[4,t,4] <- pDay*pNight 
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