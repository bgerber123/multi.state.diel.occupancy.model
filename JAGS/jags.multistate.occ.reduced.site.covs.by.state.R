
    model  { 
    
    # Priors
    pNight ~ dunif(0, 1)
    pDay ~ dunif(0, 1)
    
    #Define priors for alpha parameters on logit-scale
    for (i in 1:K.day) {
    alpha.day[i] ~ dlogis(0,1)
    }

    for (i in 1:K.night) {
    alpha.night[i] ~ dlogis(0,1)
    }


    # Define state vector
    for (s in 1:R){

    logit(aa[s]) <- inprod(Xday[s,],alpha.day)
    logit(bb[s])   <- inprod(Xnight[s,],alpha.night)

    prob[s,2] <- aa[s]*(1-bb[s])
    prob[s,3] <- bb[s]*(1-aa[s])
    prob[s,4] <- aa[s]*bb[s]
    prob[s,1] <- 1-prob[s,2]-prob[s,3]-prob[s,4]
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
    p[4,t,1] <- (1-pNight)*(1-pDay)
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
    
