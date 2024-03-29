    model  { 
    # Priors
    pNight ~ dbeta(1, 1)
    pDay ~ dbeta(1, 1)
    for(q in 1:Q){
      alpha[q]~dlogis(0,1)
    }
    # Define state vector
    for (i in 1:N){
    logit(psiDay[i]) <- alpha[1]+alpha[2]*x[i]
    logit(psiNight[i]) <- alpha[3]+alpha[4]*x[i]
    PSI[i,1] <- (1-psiNight[i])*(1-psiDay[i])
    PSI[i,2] <- psiDay[i]*(1-psiNight[i])
    PSI[i,3] <- psiNight[i]*(1-psiDay[i])
    PSI[i,4] <- psiDay[i]*psiNight[i]
    }
    # Define observation matrix
    # Order of indices: true state, time, observed state
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
    p[4,j,1] <- (1-pNight)*(1-pDay)
    p[4,j,2] <- pDay*(1-pNight)
    p[4,j,3] <- pNight*(1-pDay)
    p[4,j,4] <- pDay*pNight 
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