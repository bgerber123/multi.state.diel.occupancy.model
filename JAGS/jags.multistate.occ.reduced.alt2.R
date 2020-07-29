    model  { 
    # Priors
    beta1 ~ dlogis(0, 1)
    beta2 ~ dlogis(0, 1)
    psiNight ~ dbeta(1, 1)
    psiDay ~ dbeta(1, 1)
    # Define state vector
    for (i in 1:N){
    PSI[i,1] <-  (1-psiDay)*(1-psiNight)                  #1-(psiDay-psiDay*psiNight)-(psiNight-psiDay*psiNight)-psiDay*psiNight 
    PSI[i,2] <-  psiDay*(1-psiNight)  
    PSI[i,3] <-  psiNight*(1-psiDay)
    PSI[i,4] <-  psiDay*psiNight
    }
    # Define observation matrix
    # Order of indices: true state, survey occ, observed state
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
    q[4,j,4] <- exp(beta1+beta2)
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
    pDay <- q[2,1,2]/sum(q[2,1,])
    pNight <- q[3,1,3]/sum(q[3,1,])
    pNDND <- q[4,1,4]/sum(q[4,1,])
}