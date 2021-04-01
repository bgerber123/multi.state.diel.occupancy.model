model{
######
# Detection model
######
for(i in 1:nsite) {
  for(ti in 1:nseason) {
    for(survey in 1:nsurvey) {
      y[i, ti, survey] ~ dcat( rdm[i,  , z[i, ti], ti ] )
    }
  }
}
#######
# Latent state model
#######
for(i in 1:nsite) {
  # for first season
  z[i, 1] ~ dcat( PSI[i, ] )
  for(t in 2:nseason) {
    z[i, t] ~ dcat( tpm[i,  , z[ i, t-1], t-1] )
  }
}
for( i in 1:nsite ) {
######
# Fill in all of the transition probabilities
######
######
# Latent state for first season, fsm = transition vector for season 1
######
# First season probabilities for each state
PSI[i, 1] <- 1 #--------------------------------------------------|U
PSI[i, 2] <- exp( psi[1, i] ) #-----------------------------------|D
PSI[i, 3] <- exp( psi[2, i] ) #-----------------------------------|N
PSI[i, 4] <- exp( psi_one[i] ) #----------------------------------|DN
for(t in 2:nseason) {
######
# Latent state for dynamic part of model
# tpm = transition probability matrix. All columns sum to 1.
# dim(tpm)[1] = site j 
# dim(tpm)[2] = state at time t
# dim(tpm)[3] = state at time t-1
# dim(tpm)[4] = time t (for temporally varying covariates)
######
# U to ...
tpm[i, 1, 1, t-1] <- 1 #----------------------------------------------------|U
tpm[i, 2, 1, t-1] <- exp( gam[1, i, t-1] ) #--------------------------------|D
tpm[i, 3, 1, t-1] <- exp(                  gam[2, i, t-1] ) #---------------|N
tpm[i, 4, 1, t-1] <- exp( gam[1, i, t-1] + gam[2, i, t-1] ) #---------------|DN
# D to ...
tpm[i, 1, 2, t-1] <- exp( eps[1, i, t-1] ) #--------------------------------|U
tpm[i, 2, 2, t-1] <- 1 #----------------------------------------------------|D
tpm[i, 3, 2, t-1] <- exp( eps[1, i, t-1] + gam_one[2, 1, i, t-1] ) #--------|N
tpm[i, 4, 2, t-1] <- exp(                  gam_one[2, 1, i, t-1] ) #--------|DN
# N to ...
tpm[i, 1, 3, t-1] <- exp(                         eps[2, i, t-1] ) #--------|U
tpm[i, 2, 3, t-1] <- exp( gam_one[1, 2, i, t-1] + eps[2, i, t-1] ) #--------|D
tpm[i, 3, 3, t-1] <- 1 #----------------------------------------------------|N
tpm[i, 4, 3, t-1] <- exp( gam_one[1, 2, i, t-1] ) #-------------------------|DN
# DN to ..
tpm[i, 1, 4, t-1] <- exp( eps_one[1, 2, i, t-1] + eps_one[2, 1, i, t-1] ) #-|U
tpm[i, 2, 4, t-1] <- exp( eps_one[2, 1, i, t-1] ) #-------------------------|D
tpm[i, 3, 4, t-1] <- exp(                         eps_one[1, 2, i, t-1] ) #-|N
tpm[i, 4, 4, t-1] <- 1 #----------------------------------------------------|DN
} # close t year loop
######
# detection matrix (OS = observed state, TS = true state)
# rdm = rho detection matrix. Each row sums to 1.
# OS along rows, TS along columns
######
# TS = U
for(ti in 1:nseason) {
rdm[i, 1, 1, ti] <- 1 #-------------------------------------------------|OS = U
rdm[i, 2, 1, ti] <- 0 #-------------------------------------------------|OS = D
rdm[i, 3, 1, ti] <- 0 #-------------------------------------------------|OS = N
rdm[i, 4, 1, ti] <- 0 #-------------------------------------------------|OS = DN
# TS = D
rdm[i, 1, 2, ti] <- 1 #-------------------------------------------------|OS = U
rdm[i, 2, 2, ti] <- exp( rho[1, i, ti] ) #------------------------------|OS = D
rdm[i, 3, 2, ti] <- 0 #-------------------------------------------------|OS = N
rdm[i, 4, 2, ti] <- 0 #-------------------------------------------------|OS = DN
# TS = N
rdm[i, 1, 3, ti] <- 1 #-------------------------------------------------|OS = U
rdm[i, 2, 3, ti] <- 0 #-------------------------------------------------|OS = D
rdm[i, 3, 3, ti] <- exp( rho[2, i, ti] ) #------------------------------|OS = N
rdm[i, 4, 3, ti] <- 0 #-------------------------------------------------|OS = DN
# TS = DN
rdm[i, 1, 4, ti] <- 1 #-------------------------------------------------|OS = U
rdm[i, 2, 4, ti] <- exp( rho[1+to_add, i, ti] ) #-----------------------|OS = D
rdm[i, 3, 4, ti] <- exp(                  rho[2+to_add, i, ti] ) #------|OS = N
rdm[i, 4, 4, ti] <- exp( rho[1+to_add, i, ti] +  rho[2+to_add, i, ti] )#|OS = DN
} # close ti year
######
# Fill in the linear predictors for the transition matrices
######
for( s in 1:ncat ) {
# base occupancy
psi[s ,i] <- inprod( a[s, ], psi_cov[i, ] )
}
for(t in 2:nseason) {
 for(s in 1:ncat){
# base colonization
gam[s, i, t-1] <- inprod( b[s, ], gam_cov[i, ,t-1] ) 
# base extinction
eps[s, i, t-1] <- inprod( d[s, ], eps_cov[i, , t-1] )
# species interactions
 pi[s, s, i, t-1] <- 0 # fill inxs diagonal with 0's
tau[s, s, i, t-1] <- 0 # fill inxs diagonal with 0's
# These are the matrices to hold conditional linear predictors for gamma and eps
#  when one other species is present at t-1 (hence the _one after the name)
gam_one[s, s, i, t-1] <- 0 # fill inxs diagonal with 0's
eps_one[s, s, i, t-1] <- 0 # fill inxs diagonal with 0's
}
# 
# In the model we put the inxs in a matrix
#  these matrices get filled in this order:
#######
#0|4|5#
#-|-|-#
#1|0|6# Zeros along the diagonal
#-|-|-#
#2|3|0#
#######
# We have set it up in the model such that the species in the off-diagonal
#  influence the species along the diagonal. For example, pi[3,2] is the 
#  influence species 2 has on the colonization rate of species 3 while 
#  pi[2, 1] is the influence species 1 has on species 2. More generally it
#  is pi[species of interest, species conditioning on]. 
  for(k in 1:n_inxs){
    # inxs on colonization 
     pi[rows_vec[k], cols_vec[k], i, t-1] <- inprod( g[k, ], pi_cov[i, ,t-1] )
    # inxs on extinction 
    tau[rows_vec[k], cols_vec[k], i, t-1] <- inprod( h[k, ], tau_cov[i, ,t-1] )
    # linear predictor for colonization | one species present at t-1
    gam_one[rows_vec[k], cols_vec[k], i, t-1] <- 
      inprod( b[rows_vec[k], ], gam_cov[i, ,t-1] ) + pi[rows_vec[k], cols_vec[k], i, t-1]
    # linear predictor for extinction | one species present at t-1
    eps_one[rows_vec[k], cols_vec[k], i, t-1] <- 
      inprod( d[rows_vec[k], ], eps_cov[i, ,t-1] ) + tau[rows_vec[k], cols_vec[k], i, t-1]
      } # close inxs
    } # close t
for(srho in 1:(ncat+to_add)){
  for(ti in 1:nseason){
    # base detection probability
    rho[srho, i, ti] <- inprod( f[srho, ], rho_cov[i, , ti] ) 
  }
} # close category
psi_one[i] <- psi[1,i] + psi[2,i] + inprod(a_inxs, psi_inxs_cov[i,] )
} # closes for loop for i (sites) all the way up at the top of the model
#####
# Priors
######
  lambda ~ dunif(0.001, 20)
  lambda_det ~ dunif(0.001, 20)
  for(psi_inxsp in 1:ncov_psi_inxs){
    a_inxs0[psi_inxsp] ~ dlogis(0, 1)
    a_inxs1[psi_inxsp] ~ ddexp(0, lambda)
    a_inxsind[psi_inxsp] <- ifelse(psi_inxsp == 1, 1, 0)
    a_inxs[psi_inxsp] <- (a_inxsind[psi_inxsp] * a_inxs0[psi_inxsp]) +
      ((1-a_inxsind[psi_inxsp]) * a_inxs1[psi_inxsp])
  }
  for(s in 1:ncat){
    # Initial Occupancy
    for( psip in 1:ncov_psi ){
      a0[s,psip] ~ dlogis(0, 1)
      a1[s, psip] ~ ddexp(0, lambda)
      aind[s, psip] <- ifelse(psip == 1, 1, 0)
      a[s, psip] <- (aind[s,psip] * a0[s,psip]) +
        ((1-aind[s,psip]) * a1[s,psip])
    }
    # Colonization
    for( gamp in 1:ncov_gam ){
      b0[s,gamp] ~ dlogis(0, 1)
      b1[s, gamp] ~ ddexp(0, lambda)
      bind[s, gamp] <- ifelse(gamp == 1, 1, 0)
      b[s, gamp] <- (bind[s,gamp] * b0[s,gamp]) +
        ((1-bind[s,gamp]) * b1[s,gamp])
    }  
    # Extinction
    for( epsp in 1:ncov_eps ){
      d0[s,epsp] ~ dlogis(0, 1)
      d1[s, epsp] ~ ddexp(0, lambda)
      dind[s, epsp] <- ifelse(epsp == 1, 1, 0)
      d[s, epsp] <- (dind[s,epsp] * d0[s,epsp]) +
        ((1-dind[s,epsp]) * d1[s,epsp])
    }
  }
  # Detection
  for(srho in 1:(ncat+to_add)){
    for( rhop in 1:ncov_rho ){
      f0[srho, rhop] ~ dlogis(0, 1)
      f1[srho, rhop] ~ ddexp(0, lambda_det)
      find[srho, rhop] <- ifelse(rhop == 1, 1, 0)
      f[srho, rhop] <- (find[srho, rhop] * f0[srho, rhop]) +
        ((1 - find[srho, rhop]) * f1[srho, rhop])
    }
  }
  for( k in 1:n_inxs ){
    # Inxs on colonization
    for( pip in 1:ncov_pi ){
      g0[k, pip] ~ dlogis(0, 1)
      g1[k, pip] ~ ddexp(0, lambda)
      gind[k, pip] <- ifelse(pip == 1, 1, 0)
      g[k, pip] <- (gind[k, pip] * g0[k, pip]) + 
        ((1 - gind[k, pip]) * g1[k, pip])
    }
    # Inxs on extinction
    for( taup in 1:ncov_tau ){
      h0[k, taup] ~ dlogis(0, 1)
      h1[k, taup] ~ ddexp(0, lambda)
      hind[k, taup] <- ifelse(taup == 1, 1, 0)
      h[k, taup] <- (hind[k, taup] * h0[k, taup]) + 
        ((1 - hind[k, taup]) * h1[k, taup])
    }
  }
}