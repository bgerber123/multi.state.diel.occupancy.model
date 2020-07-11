model{
######
# Detection model
######
for(j in 1:nsite) {
  for(ti in 1:nyear) {
    for(week in 1:nsurvey) {
      y[j, ti, week] ~ dcat( rdm[j, ( 1:nout ) , z[j, ti], ti ] )
    }
  }
}
#######
# Latent state model
#######
for(j in 1:nsite) {
  # for first season
  z[j, 1] ~ dcat( fsm[j, ( 1:nout )] )
  for(t in 2:nyear) {
    z[j, t] ~ dcat( tpm[j, ( 1:nout ) , z[ j, t-1], t-1] )
  }
}
for( j in 1:nsite ) {
######
# Fill in all of the transition probabilities
######
######
# Latent state for first season, fsm = transition vector for season 1
######
# First season probabilities for each state
fsm[j, 1] <- 1 #------------------------------------------------------------|U
fsm[j, 2] <- exp( psi[1, j] ) #------------------------------------------|A
fsm[j, 3] <- exp( psi[2, j] ) #------------------------------------------|B
fsm[j, 4] <- exp( psi[1, j] + psi[2, j] ) #---------------------------|AB
for(t in 2:nyear) {
######
# Latent state for dynamic part of model
# tpm = transition probability matrix. All columns sum to 1.
# dim(tpm)[1] = site j 
# dim(tpm)[2] = state at time t
# dim(tpm)[3] = state at time t-1
# dim(tpm)[4] = time t (for temporally varying covariates)
######
# U to ...
tpm[j, 1, 1, t-1] <- 1 #----------------------------------------------------|U
tpm[j, 2, 1, t-1] <- exp( gam[1, j, t-1] ) #--------------------------------|A
tpm[j, 3, 1, t-1] <- exp(                  gam[2, j, t-1] ) #---------------|B
tpm[j, 4, 1, t-1] <- exp( gam[1, j, t-1] + gam[2, j, t-1] ) #---------------|AB
# A to ...
tpm[j, 1, 2, t-1] <- exp( eps[1, j, t-1] ) #--------------------------------|U
tpm[j, 2, 2, t-1] <- 1 #----------------------------------------------------|A
tpm[j, 3, 2, t-1] <- exp( eps[1, j, t-1] + gam_one[2, 1, j, t-1] ) #--------|B
tpm[j, 4, 2, t-1] <- exp(                  gam_one[2, 1, j, t-1] ) #--------|AB
# B to ...
tpm[j, 1, 3, t-1] <- exp(                         eps[2, j, t-1] ) #--------|U
tpm[j, 2, 3, t-1] <- exp( gam_one[1, 2, j, t-1] + eps[2, j, t-1] ) #--------|A
tpm[j, 3, 3, t-1] <- 1 #----------------------------------------------------|B
tpm[j, 4, 3, t-1] <- exp( gam_one[1, 2, j, t-1] ) #-------------------------|AB
# AB to ..
tpm[j, 1, 4, t-1] <- exp( eps_one[1, 2, j, t-1] + eps_one[2, 1, j, t-1] ) #-|U
tpm[j, 2, 4, t-1] <- exp( eps_one[1, 2, j, t-1] ) #-------------------------|A
tpm[j, 3, 4, t-1] <- exp(                         eps_one[2, 1, j, t-1] ) #-|B
tpm[j, 4, 4, t-1] <- 1 #----------------------------------------------------|AB
} # close t year loop
######
# detection matrix (OS = observed state, TS = true state)
# rdm = rho detection matrix. Each row sums to 1.
# OS along rows, TS along columns
######
# TS = U
for(ti in 1:nyear) {
rdm[j, 1, 1, ti] <- 1 #--------------------------------------------------|OS = U
rdm[j, 2, 1, ti] <- 0 #--------------------------------------------------|OS = A
rdm[j, 3, 1, ti] <- 0 #--------------------------------------------------|OS = B
rdm[j, 5, 1, ti] <- 0 #--------------------------------------------------|OS = AB
# TS = A
rdm[j, 1, 2, ti] <- 1 #--------------------------------------------------|OS = U
rdm[j, 2, 2, ti] <- exp( rho[1, j, ti] ) #-------------------------------|OS = A
rdm[j, 3, 2, ti] <- 0 #--------------------------------------------------|OS = B
rdm[j, 4, 2, ti] <- 0 #--------------------------------------------------|OS = AB
# TS = B
rdm[j, 1, 3, ti] <- 1 #--------------------------------------------------|OS = U
rdm[j, 2, 3, ti] <- 0 #--------------------------------------------------|OS = A
rdm[j, 3, 3, ti] <- exp( rho[2, j, ti] ) #-------------------------------|OS = B
rdm[j, 4, 3, ti] <- 0 #--------------------------------------------------|OS = AB
# TS = AB
rdm[j, 1, 4, ti] <- 1 #--------------------------------------------------|OS = U
rdm[j, 2, 4, ti] <- exp( rho[1, j, ti] ) #-------------------------------|OS = A
rdm[j, 3, 4, ti] <- exp(                  rho[2, j, ti] ) #--------------|OS = B
rdm[j, 4, 4, ti] <- exp( rho[1, j, ti] +  rho[2, j, ti] ) #--------------|OS = AB
} # close ti year
######
# Fill in the linear predictors for the transition matrices
######
for( i in 1:ncat ) {
# base occupancy
psi[i ,j] <- inprod( a[i, ], psi_cov[j, ] )
for(t in 2:nyear) {
# base colonization
gam[i, j, t-1] <- inprod( b[i, ], gam_cov[j, ,t-1] ) 
# base extinction
eps[i, j, t-1] <- inprod( d[i, ], eps_cov[j, , t-1] )
# species interactions
 pi[i, i, j] <- 0 # fill inxs diagonal with 0's
tau[i, i, j] <- 0 # fill inxs diagonal with 0's
# These are the matrics to hold conditional linear predictors for gamma and eps
#  when one other species is present at t-1 (hence the _one after the name)
gam_one[i, i, j, t-1] <- 0 # fill inxs diagonal with 0's
eps_one[i, i, j, t-1] <- 0 # fill inxs diagonal with 0's
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
     pi[rows_vec[k], cols_vec[k], j, t-1] <- inprod( g[k, ], pi_cov[j, ,t-1] )
    # inxs on extinction 
    tau[rows_vec[k], cols_vec[k], j, t-1] <- inprod( h[k, ], tau_cov[j, ,t-1] )
    # linear predictor for colonization | one species present at t-1
    gam_one[rows_vec[k], cols_vec[k], j, t-1] <- 
      inprod( b[rows_vec[k], ], gam_cov[j, ,t-1] ) + pi[rows_vec[k], cols_vec[k], j, t-1]
    # linear predictor for extinction | one species present at t-1
    eps_one[rows_vec[k], cols_vec[k], j, t-1] <- 
      inprod( d[rows_vec[k], ], eps_cov[j, ,t-1] ) + tau[rows_vec[k], cols_vec[k], j, t-1]
      } # close inxs
    } # close t
for(ti in 1:nyear){
  # base detection probability
  rho[i, j, t] <- inprod( f[i, ], rho_cov[j, , t] ) 
}

  } # close category
} # closes for loop for j (sites) all the way up at the top of the model
#####
# Priors
######
  for(i in 1:nspec){
    # Initial Occupancy
    for( psip in 1:ncov_psi ){
      a[i, psip] ~ dlogis(0, 1)
    }
    # Colonization
    for( gamp in 1:ncov_gam ){
      b[i, gamp] ~ dlogis(0, 1)
      #b0[i, gamp] ~ dlogis(0, 1)
    }  
    # Extinction
    for( epsp in 1:ncov_eps ){
      d[i, epsp] ~ dlogis(0, 1)
      #d0[i, epsp] ~ dlogis(0, 1)
    }
    # Detection
    for( rhop in 1:ncov_rho ){
      f[i, rhop] ~ dlogis(0, 1)
    }
  }
  for( k in 1:n_inxs ){
    # Inxs on colonization
    for( pip in 1:ncov_pi ){
      g[k, pip] ~ dlogis(0, 1)
    }
    # Inxs on extinction
    for( taup in 1:ncov_tau ){
      h[k, taup] ~ dlogis(0, 1)
    }
  }
}