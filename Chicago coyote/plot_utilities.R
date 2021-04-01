
# calculates the transitions probability matrix
#  for the urbanization covariate
predict_latent_conditional <- function(
  mcmc,
  data_list
){
  
  # initial occupancy
  a <- matrix(
    mcmc[grep("a\\[", names(mcmc))], 
    ncol = data_list$ncov_psi,
    nrow = data_list$ncat
  )
  # a_inxs 
  a_inxs <- matrix(
    mcmc[grep("a_inxs", names(mcmc))],
    ncol = data_list$ncov_psi_inxs,
    nrow = 1
  )
  # colonization 
  b <- matrix(
    mcmc[grep("b\\[", names(mcmc))],
    ncol = data_list$ncov_gam,
    nrow = data_list$ncat
  )
  # colonization 
  # extinction
  d <- matrix(
    mcmc[grep("d\\[", names(mcmc))],
    ncol = data_list$ncov_eps,
    nrow = data_list$ncat
  )
  # colonization inxs
  g <- matrix(
    mcmc[grep("g\\[", names(mcmc))],
    ncol = data_list$ncov_pi,
    nrow = data_list$n_inxs
  ) 
  if(sum(is.na(g)>0)) g[is.na(g)] <- 0 # make 0 if not in the model
  # extinciton inxs
  h <- matrix(
    mcmc[grep("h", names(mcmc))],
    ncol = data_list$ncov_tau, nrow = data_list$n_inxs
  )
  if(sum(is.na(h)>0)) h[is.na(h)] <- 0 # make 0 if not in the model
  
  
  params <- list(
    psi = a,
    gam = b,
    eps = d,
    psi_cond = a_inxs,
    gam_cond = g,
    eps_cond = h
  )
  # The urbanization covariate for predictions
  xp <- matrix(1, ncol = 2, nrow = 250)
  xp[,2] <- seq(-4, 4, length.out = 250)
  covars <- list(
    psi_cov = xp,
    gam_cov = xp,
    eps_cov = xp,
    pi_cov = matrix(1, ncol = 1, nrow = nrow(xp)),
    tau_cov = matrix(1, ncol = 1, nrow = nrow(xp)),
    nsite = nrow(xp)
  )
  #probability of each state initial occupancy
  # First season probabilities for each state
  PSI <- matrix(1, ncol = 4, nrow = covars$nsite)
  # Fill in remaining columns
  PSI[, 2] <- exp( covars$psi_cov %*% params$psi[1,] ) 
  PSI[, 3] <- exp( covars$psi_cov %*% params$psi[2,] ) 
  PSI[, 4] <- exp( 
    covars$psi_cov %*% params$psi[1,] +
      covars$psi_cov %*% params$psi[2,] +
      covars$psi_cov %*% t(params$psi_cond)
  )
  # convert to probability
  PSI <- sweep(
    PSI,
    1,
    rowSums(PSI),
    FUN = "/"
  )
  # Latent state for dynamic part of model
  # TPM = transition probability matrix. All columns sum to 1.
  # dim(TPM)[1] = site i 
  # dim(TPM)[2] = state at time t
  # dim(TPM)[3] = state at time t-1
  TPM <- array(
    1,
    dim = c(
      covars$nsite,
      4,
      4
    )
  )
  # U to ... D(2), N(3), DN(4)
  TPM[, 2, 1] <- exp( covars$gam_cov %*% params$gam[1,] ) 
  TPM[, 3, 1] <- exp( covars$gam_cov %*% params$gam[2,] ) 
  TPM[, 4, 1] <- exp( 
    covars$gam_cov %*% params$gam[1,] + 
      covars$gam_cov %*% params$gam[2,] 
  ) 
  # D to ...U(1), N(3), DN(4)
  TPM[, 1, 2] <- exp( covars$eps_cov %*% params$eps[1,] ) 
  TPM[, 3, 2] <- exp( 
    covars$eps_cov %*% params$eps[1,] +
      covars$gam_cov %*% params$gam[2,] + 
      covars$pi_cov %*% params$gam_cond[1,]
  ) 
  TPM[, 4, 2] <- exp( 
    covars$gam_cov %*% params$gam[2,] + 
      covars$pi_cov %*% params$gam_cond[1,]
  ) 
  # N to ... U(1), D(2), DN(4)
  TPM[, 1, 3] <- exp( covars$eps_cov %*% params$eps[2,] ) 
  TPM[, 2, 3] <- exp(
    covars$gam_cov %*% params$gam[1,] +
      covars$eps_cov %*% params$eps[2,] + 
      covars$pi_cov %*% params$gam_cond[2,]
  ) 
  TPM[, 4, 3] <- exp( 
    covars$gam_cov %*% params$gam[1,] + 
      covars$pi_cov %*% params$gam_cond[2,]
  ) 
  # DN to .. U(1), D(2), DN(3)
  TPM[, 1, 4] <- exp( 
    covars$eps_cov %*% params$eps[1,] + 
      covars$eps_cov %*% params$eps[2,] +
      covars$tau_cov %*% params$eps_cond[1,] + 
      covars$tau_cov %*% params$eps_cond[2,]  
  ) 
  TPM[, 2, 4] <- exp( 
    covars$eps_cov %*% params$eps[2,] +
      covars$tau_cov %*% params$eps_cond[1,]
  ) 
  TPM[, 3, 4] <- exp( 
    covars$eps_cov %*% params$eps[1,] + 
      covars$tau_cov %*% params$eps_cond[2,]
  ) 
  
  # make into a probability
  for(site in 1:covars$nsite){
    TPM[site,,] <- sweep(
      TPM[site,,],
      2,
      colSums(TPM[site,,]),
      FUN = "/"
    )
  }
  return(TPM)
}
