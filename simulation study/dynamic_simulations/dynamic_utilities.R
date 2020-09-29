#############################################
#
# Utility functions for dynamic simulations
#
# written by M. Fidino
#
#############################################

# Logit link
logit <- function(x){
  log(x / (1 - x))
}


# inverse logit link
ilogit <- function(x){ 
  exp(x) / (1 + exp(x))
}

# a little utility fucntion to add a green checkmark if the cli
#  package exists, otherwise it will just report 'complete!'
.success <- function(x){
  if("cli" %in% (.packages())){
    cat(
      x,
      cli::col_green(symbol$tick, "\n"))
  } else {
    cat(x,"Complete!\n")
  }
}

# a function used to index interactive parameters in the model.
# In the model we put the inxs in a matrix
# these matrices get filled in this order:
#######
#0|4|5#
#-|-|-#
#1|0|6# Zeros along the diagonal
#-|-|-#
#2|3|0#
#######
# We have set it up in the model such that the species in the off-diagonal
# influence the species along the diagonal. For example, pi[3,2] is the 
# influence species 2 has on the colonization rate of species 3 while 
# pi[2, 1] is the influence species 1 has on species 2. More generally it
# is pi[species of interest, species conditioning on]. 
make_inxs <- function(nspec = NULL){
  ns <- nspec
  col_l <- row_l <- col_u <- row_u <-numeric((ns * ns - ns)/2)
  for(i in 1:(ns-1)){
    # algorithm to fill lower columns
    col_l[min(which(col_l ==0)):(min(which(col_l ==0)) + 
                                   ns-1-i)] <-  rep(i, ns-i)
    # algorithm to fill lower rows
    row_l[min(which(row_l == 0)):(min(which(row_l==0))+ 
                                    ns-1-i)] <- (i+1):ns
    # algorithm to fill upper rows
    row_u[min(which(row_u == 0)):(min(which(row_u==0))+ 
                                    ns-(ns+1)+i)] <- 1:i
    # algorithm to fill upper columns
    col_u[min(which(col_u == 0)):(min(which(col_u==0))+ 
                                    ns-(ns+1)+i)] <- rep(i + 1, i)
  }
  row_vec <- c(row_l, row_u)
  col_vec <- c(col_l, col_u)
  return(list(rows_vec = row_vec,
              cols_vec = col_vec))
}


# simulate covariates for analysis. 
#  Assumptions:
#  1. All covaraites are continuous
#  2. The covariate count INCLUDES the intercept.
#  3. Covariates apply to all states in model for a given process
sim_covariates <- function(
  nsite = 100,
  nseason = 5,
  nrep = 4,
  ncov_psi = 1,
  ncov_gam = 1,
  ncov_eps = 1,
  ncov_rho = 1,
  my_seed = NULL
) {
  if(!is.null(my_seed)){
    set.seed(my_seed)
  } else {
    my_seed <- sample(
      1e6,
      1
    )
    set.seed(my_seed)
  }
  # covariates for initial occupancy
  psi_cov <- matrix(
    1,
    ncol = ncov_psi,
    nrow = nsite
  )
  # fill-in other columns as needed
  if(ncol(psi_cov)>1){
    psi_cov[,-1] <- rnorm(
      nsite * (ncol(psi_cov) - 1)
    )
  }
  # covariates for colonization
  gam_cov <- matrix(
    1,
    ncol = ncov_gam,
    nrow = nsite
  )
  if(ncol(gam_cov)>1){
    gam_cov[,-1] <- rnorm(
      nsite * (ncol(gam_cov) - 1)
    )
  }
  # covariates for extinction
  eps_cov <- matrix(
    1,
    ncol = ncov_eps,
    nrow = nsite
  )
  if(ncol(eps_cov)>1){
    eps_cov[,-1] <- rnorm(
      nsite * (ncol(eps_cov) - 1)
    )
  }
  # covariates for detection
  rho_cov <- matrix(
    1,
    ncol = ncov_rho,
    nrow = nsite
  )
  if(ncol(psi_cov)>1){
    rho_cov[,-1] <- rnorm(
      nsite * (ncol(rho_cov) - 1)
    )
  }
  # make a list object and return.
  to_return <- list(
    psi_cov = psi_cov,
    gam_cov = gam_cov,
    eps_cov = eps_cov,
    rho_cov = rho_cov,
    nsite = nsite,
    nseason = nseason,
    nrep = nrep,
    covar_seed = my_seed
  )
  # provide some info to command line
  .success("Covariates simulated:")
  # gimme back the data
  return(
    to_return
  )
}

# function to make the parmaters for each process
sim_params <- function(
  cov_mat,
  indep_detections = FALSE
){
  # this will determine if we estimate each state seperately
  #  for detection. If so, we need an extra linear predictor.
  nstate <- ifelse(indep_detections, 4, 2)
    
    generic_mat <- matrix(
      runif(
        n = nstate * ncol(cov_mat),
        min = -1,
        max = 1
      ),
      nrow = nstate,
      ncol = ncol(cov_mat)
    )
  
  return(generic_mat)
}

# simulate null model
sim_null <- function(
  params,
  covars,
  indep_detections = FALSE
){
  #probability of each state initial occupancy
  # First season probabilities for each state
  PSI <- matrix(1, ncol = 4, nrow = covars$nsite)
  # Fill in remaining columns
  PSI[, 2] <- exp( covars$psi_cov %*% params$psi[1,] ) 
  PSI[, 3] <- exp( covars$psi_cov %*% params$psi[2,] ) 
  PSI[, 4] <- exp( 
    covars$psi_cov %*% params$psi[1,] +
    covars$psi_cov %*% params$psi[2,] 
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
    covars$gam_cov %*% params$gam[2,]
  ) 
  TPM[, 4, 2] <- exp( covars$gam_cov %*% params$gam[2,] ) 
  # N to ... U(1), D(2), DN(4)
  TPM[, 1, 3] <- exp( covars$eps_cov %*% params$eps[2,] ) 
  TPM[, 2, 3] <- exp(
      covars$gam_cov %*% params$gam[1,] +
      covars$eps_cov %*% params$eps[2,]
  ) 
  TPM[, 4, 3] <- exp( covars$gam_cov %*% params$gam[1,] ) 
  # DN to .. U(1), D(2), DN(3)
  TPM[, 1, 4] <- exp( 
    covars$eps_cov %*% params$eps[1,] + 
    covars$eps_cov %*% params$eps[2,]
  ) 
  TPM[, 2, 4] <- exp( covars$eps_cov %*% params$eps[2,] ) 
  TPM[, 3, 4] <- exp( covars$eps_cov %*% params$eps[1,] ) 
  
  # make into a probability
  for(site in 1:covars$nsite){
    TPM[site,,] <- sweep(
      TPM[site,,],
      2,
      colSums(TPM[site,,]),
      FUN = "/"
    )
  }
# Create detection matrix
  
  RDM <- array(
  0,
  dim = dim(TPM)
  )
  to_add <- ifelse(indep_detections, 2, 0)
  
  # Detection for state U 
  RDM[, 1, 1] <- 1
  # Detection for state D
  RDM[, 1, 2] <- 1
  RDM[, 2, 2] <- exp( covars$rho_cov %*% params$rho[1,] )
  # Detection for state N
  RDM[, 1, 3] <- 1 
  RDM[, 3, 3] <- exp( covars$rho_cov %*% params$rho[2,] ) 
  # Detection for State DN
  RDM[, 1, 4] <- 1
  RDM[, 2, 4] <- exp( covars$rho_cov %*% params$rho[1+to_add,] ) 
  RDM[, 3, 4] <- exp( covars$rho_cov %*% params$rho[2+to_add,] ) 
  RDM[, 4, 4] <- exp( 
    covars$rho_cov %*% params$rho[1+to_add,] +
    covars$rho_cov %*% params$rho[2+to_add,]
  )
  # make into a probability
  for(site in 1:covars$nsite){
    RDM[site,,] <- sweep(
      RDM[site,,],
      2,
      colSums(RDM[site,,]),
      FUN = "/"
    )
  }
  # simulate the data!
  z <- matrix(
    NA,
    ncol = covars$nseason,
    nrow = covars$nsite
  )
  # simulate first season
  for(site in 1:covars$nsite){
    z[site,1] <- sample(
      1:4,
      1,
      FALSE,
      PSI[site,]
    )
    # now the rest of the seasons
    for(year in 2:covars$nseason){
      z[site, year] <- sample(
        1:4,
        1,
        FALSE,
        TPM[site,,z[site,year-1]]
      )
    }
  }
  # simulate the detection data
  y <- array(
    NA,
    dim = c(covars$nsite, covars$nseason, covars$nrep)
  )
  for(site in 1:covars$nsite){
    for(year in 1:covars$nseason){
      y[site,year,] <- 
        sample(
          1:4,
          covars$nrep,
          TRUE,
          RDM[site, , z[site,year]]
        )
    }
  }
  to_return <- list(
    z = z,
    y = y
  )
  return(to_return)
}

# simulate conditional model
sim_conditional <- function(
  params,
  covars,
  indep_detections = FALSE
){
  #probability of each state initial occupancy
  # First season probabilities for each state
  PSI <- matrix(1, ncol = 4, nrow = covars$nsite)
  # Fill in remaining columns
  PSI[, 2] <- exp( covars$psi_cov %*% params$psi[1,] ) 
  PSI[, 3] <- exp( covars$psi_cov %*% params$psi[2,] ) 
  PSI[, 4] <- exp( 
    covars$psi_cov %*% params$psi[1,] +
      covars$psi_cov %*% params$psi[2,] +
      covars$psi_cov %*% params$psi_cond
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
    covars$gam_cov %*% params$gam_cond[1,]
  ) 
  TPM[, 4, 2] <- exp( 
    covars$gam_cov %*% params$gam[2,] + 
    covars$gam_cov %*% params$gam_cond[1,]
  ) 
  # N to ... U(1), D(2), DN(4)
  TPM[, 1, 3] <- exp( covars$eps_cov %*% params$eps[2,] ) 
  TPM[, 2, 3] <- exp(
    covars$gam_cov %*% params$gam[1,] +
    covars$eps_cov %*% params$eps[2,] + 
    covars$gam_cov %*% params$gam_cond[2,]
  ) 
  TPM[, 4, 3] <- exp( 
    covars$gam_cov %*% params$gam[1,] + 
    covars$gam_cov %*% params$gam_cond[2,]
  ) 
  # DN to .. U(1), D(2), DN(3)
  TPM[, 1, 4] <- exp( 
    covars$eps_cov %*% params$eps[1,] + 
    covars$eps_cov %*% params$eps[2,] +
    covars$eps_cov %*% params$eps_cond[1,] + 
    covars$eps_cov %*% params$eps_cond[2,]  
  ) 
  TPM[, 2, 4] <- exp( 
    covars$eps_cov %*% params$eps[2,] +
    covars$eps_cov %*% params$eps_cond[1,]
  ) 
  TPM[, 3, 4] <- exp( 
    covars$eps_cov %*% params$eps[1,] + 
    covars$eps_cov %*% params$eps_cond[2,]
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
  # Create detection matrix
  
  RDM <- array(
    0,
    dim = dim(TPM)
  )
  to_add <- ifelse(indep_detections, 2, 0)
  
  # Detection for state U 
  RDM[, 1, 1] <- 1
  # Detection for state D
  RDM[, 1, 2] <- 1
  RDM[, 2, 2] <- exp( covars$rho_cov %*% params$rho[1,] )
  # Detection for state N
  RDM[, 1, 3] <- 1 
  RDM[, 3, 3] <- exp( covars$rho_cov %*% params$rho[2,] ) 
  # Detection for State DN
  RDM[, 1, 4] <- 1
  RDM[, 2, 4] <- exp( covars$rho_cov %*% params$rho[1+to_add,] ) 
  RDM[, 3, 4] <- exp( covars$rho_cov %*% params$rho[2+to_add,] ) 
  RDM[, 4, 4] <- exp( 
    covars$rho_cov %*% params$rho[1+to_add,] +
    covars$rho_cov %*% params$rho[2+to_add,]
  )
  # make into a probability
  for(site in 1:covars$nsite){
    RDM[site,,] <- sweep(
      RDM[site,,],
      2,
      colSums(RDM[site,,]),
      FUN = "/"
    )
  }
  # simulate the data!
  z <- matrix(
    NA,
    ncol = covars$nseason,
    nrow = covars$nsite
  )
  # simulate first season
  for(site in 1:covars$nsite){
    z[site,1] <- sample(
      1:4,
      1,
      FALSE,
      PSI[site,]
    )
    # now the rest of the seasons
    for(year in 2:covars$nseason){
      z[site, year] <- sample(
        1:4,
        1,
        FALSE,
        TPM[site,,z[site,year-1]]
      )
    }
  }
  # simulate the detection data
  y <- array(
    NA,
    dim = c(covars$nsite, covars$nseason, covars$nrep)
  )
  for(site in 1:covars$nsite){
    for(year in 1:covars$nseason){
      y[site,year,] <- 
        sample(
          1:4,
          covars$nrep,
          TRUE,
          RDM[site, , z[site,year]]
        )
    }
  }
  to_return <- list(
    z = z,
    y = y
  )
  return(to_return)
}
# A function to simulate all of the data for the varying models.
#  Assumptions:
#  1. Processes vary spatial not temporal (i.e., covariate matrices
#     are two dimensional).
#  2. You simulated the data with sim_covariates()

sim_data <- function(
  covar_list,
  conditional_params = FALSE,
  independent_detections = FALSE
){
  # set seed based on what is in the covar_list
  if('covar_seed' %in% names(covar_list)){
    my_seed <- covar_list$covar_seed + 1
  } else {
    my_seed <- sample(
      1e6,
      1
    )
  }
  set.seed(my_seed)
  # generate parameters for initial occupancy
  psi <- sim_params(
    covar_list$psi_cov
  )
  # generate parameters for colonization
  gam <- sim_params(
    covar_list$gam_cov
  )
  # generate parameters for extinction
  eps <- sim_params(
    covar_list$eps_cov
  )
  # generate parmeters for detections
  rho <- sim_params(
    covar_list$rho,
    indep_detections = independent_detections
  )
  if(conditional_params){
    # we only need the first row here
    psi_cond <- sim_params(
      covar_list$psi
    )[1,]
    gam_cond <- sim_params(
      covar_list$gam_cov
    )
    eps_cond <- sim_params(
      covar_list$eps_cov
    )
  }
  .success("Parameters simulated:")
  # simulate states for null model
  if(!conditional_params){
    pars_list <- list(
      psi = psi,
      gam = gam,
      eps = eps,
      rho = rho
    )
    data_list <- sim_null(
      pars_list,
      covar_list,
      indep_detections = independent_detections
    )
  } else {
    pars_list <- list(
      psi = psi,
      gam = gam,
      eps = eps,
      rho = rho,
      psi_cond = psi_cond,
      gam_cond = gam_cond,
      eps_cond = eps_cond
    )
    data_list <- sim_conditional(
      pars_list,
      covar_list,
      indep_detections = independent_detections
    )
  }
  .success("Data simulated:")
  to_return <- list(
    parameters = pars_list,
    data = data_list
  )
  return(to_return)
}

# this function just combines the other lists we created.
jags_prep <- function(
  covar_list,
  data_list,
  has_conditionals = FALSE,
  independent_detections = FALSE
){
  if(independent_detections){
    to_add <- 2
  } else {
    to_add <- 0
  }
  # make arrays if needed
  if(length(dim(covar_list$gam_cov)) == 2){
    covar_list$gam_cov <- array(
      as.numeric(covar_list$gam_cov),
      dim = c(
        covar_list$nsite,
        dim(covar_list$gam_cov)[2],
        (covar_list$nseason-1)
      )
    )
  }
  if(length(dim(covar_list$eps_cov)) == 2){
    covar_list$eps_cov <- array(
      as.numeric(covar_list$eps_cov),
      dim = c(
        covar_list$nsite,
        dim(covar_list$eps_cov)[2],
        (covar_list$nseason-1)
      )
    )
  }
  if(length(dim(covar_list$rho_cov)) == 2){
    covar_list$rho_cov <- array(
      as.numeric(covar_list$rho_cov),
      dim = c(
        covar_list$nsite,
        dim(covar_list$rho_cov)[2],
        covar_list$nseason
      )
    )
  }
  
  to_return<- list(
    # count of number of paramters
    ncov_psi = ncol(covar_list$psi_cov),
    ncov_gam = dim(covar_list$gam_cov)[2],
    ncov_eps = ncol(covar_list$eps_cov),
    ncov_rho = ncol(covar_list$rho_cov),
    # number of categories, going to be 2 regardless (Day and Night)
    ncat = 2,
    # max possible state in model, which is always going to be 4
    max_state = 4,
    # design matrices
    psi_cov = covar_list$psi_cov,
    gam_cov = covar_list$gam_cov,
    eps_cov = covar_list$eps_cov,
    rho_cov = covar_list$rho_cov,
    # count of sites, seasons, surveys
    nsite = covar_list$nsite,
    nseason = covar_list$nseason,
    nsurvey = covar_list$nrep,
    # observed data
    y = data_list$data$y,
    # 0 if non-independent detection, 2 if independent detection
    to_add = to_add
  )
  if(has_conditionals){
    to_return$n_inxs <- 2
    to_return$ncov_psi_inxs <- to_return$ncov_psi
    to_return$ncov_pi <- to_return$ncov_gam
    to_return$ncov_tau <- to_return$ncov_eps
    tmp <- make_inxs(2)
    to_return$rows_vec <- tmp$rows_vec
    to_return$cols_vec <- tmp$cols_vec
    to_return$psi_inxs_cov <- to_return$psi_cov
    to_return$pi_cov <- to_return$gam_cov
    to_return$tau_cov <- to_return$eps_cov
  }
  
  return(to_return)
}

compare_ests <- function(
  msum,
  pars,
  conditional_params = FALSE
){

  if(conditional_params){
    pnames <- list(
      psi = "a",
      gam = "b",
      eps = "d",
      rho = "f",
      psi_cond = "a_inxs",
      gam_cond = "g",
      eps_cond = "h"
    )
  } else {
    pnames <- list(
      psi = "a",
      gam = "b",
      eps = "d",
      rho = "f"
    )
  }
  par_split <- vector("list", length= length(pnames))
  for(i in 1:length(pnames)){
    if(pnames[[i]] == "a"){
      par_split[[i]] <- msum[grep(
        paste0(pnames[i],"\\["), row.names(msum)),1:3]
    } else {
      par_split[[i]] <- msum[grep(
        pnames[i], row.names(msum)),1:3]
    }
    
  }
  names(par_split) <- names(pnames)
  
  # generate plot of estimates compares to true values
  windows(10,10)
  par(mfrow = c(3,3))
  
  for(i in 1:length(pnames)){
    if(is.null(nrow(par_split[[i]]))){
      my_y <- 1
      par_split[[i]] <- t(par_split[[i]])
    } else {
      my_y <- seq_len(nrow(par_split[[i]]))
    }
    plot(
      x = par_split[[i]][,2],
      y = my_y,
      pch = 15,
      col = "blue",
      xlim = c(-2,2),
      cex = 2,
      bty = 'l',
      xlab = "Estimate",
      ylab = "",
      yaxt = "n",
      main = names(pnames)[i]
    )
    for(j in 1:nrow(par_split[[names(pnames)[i]]])){
      lines(
        y = rep(j,2),
        x = par_split[[i]][j,c(1,3)],
        col = "blue",
        lwd = 2
      )
    }
    points(
      x = pars[[names(pnames)[i]]],
      y = my_y,
      col = "grey50",
      pch = 19,
      cex = 1.5
    )
    mtext(
      row.names(par_split[[i]]),
      2,
      at = seq_len(nrow(par_split[[i]])),
      las = 1,
      line = 0.75
    )
    if(i == length(pnames)){
    legend(
      "right",
      col = c("blue", "grey50"),
      bty = "n",
      legend = c("Estimate", "Truth"),
      pch = c(15,19),
      pt.cex = c(2,1.5)
    )
    }
  }
}


calc_binco <- function(
  covar_list,
  data_list
){
  # this is the number of times each state was observed
  # per season. nsite x nseason x total number of states. Each cell holds
  # the number of times that state was observed. This is used
  # to generate the binomial coefficient for model selection.
  state_count <- array(
    0,
    dim = c(
      covar_list$nsite,
      covar_list$nseason,
      4
    )
  )
  for(i in 1:covar_list$nsite){
    for(j in 1:covar_list$nseason){
      a <- table(data_list$data$y[i,j,])
      if(length(a)==0) {
        next 
      } else {
        state_count[i,j, as.numeric(names(a))] <- as.numeric(a)
      }
    }
  }

  # The number of observations that occured at a site and season
  n_obs <- apply(
    state_count,
    c( 1, 2 ),
    sum 
  )
  # The numerator of the binomial coefficient
  numerator <- factorial( n_obs )
  # an indicator variable that will convert an NA to 0
  # Which will, in turn zero out the binomial coefficent
  # when no samples were taken so that estimated states do not influence the 
  # likelihood of the model.
  ind_var <- n_obs 
  ind_var[ind_var > 0] <- 1
  # zero out numerator if NA
  numerator[is.na(numerator)] <- 0
  # give denominator the same dimensions as numerator
  denominator <- n_obs
  
  # the product of the factorial of state_count is the denominator
  for(i in 1:covar_list$nsite){
    for( j in 1:covar_list$nseason){
      denominator[i, j] <- prod( factorial(state_count[i,j, ] ) )
    }
  }

  # zero it out if no samples were taken
  binco <- (numerator / denominator) * ind_var
  
  to_return <- list(
    state_count = state_count,
    binco = binco
  )
  
  return(to_return)

}

# a function to calculate the CPO for the null model

calculate_cpo <- function(
  mm,
  list2jags,
  binco_list
){
  if(!"n_inxs" %in% list2jags) list2jags$n_inxs <- list2jags$ncat^2-list2jags$ncat
  
  
  # ensure mm is a matrix
  if(!is.matrix(mm)) {
    stop("The object mm must be a matrix")
  }
  
  # This is the likelihood matrix, which stores the likelihood of each
  # observations based off the parameters in the model for each step
  # of the mcmc chain.
  lik_mat <- matrix(
    0, 
    ncol = list2jags$nseason * list2jags$nsite,
    nrow = nrow(mm)
  )
  state_count <- binco_list$state_count
  binco <- binco_list$binco
  for(o in 1:nrow(mm)){ # going through each step of the mcmc chain
    
    # initial occupancy
    a <- matrix(
      mm[o,grep("a\\[", colnames(mm))],
      ncol = list2jags$ncov_psi,
      nrow = list2jags$ncat
    )
    # colonization 
    b <- matrix(
      mm[o,grep("b\\[", colnames(mm))],
      ncol = list2jags$ncov_gam,
      nrow = list2jags$ncat
    )
    # extinction
    d <- matrix(
      mm[o,grep("d\\[", colnames(mm))],
      ncol = list2jags$ncov_eps,
      nrow = list2jags$ncat
    )
    # detection
    f <- matrix(
      mm[o,grep("f", colnames(mm))],
      ncol = list2jags$ncov_rho,
      nrow = list2jags$ncat + list2jags$to_add
    )
    #psi inxs
    a_inxs <- matrix(
      mm[o, grep("a_inxs", colnames(mm))],
      ncol = list2jags$ncov_psi_inxs,
      nrow = 1
    )
    if(sum(is.na(a_inxs)>0)) a_inxs[is.na(a_inxs)] <- 0
    # colonization inxs
    g <- matrix(
      mm[o,grep("g", colnames(mm))],
      ncol = list2jags$ncov_pi,
      nrow = list2jags$n_inxs
    ) 
    if(sum(is.na(g)>0)) g[is.na(g)] <- 0 # make 0 if not in the model
    # extinciton inxs
    h <- matrix(
      mm[o,grep("h", colnames(mm))],
      ncol = list2jags$ncov_tau,
      nrow = list2jags$n_inxs
    )
    if(sum(is.na(h)>0)) h[is.na(h)] <- 0 # make 0 if not in the model
    # detection inxs (remove -1 if doing all species)
    #l <- matrix(mm[o,grep("l", colnames(mm))], ncol = ncov_delta, nrow = nspec - 1)
    #if(sum(is.na(l)>0)) l[is.na(l)] <- 0 # make 0 if not in the model
    # estimated community state at time t and site j
    z <- matrix(
      mm[o,grep("z", colnames(mm))],
      ncol = list2jags$nseason,
      nrow = list2jags$nsite
    )
    
    # occupancy linear predictor
    psinit  <- a %*% t(list2jags$psi_cov)
    # occupacny at fourth state
    psi_one <- colSums(psinit) + a_inxs %*% t(list2jags$psi_inxs_cov)
    # colonization linear predictor
    gam <- eps <- array(
      NA,
      dim = c(
        list2jags$ncat,
        list2jags$nsite,
        list2jags$nseason-1
      )
    )
    for(yr in 2:list2jags$nseason){
      gam[,,yr-1] <- b %*% list2jags$gam_cov[,,yr-1]
      eps[,,yr-1] <- d %*% list2jags$eps_cov[,,yr-1]
    }
    rho <- array(
      NA,
      dim = c(
        list2jags$ncat + list2jags$to_add,
        list2jags$nsite,
        list2jags$nseason
      )
    )
    # detection linear predictor
    for(yr in 1:list2jags$nseason){
    rho[,,yr] <-f %*% list2jags$rho_cov[,,yr]
    }
    # set up arrays for the species interactions.
    # these arrays get filled in this order:
    #######
    #0|4|5#
    #-|-|-#
    #1|0|6# Zeros along the diagonal
    #-|-|-#
    #2|3|0#
    #######
    # We have set it up in the model such that the species in the off-diagonal
    # influence the species along the diagonal. For example, pi[3,2] is the 
    # influence species 2 has on the colonization rate of species 3 while 
    # pi[2, 1] is the influence species 1 has on species 2. More generally it
    # is pi[species of interest, species conditioning on]. 
    pi <- tau <-  gam_one <- eps_one <- array(
      0,
      dim = c(
        list2jags$ncat,
        list2jags$ncat,
        list2jags$nsite,
        list2jags$nseason-1
      )
    )
    for(yr in 2:list2jags$nseason){
    for(k in 1:list2jags$n_inxs){
      # inxs on colonization 
      pi[list2jags$rows_vec[k], list2jags$cols_vec[k], ,yr-1] <-  
        list2jags$pi_cov[j, ,yr-1] %*% g[k, ]
      # inxs on extinction 
      tau[list2jags$rows_vec[k], list2jags$cols_vec[k], ,yr-1] <-
        list2jags$tau_cov[j, ,yr-1] %*% h[k, ] 
      # linear predictor for colonization | one species present at t-1
      gam_one[list2jags$rows_vec[k], list2jags$cols_vec[k], ,yr-1] <- 
        gam[list2jags$rows_vec[k],,yr-1] + 
        pi[list2jags$rows_vec[k], list2jags$cols_vec[k], ,yr-1]
      # linear predictor for extinction | one species present at t-1
      eps_one[list2jags$rows_vec[k], list2jags$cols_vec[k], ,yr-1] <- 
        eps[list2jags$rows_vec[k],,yr-1] +
        tau[list2jags$rows_vec[k], list2jags$cols_vec[k],,yr-1 ]
    }
    }
    
    PSI <- matrix(
      1, 
      ncol = list2jags$max_state,
      nrow = list2jags$nsite
    )
    ######
    # Fill in all of the transition probabilities
    ######
    ######
    # Latent state for first season, fsm = transition vector for season 1
    ######
    # First season probabilities for each state
    PSI[, 1] <- 1 #----------------------------------------------------|U
    PSI[, 2] <- exp( psinit[1, ] ) #-----------------------------------|D
    PSI[, 3] <- exp( psinit[2, ] ) #-----------------------------------|N
    PSI[, 4] <- exp( psi_one[i]  ) #------------------------------------|DN
    
    tpm <- array(
      1,
      dim = c(
        list2jags$nsite,
        list2jags$max_state,
        list2jags$max_state,
        list2jags$nseason - 1
      )
    )
    for(yr in 2:list2jags$nseason) {
      ######
      # Latent state for dynamic part of model
      # tpm = transition probability matrix. All columns sum to 1.
      # dim(tpm)[1] = site j 
      # dim(tpm)[2] = state at time t
      # dim(tpm)[3] = state at time t-1
      # dim(tpm)[4] = time t (for temporally varying covariates)
      ######
      # U to ...
      tpm[, 1, 1, yr-1] <- 1 #-----------------------------------------------|U
      tpm[, 2, 1, yr-1] <- exp( gam[1, , yr-1] ) #---------------------------|D
      tpm[, 3, 1, yr-1] <- exp(                  gam[2, , yr-1] ) #----------|N
      tpm[, 4, 1, yr-1] <- exp( gam[1, , yr-1] + gam[2, , yr-1] ) #----------|DN
      # D to ...
      tpm[, 1, 2, yr-1] <- exp( eps[1, , yr-1] ) #---------------------------|U
      tpm[, 2, 2, yr-1] <- 1 #-----------------------------------------------|D
      tpm[, 3, 2, yr-1] <- exp( eps[1, , yr-1] + gam_one[2, 1, , yr-1] ) #---|N
      tpm[, 4, 2, yr-1] <- exp(                  gam_one[2, 1, , yr-1] ) #---|DN
      # N to ...
      tpm[, 1, 3, yr-1] <- exp(                         eps[2, , yr-1] ) #---|U
      tpm[, 2, 3, yr-1] <- exp( gam_one[1, 2, , yr-1] + eps[2, , yr-1] ) #---|D
      tpm[, 3, 3, yr-1] <- 1 #-----------------------------------------------|N
      tpm[, 4, 3, yr-1] <- exp( gam_one[1, 2, , yr-1] ) #--------------------|DN
      # DN to ..
      tpm[, 1, 4, yr-1] <- exp( eps_one[1, 2, , yr-1] + eps_one[2, 1, , yr-1] ) #-|U
      tpm[, 2, 4, yr-1] <- exp( eps_one[2, 1, , yr-1] ) #-------------------------|D
      tpm[, 3, 4, yr-1] <- exp(                         eps_one[1, 2, , yr-1] ) #-|N
      tpm[, 4, 4, yr-1] <- 1 #----------------------------------------------------|DN
    } # close yr year loop
    ######
    # detection matrix (OS = observed state, TS = true state)
    # rdm = rho detection matrix. Each row sums to 1.
    # OS along rows, TS along columns
    ######
    # TS = U
    rdm <- array(
      0,
      dim = c(
        list2jags$nsite,
        list2jags$max_state,
        list2jags$max_state,
        list2jags$nseason
      )
    )
    for(ti in 1:list2jags$nseason) {
      rdm[, 1, 1, ti] <- 1 #-------------------------------------------------|OS = U
      rdm[, 2, 1, ti] <- 0 #-------------------------------------------------|OS = D
      rdm[, 3, 1, ti] <- 0 #-------------------------------------------------|OS = N
      rdm[, 4, 1, ti] <- 0 #-------------------------------------------------|OS = DN
      # TS = D
      rdm[, 1, 2, ti] <- 1 #-------------------------------------------------|OS = U
      rdm[, 2, 2, ti] <- exp( rho[1, , ti] ) #------------------------------|OS = D
      rdm[, 3, 2, ti] <- 0 #-------------------------------------------------|OS = N
      rdm[, 4, 2, ti] <- 0 #-------------------------------------------------|OS = DN
      # TS = N
      rdm[, 1, 3, ti] <- 1 #-------------------------------------------------|OS = U
      rdm[, 2, 3, ti] <- 0 #-------------------------------------------------|OS = D
      rdm[, 3, 3, ti] <- exp( rho[2, , ti] ) #------------------------------|OS = N
      rdm[, 4, 3, ti] <- 0 #-------------------------------------------------|OS = DN
      # TS = DN
      rdm[, 1, 4, ti] <- 1 #-------------------------------------------------|OS = U
      rdm[, 2, 4, ti] <- exp( rho[1+list2jags$to_add, , ti] ) #-----------------------|OS = D
      rdm[, 3, 4, ti] <- exp( rho[2+list2jags$to_add, , ti] ) #------|OS = N
      rdm[, 4, 4, ti] <- exp( rho[1+list2jags$to_add, , ti] +
                              rho[2+list2jags$to_add, , ti] )#|OS = DN
    } # close ti year
    
    # This is the probability of detecting each community state per site
    pstate <- array(
      0,
      dim = c(
        list2jags$nsite,
        list2jags$max_state,
        list2jags$nseason
        )
    )
    # A likelihood matrix to feed back into the full MCMC chain lik_mat
    lik <- matrix(
      0,
      ncol = list2jags$nseason,
      nrow = list2jags$nsite
    )
    
    for(j in 1:list2jags$nsite){
      
      # Product of this times binomial coefficient = 
      # multinomial likelihood for detection
      # This is for the first year, which uses the fsm matrix below
      # divide by the appropriate column sum (indexed by z matrix)
      # to get the probability, which is then raised to the number
      # of times each state was observed.
      pstate[j, , 1] <- (rdm[j, ,z[j, 1] ,1] /
                           sum(rdm[j, ,z[j, 1],1 ] ) )^state_count[j, 1, ]
      
      # fsm indexed to z to get likelihood of latent state
      # at first season.
      
      lik[j,1] <- ( binco[j, 1] * prod( pstate[j, ,1] ) ) * 
        ( PSI[j, z[j, 1]] / sum(PSI[j, ] ) )
      for( yr in 2:nseason ){
        # Product of this x binomial coefficient = 
        # multinomial likelihood for detection
        pstate[j, , yr] <-  (rdm[j, , z[j, yr],yr ]/ 
                              sum( rdm[j, , z[j, yr],yr ]) )^ state_count[j, yr, ] 
        # tpm indexed to site, state at times t, state at t-1 to get likelihood
        # for latent state
        lik[j, yr] <- ( binco[j, yr] * prod( pstate[j, ,yr] ) ) * 
          ( tpm[ j, z[ j, yr ] , z[ j, yr-1 ],yr ] / # numerator smax
              sum( tpm[ j, ,z[j, yr-1 ],yr ] ) ) # denominator smax
      } # close t (year)
    } # close j (site)
    lik_mat[o,] <- as.numeric(lik) # put lik into likelihood matrix
  }
  
  # these are sites by seasons that sampling did not happen
  # we do not want them to contribute towards the likelihood
  to_zero <- rep(0, ncol(lik_mat))
  
  for(i in 1:length(to_zero)){
    to_zero[i] <- sum(lik_mat[,i]==0)
  }
  # these are all unobserved sites
  to_go <- which(to_zero == nrow(mm))
  # remove them from the likelihood matrix
  lik_mat <- lik_mat[,-to_go]
  
  # calculate cpo
  CPO <--sum(log(nrow(mm) /(apply(1/lik_mat, 2, sum))))
  # return cpo
  return(CPO)
}