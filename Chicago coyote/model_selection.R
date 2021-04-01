################################################
#
# Do model selection with the fitted models
#
#   Written by M. Fidino
#
################################################

library(runjags)
library(coda)

# The simulation code has some helpful utility functions we can
#  use to set up the analysis.

source("./simulation study/dynamic_simulations/dynamic_utilities.R")

# and there are also a couple ones for this analysis here

source("./Chicago coyote/analysis_utilities.R")

# This is the code that sets up the models we are going to fit
#  and brings in the data.

source("./Chicago coyote/prep_objects_for_model.R")
#---Objects created by prep_objects_for_model.R---#

#  1. covs: The covariate data
#  2. y: The coyote detection non-detection data
#  3. models: The linear predictors for the 20 different models we are fitting.
#  4. nseason: number of seasons sampled
#  5. nsite: number of sites sampled
#  6. nsurvey: Number of secondary sampling weeks within a season.

# Calculate the binomial coefficient for model selection. This function
#  was built for the simulation study so we are creating the lists
#  in the format it needs.

coyote_binco <- calc_binco(
  list(
    nsite = nsite,
    nseason = nseason
  ),
  list(
    data = list(
      y = y
    )
  )
)

# get model outputs
mod_paths <- list.files(
  "./Chicago coyote/model_output/",
  full.names = TRUE
)
mod_paths[9] <- "./Chicago coyote/simpler_models/full_model_urb.RDS"
models$CPO <- NA


# loop through them
for(mod in c(9, 14, 19)){
  cat(paste0("\nCPO for model ", mod, "\n"))
  
  # read in model matrix
  mod_mcmc <- as.matrix(
    as.mcmc.list(
      readRDS(
        mod_paths[mod]
      )
    ),
    chains = TRUE
  )[,-1]

  
  gam_cov <- eps_cov <- pi_cov <- tau_cov <- make_model_matrix(
    formula = models$formula[mod],
    df = coyote,
    rho = FALSE,
    comp_winter = TRUE
  )
  pi_cov <- tau_cov <-  array(1, dim = c(105,1,12))

  rho_cov <- make_model_matrix(
    formula = models$formula[mod],
    df = coyote,
    rho = TRUE,
    comp_winter = TRUE
  )
  has_urb <- length(grep("urb", models$formula[mod])) > 0
  if(has_urb){
    cols_to_use <- c(1, which(dimnames(rho_cov)[[2]] == "urb"))
    psi_cov  <- rho_cov[,cols_to_use,1]
  } else {
    psi_cov  <- matrix(
      rho_cov[,1,1],
      ncol = 1,
      nrow = nsite
    )
  }
  
  model_list <- list(
    ncov_psi = dim(psi_cov)[2],
    ncov_gam = dim(gam_cov)[2],
    ncov_eps = dim(eps_cov)[2],
    ncov_rho = dim(rho_cov)[2],
    ncat = 2,
    max_state = 4,
    psi_cov = psi_cov,
    gam_cov = gam_cov,
    eps_cov = eps_cov,
    rho_cov = rho_cov,
    nsite = nsite,
    nseason = nseason,
    nsurvey = nsurvey,
    y = y,
    to_add = ifelse(models$indep_rho[mod], 2, 0) , # 2 because we have independent rho
    n_inxs = 2,
    ncov_psi_inxs = dim(psi_cov)[2],
    ncov_pi = dim(pi_cov)[2],
    ncov_tau = dim(tau_cov)[2],
    rows_vec = make_inxs(2)$rows_vec,
    cols_vec = make_inxs(2)$cols_vec,
    psi_inxs_cov = psi_cov,
    pi_cov = pi_cov,
    tau_cov = tau_cov
  )
  
  # if(models$inxs[mod]){
  #   
  #   the_model <- "./JAGS/jags.dynamic.multistate.covars.R"
  # } else{
  # 
  #   
  #   the_model <- "./JAGS/jags.dynamic.multistate.null.R"
  # }
  # 
  # if(models$true_null[mod]){
  #   
  #   the_model <- "./JAGS/jags.dynamic.fake.multistate.R"
  #   
  # }
  if(!models$true_null[mod]){
  models$CPO[mod] <- calculate_cpo(
    mm = mod_mcmc,
    list2jags = model_list,
    binco_list = coyote_binco
  )} else{
    models$CPO[mod] <- calculate_cpo_null(
      mm = mod_mcmc,
      list2jags = model_list,
      binco_list = coyote_binco
    )
  }


}
hm <- readRDS("./Chicago coyote/simpler_models/full_model_urb.RDS")

yo <-summary(mout, vars = c("a", "a_inxs", "b", "d", "f", "g", "h"))
yo <-summary(hm, vars = c("a" , "b", "d", "f"))
round(yo, 2)

aa <- as.matrix(as.mcmc.list(hm))

windows(20,20)

pairs(test[,20:27])
