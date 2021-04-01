################################################
#
# Fit varying models to the chicago coyote data
#
#   Written by M. Fidino
#
################################################

library(runjags)

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

# loop through them
for(mod in 1:20){
  
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
  
  if(models$inxs[mod]){
    inits <- function(chain){
      gen_list <- function(chain = chain){
        list( 
          z = matrix(
            model_list$max_state,
            ncol = model_list$nseason,
            nrow = model_list$nsite
          ),
          a = matrix(rnorm(model_list$ncat * model_list$ncov_psi),
                     ncol = model_list$ncov_psi, nrow = model_list$ncat),
          a_inxs = rnorm(model_list$ncov_psi_inxs),
          b = matrix(rnorm(model_list$ncat * model_list$ncov_gam),
                     ncol = model_list$ncov_gam, nrow = model_list$ncat),
          d = matrix(rnorm(model_list$ncat * model_list$ncov_eps),
                     ncol = model_list$ncov_eps, nrow = model_list$ncat),
          f = matrix(rnorm((model_list$ncat+model_list$to_add) * model_list$ncov_rho),
                     ncol = model_list$ncov_rho, 
                     nrow = model_list$ncat + model_list$to_add),
          g = matrix(rnorm(model_list$n_inxs * model_list$ncov_gam),
                     ncol = model_list$ncov_pi, nrow = model_list$n_inxs),
          h = matrix(rnorm(model_list$n_inxs * model_list$ncov_eps),
                     ncol = model_list$ncov_tau, nrow = model_list$n_inxs),
          .RNG.name = switch(chain,
                             "1" = "base::Wichmann-Hill",
                             "2" = "base::Marsaglia-Multicarry",
                             "3" = "base::Super-Duper",
                             "4" = "base::Mersenne-Twister",
                             "5" = "base::Wichmann-Hill",
                             "6" = "base::Marsaglia-Multicarry",
                             "7" = "base::Super-Duper",
                             "8" = "base::Mersenne-Twister"),
          .RNG.seed = sample(1:1e+06, 1)
        )
      }
      return(switch(chain,           
                    "1" = gen_list(chain),
                    "2" = gen_list(chain),
                    "3" = gen_list(chain),
                    "4" = gen_list(chain),
                    "5" = gen_list(chain),
                    "6" = gen_list(chain),
                    "7" = gen_list(chain),
                    "8" = gen_list(chain)
      )
      )
    }
    
    the_model <- "./JAGS/jags.dynamic.multistate.covars.R"
  } else{
    inits <- function(chain){
      gen_list <- function(chain = chain){
        list( 
          z = matrix(
            model_list$max_state,
            ncol = model_list$nseason,
            nrow = model_list$nsite
          ),
          a = matrix(rnorm(model_list$ncat * model_list$ncov_psi),
                     ncol = model_list$ncov_psi, nrow = model_list$ncat),
          b = matrix(rnorm(model_list$ncat * model_list$ncov_gam),
                     ncol = model_list$ncov_gam, nrow = model_list$ncat),
          d = matrix(rnorm(model_list$ncat * model_list$ncov_eps),
                     ncol = model_list$ncov_eps, nrow = model_list$ncat),
          f = matrix(rnorm((model_list$ncat+model_list$to_add) * model_list$ncov_rho),
                     ncol = model_list$ncov_rho, 
                     nrow = model_list$ncat + model_list$to_add),
          .RNG.name = switch(chain,
                             "1" = "base::Wichmann-Hill",
                             "2" = "base::Marsaglia-Multicarry",
                             "3" = "base::Super-Duper",
                             "4" = "base::Mersenne-Twister",
                             "5" = "base::Wichmann-Hill",
                             "6" = "base::Marsaglia-Multicarry",
                             "7" = "base::Super-Duper",
                             "8" = "base::Mersenne-Twister"),
          .RNG.seed = sample(1:1e+06, 1)
        )
      }
      return(switch(chain,           
                    "1" = gen_list(chain),
                    "2" = gen_list(chain),
                    "3" = gen_list(chain),
                    "4" = gen_list(chain),
                    "5" = gen_list(chain),
                    "6" = gen_list(chain),
                    "7" = gen_list(chain),
                    "8" = gen_list(chain)
      )
      )
    }
    
    the_model <- "./JAGS/jags.dynamic.multistate.null.R"
  }
  
  if(models$true_null[mod]){
    inits <- function(chain){
      gen_list <- function(chain = chain){
        list( 
          z = matrix(
            model_list$max_state,
            ncol = model_list$nseason,
            nrow = model_list$nsite
          ),
          a = rnorm(model_list$ncov_psi),
          b = rnorm(model_list$ncov_gam),
          d = rnorm(model_list$ncov_eps),
          f = rnorm(model_list$ncov_rho),
          .RNG.name = switch(chain,
                             "1" = "base::Wichmann-Hill",
                             "2" = "base::Marsaglia-Multicarry",
                             "3" = "base::Super-Duper",
                             "4" = "base::Mersenne-Twister",
                             "5" = "base::Wichmann-Hill",
                             "6" = "base::Marsaglia-Multicarry",
                             "7" = "base::Super-Duper",
                             "8" = "base::Mersenne-Twister"),
          .RNG.seed = sample(1:1e+06, 1)
        )
      }
      return(switch(chain,           
                    "1" = gen_list(chain),
                    "2" = gen_list(chain),
                    "3" = gen_list(chain),
                    "4" = gen_list(chain),
                    "5" = gen_list(chain),
                    "6" = gen_list(chain),
                    "7" = gen_list(chain),
                    "8" = gen_list(chain)
      )
      )
    }
    
    the_model <- "./JAGS/jags.dynamic.fake.multistate.R"
    
  }
  
  mout <- run.jags(
    model = the_model,
    monitor = c("a", "a_inxs", "b", "d", "f", "g", "h","z"), 
    data = model_list,
    n.chains = 4,
    inits = inits,
    adapt = 400,
    burnin = 40000,
    sample = ceiling(40000/4),
    thin = 6,
    summarise = FALSE,
    plots = FALSE,
    method = "parallel"
  )
  
  saveRDS(
    mout, 
    paste0("./Chicago coyote/model_output/coyote_model",
           stringr::str_pad(mod, 2, pad = "0"),
           "winter.RDS")
  )
}



t2 <- extend.jags(t2, sample = 500, method = "parallel", adapt = 400)

t2sum <- summary(t2, vars = c("a", "a_inxs", "b", "d", "f", "g", "h"))

round(t2sum, 2)

test <- as.matrix(as.mcmc.list(hm), chains = TRUE)

plot(test[,27] ~ test[,25], type = "p")
hey <- 
plot(test[test[,1]==3,37], type = 'l')

longshot <- as.mcmc.list(t2)

niter(longshot)
windowed.object <- window(longshot, start=70401)

wo <- combine.mcmc(windowed.object)
summary(t3, vars = c("a", "a_inxs", "b", "d", "f", "g", "h"))


t3 <- t2

t3$mcmc <- windowed.object
