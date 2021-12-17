################################################
#
# Dynamic simulations for null model
#
# written by M. Fidino
#
################################################

cat("Simulation for a dynamic multi-state null model\n")
cat("-----------------------------------------------\n")

cat("Loading run.jags...\n")
library(runjags)

# load utility functions
source(
  "./simulation study/dynamic_simulations/dynamic_utilities.R"
)

cat("Simulating covariates ...\n")
# create the covariates for the model.
my_covars <- sim_covariates()

# to see default values used to generate data:
# args(sim_covariates)

# use those covariates to simulate some data. We are not trying to do an
#  exhaustive simulation study here, so keeping species relatively easy
#  to detect and varying states (and transitions) are relatively common.
cat("Simulating parameters and data...\n")
my_data <- sim_data(
  my_covars
)


cat("Generating initial values for JAGS...\n")


# make the data_list we are going to supply to JAGS
jags_list <- jags_prep(
  my_covars,
  my_data
)

# here is the function for initial values

inits_no_inxs <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = matrix(
        jags_list$max_state,
        ncol = jags_list$nyear,
        nrow = jags_list$nsite
      ),
      a = matrix(rnorm(jags_list$ncat * jags_list$ncov_psi),
                 ncol = jags_list$ncov_psi, nrow = jags_list$ncat),
      b = matrix(rnorm(jags_list$ncat * jags_list$ncov_gam),
                 ncol = jags_list$ncov_gam, nrow = jags_list$ncat),
      d = matrix(rnorm(jags_list$ncat * jags_list$ncov_eps),
                 ncol = jags_list$ncov_eps, nrow = jags_list$ncat),
      f = matrix(rnorm(jags_list$ncat * jags_list$ncov_rho),
                 ncol = jags_list$ncov_rho, nrow = jags_list$ncat),
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
.success("Initial values generated:")

cat("Fitting model..\n")

mout <- run.jags(
  model = "./JAGS/jags.dynamic.multistate.null.R",
  monitor = c("a", "b", "d", "f"), 
  data = jags_list,
  n.chains = 4,
  inits = inits_no_inxs,
  adapt = 400,
  burnin = 50000,
  sample = ceiling(50000/4),
  thin = 10,
  summarise = FALSE,
  plots = FALSE,
  method = "parallel"
)

msum <- summary(mout)

compare_ests(
  msum = msum,
  pars = my_data$parameters
)

