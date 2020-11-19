############################################
#
# Simulating data for the null dynamic model
#
# Written by M. Fidino 8/20/2020
#
############################################

library(extraDistr)

# The model parameters we need to simulate.
# a = psi_covs [1:2, ncov_psi] = initial occupancy
# b = gam_covs [1:2, ncov_gam] = colonization
# d = eps_covs [1:2, ncol_eps] = extinction
# f = rho_covs [1:2, ncov_rho] = detection

#########################
# Some utility functions
logit <- function(p) log(p / (1-p))
cp <- function(p) 1 - p
#########################

# Number of parameters for each linear predictor, including 
#  intercept.
ncov_psi <- 2
ncov_gam <- 2
ncov_eps <- 2
ncov_rho <- 2

# number of sites, years, and surveys per site (per year)
nsite <- 200
nyear <- 5
nsurvey <- 4

# for reproducibility
set.seed(155)

# simulate the covariates
a <- matrix(
  logit(
    rbeta(2 * ncov_psi, 4, 4)
  ),
  ncol = ncov_psi,
  nrow = 2
  )

b <- matrix(
  logit(
    rbeta(2 * ncov_gam, 4, 4)
  ),
  ncol = ncov_gam,
  nrow = 2
)

d  <- matrix(
  logit(
    rbeta(2 * ncov_eps, 4, 4)
  ),
  ncol = ncov_eps,
  nrow = 2
)

f <- matrix(
  logit(
    rbeta(2 * ncov_rho, 4, 4)
  ),
  ncol = ncov_rho,
  nrow = 2
)

# using the same covariate throughout all linear predictors,
#  assume it's temporally varying.

x <- array(
  1,
  dim = c(nsite, 2, nyear)
)
x[,2,] <- rnorm(nsite * nyear)

# psi_cov is first year, eps_cov and gam_cov 2:4, rho_cov all
psi_cov <- x[,,1]
gam_cov <- x[,,2:nyear]
eps_cov <- x[,,2:nyear]
rho_cov <- x

# Simulate first year.
PSI <- matrix(
  NA,
  ncol = 4,
  nrow = nsite
)
tpm <- array(
  NA,
  dim = c(nsite, 4, 4, nyear-1)
)

PSI[,1] <- 1
PSI[,2:3] <- exp( 
  psi_cov %*% t(a)
)
PSI[,4] <- exp(
  rowSums( 
    psi_cov %*% t(a)
  )
)

# Fill out transition probabilities
for(i in 2:nyear){
  # U to...
  tpm[,1,1,i-1]  <- 1
  tpm[,2:3,1,i-1] <- exp(
    gam_cov[,,i-1] %*% t(b)
  )
  tpm[,4,1,i-1] <- exp(
    rowSums(
      gam_cov[,,i-1] %*% t(b)
    )
   )
  # D to ...
  tpm[, 1, 2, i-1] <- exp(
    eps_cov[,,i-1] %*% b[1,] 
  ) 
  tpm[, 2, 2, i-1] <- 1 
  tpm[, 3, 2, i-1] <- exp( 
    eps_cov[,,i-1] %*% d[1,] + gam_cov[,,i-1] %*% b[2,]
  )
  tpm[, 4, 2, i-1] <- exp( 
    gam_cov[,,i-1] %*% b[2,] 
  )
  # N to ...
  tpm[, 1, 3, i-1] <- exp(
    eps_cov[,,i-1] %*% b[2,] 
  )
  tpm[, 2, 3, i-1] <- exp( 
    gam_cov[,,i-1] %*% b[1,] + eps_cov[,,i-1] %*% d[1,]
  ) 
  tpm[, 3, 3, i-1] <- 1 
  tpm[, 4, 3, i-1] <- exp(
    gam_cov[,,i-1] %*% b[1,]
  ) 
  # DN to ..
  tpm[, 1, 4, i-1] <- exp(
    rowSums(
      eps_cov[,,i-1] %*% t(d)
    )
  )
  tpm[, 2, 4, i-1] <- exp(
    eps_cov[,,i-1] %*% d[2,]
  ) 
  tpm[, 3, 4, i-1] <- exp( 
    eps_cov[,,i-1] %*% d[1,]
  )
  tpm[, 4, 4, i-1] <- 1
}
  
# Complete softmax for PSI
PSI <- sweep(
  PSI,
  1,
  rowSums(PSI),
  "/"
)



z <- matrix(
  NA,
  ncol = nyear,
  nrow = nsite
)
# Sample what state each site begins
z[,1] <- apply(
  PSI,
  1,
  function(x) extraDistr::rcat(1,x)
)

# Do the rest of the sampling for the latent state
for(i in 2:nyear){
 for(s in 1:nsite){
   z[s,i] <- rcat(
     1,
     tpm[s, , z[s,i-1], i-1] / sum(tpm[s, , z[s,i-1], i-1])
   )
 }
}

# Set up rho detection matrix
rdm <- array(
  NA,
  dim = c(nsite, 4, 4, nyear)
)

for(i in 1:nyear){
  # U
  rdm[,1,1,i] <- 1
  rdm[,2:4,1,i] <- 0
  # D
  rdm[,1,2,i] <- 1
  rdm[,2,2,i] <- exp(
    rho_cov[,,i] %*% f[1,]
  )
  rdm[,3:4,2,i] <- 0
  #N
  rdm[,1,3,i] <- 1
  rdm[,2,3,i] <- 0
  rdm[,3,3,i] <- exp(
    rho_cov[,,i] %*% f[2,]
  )
  rdm[,4,3,i] <- 0
  # DN
  rdm[,1,4,i] <- 1
  rdm[,2,4,i] <- exp(
    rho_cov[,,i] %*% f[1,]
  )
  rdm[,3,4,i] <- exp(
    rho_cov[,,i] %*% f[2,]
  )
  rdm[,4,4,i] <- exp(
    rowSums(
      rho_cov[,,i] %*% t(f)
    )
  )
}
y <- array(
  NA,
  dim = c(nsite, nyear, nsurvey)
)
# simulate observed data
for(i in 1:nyear){
  for(s in 1:nsite){
    y[s,i,] <- rcat(
      nsurvey,
      rdm[s,,z[s,i], i] / sum(rdm[s,,z[s,i], i] )
    )
  }
}

# Set up the data for the model
data_list <- list(
  y = y,
  psi_cov = psi_cov,
  gam_cov = gam_cov,
  eps_cov = eps_cov,
  rho_cov = rho_cov,
  ncat = 2,
  ncov_psi = ncov_psi,
  ncov_rho = ncov_rho,
  ncov_gam = ncov_gam,
  ncov_eps = ncov_eps,
  nyear = nyear,
  nsite = nsite,
  nsurvey = nsurvey
)

library(runjags)

sp_inits <- matrix(
  4,
  nrow = nsite,
  ncol = nyear
)
inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = sp_inits,
      a = matrix(rnorm(data_list$ncat * data_list$ncov_psi),
                 ncol = data_list$ncov_psi, nrow = data_list$ncat),
      b = matrix(rnorm(data_list$ncat * data_list$ncov_gam),
                 ncol = data_list$ncov_gam, nrow = data_list$ncat),
      d = matrix(rnorm(data_list$ncat * data_list$ncov_eps),
                 ncol = data_list$ncov_eps, nrow = data_list$ncat),
      f = matrix(rnorm(data_list$ncat * data_list$ncov_rho),
                 ncol = data_list$ncov_rho, nrow = data_list$ncat),
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

mout_full <- run.jags(model = "./JAGS/jags.dynamic.multistate.null.R",
                      monitor = c("a", "b", "d", "f"), 
                      data = data_list,
                      n.chains = 2,
                      inits = inits,
                      adapt = 400,
                      burnin = 10000,
                      sample = ceiling(50000/2),
                      thin = 5,
                      modules = 'glm',
                      summarise = FALSE,
                      plots = FALSE,
                      method = "parallel")

# save the output
saveRDS(mout_full, "./model_output/model_1_full.rds")
