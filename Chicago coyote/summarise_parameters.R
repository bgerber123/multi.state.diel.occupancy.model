#####################################
#
# Summarise model parameters
#
# written by: M. Fidino
#
#####################################


library(runjags)
library(coda)
library(png)
library(grid)


source("./Chicago coyote/plot_utilities.R")
mout <- readRDS("./Chicago coyote/simpler_models/full_model_urb.RDS")


yo <-summary(mout, vars = c("a", "a_inxs", "b", "d", "f", "g", "h"))

# get the full posterior
mcmc <- as.matrix(as.mcmc.list(mout))
mcmc <- mcmc[,-grep("z", colnames(mcmc))]

# generate ex. occupancy across different states

ex_mat <- matrix(NA, ncol = 4, nrow = nrow(mcmc))

ex_mat[,1] <- 1
ex_mat[,2] <- exp(mcmc[,"a[1,1]"])
ex_mat[,3] <- exp(mcmc[,"a[2,1]"])
ex_mat[,4] <- exp(rowSums(mcmc[,c("a[1,1]", "a[2,1]", "a_inxs[1]")]))

ex_mat <- ex_mat / rowSums(ex_mat)

# get expected average occupancy for the first time step
round(
  apply(
    ex_mat,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
    ),
  2
)


# try parameters for day state
ex_mat <- matrix(NA, ncol = 4, nrow = nrow(mcmc))

ex_mat[,1] <- exp(mcmc[,"d[1,1]"])
ex_mat[,2] <- 1
ex_mat[,3] <- exp(rowSums(mcmc[,c("d[1,1]", "b[2,1]", "g[1,1]")]))
ex_mat[,4] <- exp(rowSums(mcmc[,c("b[2,1]", "g[1,1]")]))

ex_mat <- ex_mat / rowSums(ex_mat)

# get expected average occupancy for the first time step
round(
  apply(
    ex_mat,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  ),
  2
)


# and the night state
ex_mat <- matrix(NA, ncol = 4, nrow = nrow(mcmc))

ex_mat[,1] <- exp(mcmc[,"d[2,1]"])
ex_mat[,2] <- exp(rowSums(mcmc[,c("d[2,1]", "b[1,1]", "g[2,1]")]))
ex_mat[,3] <- 1
ex_mat[,4] <- exp(rowSums(mcmc[,c("b[1,1]", "g[2,1]")]))

ex_mat <- ex_mat / rowSums(ex_mat)

# get expected average occupancy for the first time step
round(
  apply(
    ex_mat,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  ),
  2
)

head(round(yo, 2), 26)[,1:4]
# And day night detection state
round(plogis(head(yo, 26)[,1:4]),2)

ex_mat <- matrix(NA, ncol = 4, nrow = nrow(mcmc))

ex_mat[,1] <- 1
ex_mat[,2] <- exp(mcmc[,"f[1,1]"])
ex_mat[,3] <- exp(mcmc[,"f[2,1]"])
ex_mat[,4] <- exp(rowSums(mcmc[,c("f[1,1]", "f[2,1]")]))

ex_mat <- ex_mat / rowSums(ex_mat)

round(
  apply(
    ex_mat,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  ),
  2
)

round(
  quantile(
    plogis(
      mcmc[,c("f[2,1]", "f[2,2]")] %*% c(1,4.8)
      ), probs = c(0.025,0.5,0.975)),
  2
)

      