################################################
#
# Dynamic simulations for null model
#
# written by M. Fidino
#
################################################

cat("Simulation for a dynamic multi-state null model\n")
cat("-----------------------------------------------\n")

# load utility functinos
source(
  "./simulation study/dynamic_simulations/dynamic_utilities.R"
)

# create the covariates for the model.
my_covars <- sim_covariates()

# to see default values used to generate data:
# args(sim_covariates)

# use those covariates to simulate some data. We are not trying to do an
#  exhaustive simulation study here, so keeping species relatively easy
#  to detect and varying states (and transitions) are relatively common.
my_data <- sim_data(
  my_covars
)
