##############################################
#
# Getting model objects ready for analysis
#
# Written by M. Fidino
#
##############################################

# I made this script to abstract out some of the code used in:

#  1. fit_models.R
#  2. model_selection.R

# Both of them had the same code at the beginning, so it made sense to 
#  have all of that code in one spot.

# I've already pre-scrubbed the detection data in:
#  file.edit("./Chicago coyote/scrub_coyote_data.R")

coyote <- read.csv(
  "./Chicago coyote/data/day_night_detections.csv"
)

# bring in urbanization covariates
covs <- read.csv(
  "./Chicago coyote/data/chicago_covars.csv"
) 
colnames(covs)[1] <- "Site"
# join to coyote

coyote <- dplyr::right_join(
  covs[,c("Site", "urb")],
  coyote,
  by = "Site"
)

# The data needs to placed in a site by season by sample array. Given
#  it's current arragement it'll be easier to place it in a sample
#  by site by season array and then we can rotate it.

nsurvey <- length(grep("^Week", colnames(coyote)))
nseason <- length(unique(coyote$Season))
nsite <- length(unique(coyote$Site))

y <- array(
  t(coyote[,grep("^Week", colnames(coyote))]),
  dim = c(nsurvey, nsite, nseason)
)

y <- aperm(
  y,
  c(2,3,1)
)


# These are the models we are going to fit
models <- data.frame(
  formula = c(
    rep(
      c(
        "y ~  season * urb",
        "y ~  season + urb",
        "y ~ season",
        "y ~ urb",
        "y ~ 1"
      ),
      4
    )
  ),
  indep_rho = c(
    rep(TRUE,5),
    rep(FALSE, 15)
    
  ),
  inxs = c(
    rep(TRUE, 10),
    rep(FALSE,10)
  ),
  true_null = c(
    rep(FALSE, 15),
    rep(TRUE, 5)
  )
)
