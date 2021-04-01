###############################
#
# calculate_steady_state
#
# Written by Mason Fidino
#
###############################

# This function calculates the expected occupancy from the transition matrix
#  in our 3 species coyote, opossum, and raccoon example. As the best fit model
#  indicated that interactions vary as a function of urbanization, we calculate
#  the expected occupancy of each commmunity state across this metric.
#  This function requires:

# Arguments for this function
# --------------------------- #
# mm = model matrix of jags output as a matrix
#  e.g., mm <- as.matrix(jags_output, chains = TRUE) via coda package
#  columns of this matrix represent different parameters while rows are
#  posterior simulations

# data_list = the same data_list supplied to JAGS to fit the DCOM
#   This was generated in fit_softmax_model.R
# ncores = number of cores to run this function on in parallel

# Note that the paramaters are given the same names as in the manuscript.
#  For example, 'a' parameters are for initial occupancy.

# This function returns a list of length MCMC samples (i.e., nrow(mm)). Each
#  element in this list is an 81 x 8 matrix. Rows correspond to the urban
#  PCA that range from -4 to 4 by steps of 0.1 while columns are the 
#  stationary occupancy probability of each state. In our example these states
#  are: No species, coyote, opossum, raccoon, coyote & opossum, coyote & 
#  raccoon, opossum & raccoon, all three species. Each element in the list
#  corresponds to these predictions for that particular step in the MCMC chain
#  from the DCOM.

# Packages required:
#  doParallel
#  markovchain
#  runjags
#  codda

library(runjags)
library(markovchain)
library(coda)

source("./Chicago coyote/plot_utilities.R")
mout <- readRDS("./Chicago coyote/simpler_models/full_model_urb.RDS")

mm <- as.matrix(coda::as.mcmc.list(mout))
mm <- mm[,-grep("^z", colnames(mm))]
rm(mout)

# get quantiles
model_quants <- apply(mm, 2, quantile, probs = c(0.025,0.5,0.975))

# plot out posterior traceplots
for(i in 1:ncol(mm)){
  tmp <- colnames(mm)[i]
  jpeg(paste0("./Chicago coyote/mcmc_plots/",tmp,".jpeg"))
  plot(mm[,i], type = "l", main = tmp)
  dev.off()
}
  
  # Run generate predictions for each step in the MCMC chain
  n <- nrow(mm)
  smat <- vector("list", n) 
  if(!file.exists("./Chicago coyote/simpler_models/steady_state.RDS")){
    pb <- txtProgressBar(1, n)
    for(o in 1:n){
    setTxtProgressBar(pb, o)
    steady_mat <- array(0, dim = c(nrow(xp), model_list$max_state))
    
    tpm <- predict_latent_conditional(
      mm[o,],
      model_list
    )
    # convert to probability, Divide by column sum of each
    for(i in 1:dim(tpm)[1]){
      
      mark <- 
        new("markovchain", states = c("U", "D", "N", "DN"),
            byrow = FALSE, tpm[i,,])
      steady_mat[i,]  <-   steadyStates(mark)
    }
    smat[[o]] <- steady_mat
    }
    saveRDS(smat, "./Chicago coyote/simpler_models/steady_state.RDS")
  } else {
    smat <- readRDS("./Chicago coyote/simpler_models/steady_state.RDS")
  }
  

# make into an array 
hm <- array(unlist(smat), dim = c(250, 4, length(smat)))


# calculate some marginal occupancy probabilities
#  This takes a lil while.
day_given_night <- apply(hm[,c(2,4),], c(1,3), sum)
night_given_day <- apply(hm[,c(3,4),], c(1,3), sum)
day_not_night <- apply(hm[,c(1,2),], c(1,3), sum)
night_not_day <- apply(hm[,c(1,3),], c(1,3), sum)
coyote_present <- apply(hm[,2:4,], c(1,3), sum)

# wrapper to get quantiles
dq <- function(x,y){
 apply(x, y, quantile, probs = c(0.025,0.5,0.975))
}

# get the quantiles from the steady state object.
quants <- dq(hm, c(1,2))

# and for plotting, the urbanization metric
xp <- matrix(1, ncol = 2, nrow = 250)
xp[,2] <- seq(-4, 4, length.out = 250)

# state 1.
plot(quants[2,,1] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "unoccupied", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(quants[1,,1] ~ xp[,2], lty = 2)
lines(quants[3,,1] ~ xp[,2], lty = 2)

# marginal prob of day (i.e., day + day & night)
day_quants <- dq(day_given_night,1)
plot(day_quants[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "day use", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(day_quants[1,] ~ xp[,2], lty = 2)
lines(day_quants[3,] ~ xp[,2], lty = 2)

# marginal prob of night (i.e., night + day & night)
night_quants <- dq(night_given_day, 1)
plot(night_quants[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "night use", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(night_quants[1,] ~ xp[,2], lty = 2)
lines(night_quants[3,] ~ xp[,2], lty = 2)

# day & night use (state 4)
plot(quants[2,,4] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "day night use", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(quants[1,,4] ~ xp[,2], lty = 2)
lines(quants[3,,4] ~ xp[,2], lty = 2)

# overall occupancy
coyote_occupancy <- dq(
  apply(hm[,2:4,], c(1,3), sum),
  1
)

plot(coyote_occupancy[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "Coyote, any state", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(coyote_occupancy[1,] ~ xp[,2], lty = 2)
lines(coyote_occupancy[3,] ~ xp[,2], lty = 2)


#calculate the probability of day use given no night use
cond_dnn <- dq(
  hm[,2,] / day_not_night,
  1
)

plot(cond_dnn[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "day use | no night use", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(cond_dnn[1,] ~ xp[,2], lty = 2)
lines(cond_dnn[3,] ~ xp[,2], lty = 2)

# conditional night use given no day use
cond_nnd <- dq(
  hm[,3,] / night_not_day,
  1
)

plot(cond_nnd[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "night use | no day use", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(cond_nnd[1,] ~ xp[,2], lty = 2)
lines(cond_nnd[3,] ~ xp[,2], lty = 2)


# The conditional probabilities GIVEN coyote presence 
daynight_cpres <- dq(
  hm[,4,] / coyote_present,
  1
)

plot(daynight_cpres[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "day night use | coyote present", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(daynight_cpres[1,] ~ xp[,2], lty = 2)
lines(daynight_cpres[3,] ~ xp[,2], lty = 2)

# day use given coyote present
day_cpres  <- dq(
  hm[,2,] / coyote_present,
  1
)

plot(day_cpres[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "day use | coyote present", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(day_cpres[1,] ~ xp[,2], lty = 2)
lines(day_cpres[3,] ~ xp[,2], lty = 2)


night_cpres <- dq(
  hm[,3,] / coyote_present,
  1
)

plot(night_cpres[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "night use | coyote present", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(night_cpres[1,] ~ xp[,2], lty = 2)
lines(night_cpres[3,] ~ xp[,2], lty = 2)

windows(5,5)
windows(5.5,11)

# plot all together?
tiff("./Chicago coyote/figures/diel_prob_of_use.tiff", height = 11, width = 5.5, units = "in",
     res = 600, compression = "lzw")
{
  par(mar = c(5,5,0.5,0.5), mfrow = c(2,1))

  plot(1~1, ylim = c(0,1), type = "n", xlim = c(-4,4),
       bty = "l", lwd = 2, xlab = "", 
       ylab = "", las = 1, xaxs = "i", yaxs = "i", cex.axis = 1.2)
  axis(1, seq(-4,4, 0.5), labels = FALSE,  tck = -0.015)
  axis(2, seq(0,1, 0.1), labels = FALSE,  tck = -0.015)
  mtext(text = "Urban intensity", 1 , at = 0, line = 3.1, cex = 1.5)
  mtext(text = "Probability of Occupancy",2, las = 0, at = 0.5, line = 3.1,
        cex = 1.5)
  # get the plotting margins to place sub-figure symbol
  u <- par("usr")
  text(x = u[1]+0.5, u[4]-0.05, labels = "(a)", cex = 1.5)
  # alpha shadeing
  alp <- 0.2
  mc <- c("#24d5f7ff", "gray50" ,"#5ee38bff")#"#ffb226ff" )
  
  
  # do the legends for the 95% CI 
  legend(x = -1,
         y = 0.95,
         legend = rep("",3),
         pch = 22,
         pt.bg = scales::alpha(mc[c(3,2,1)], alp+0.1),
         pt.cex = 3,
         y.intersp = 1.5,
         bty = "n", cex = 1.1)
  legend(x = -0.35,
         y = 0.95,
         legend = c("day use", "night use", "night & day use"),
         lty = c(3,2,1),
         lwd = 3,
         seg.len = 2,
         y.intersp = 1.5,
         bty = "n", cex = 1.1)
  par(xpd = NA)
  text(x = -0.85, y = 1, labels = "95% CI", pos = 1, cex = 1.1)
  text(y = 1, x = 0.17, labels = "E(X)", pos = 1, cex = 1.1)
  
  x1 <- xp[,2]
  x2 <- rev(x1)
  y1 <- quants[1,,4]
  y2 <- rev(quants[3,,4])
  
  
  polygon(
    c(x1, x2), c(y1, y2),
    col = scales::alpha(mc[1], alp),
    border = NA
  )
  
  
  
  
  # polygon(
  #   c(x1, x2), c(y1, y2),
  #   col = scales::alpha(mc[1], alp),
  #   border = NA
  # )
  
  y1 <- quants[1,,3]
  y2 <- rev(quants[3,,3])
  
  polygon(
    c(x1, x2), c(y1, y2),
    col = scales::alpha(mc[2], alp),
    border = NA
  )
  
  
  
  y1 <- quants[1,,2]
  y2 <- rev(quants[3,,2])
  
  polygon(
    c(x1, x2), c(y1, y2),
    col = scales::alpha(mc[3], alp),
    border = NA
  )
  lines(quants[2,,4] ~ xp[,2], lwd = 3)
  lines(quants[2,,3] ~ xp[,2], lwd = 3, lty = 2)
  lines(quants[2,,2] ~ xp[,2], lwd = 3, lty = 3)
  
plot(1~1, ylim = c(0,1), type = "n", xlim = c(-4,4),
     bty = "l", lwd = 2, xlab = "", 
     ylab = "", las = 1, xaxs = "i", yaxs = "i", cex.axis = 1.2)
axis(1, seq(-4,4, 0.5), labels = FALSE,  tck = -0.015)
axis(2, seq(0,1, 0.1), labels = FALSE,  tck = -0.015)
mtext(text = "Urban intensity", 1 , at = 0, line = 3.1, cex = 1.5)
mtext(text = "Conditional probability of use",2, las = 0, at = 0.5, line = 3.1,
      cex = 1.5)
text(x = u[1]+0.5, u[4]-0.05, labels = "(b)", cex = 1.5)
# alpha shadeing
alp <- 0.2
mc <- c("#24d5f7ff", "gray50" ,"#5ee38bff")#"#ffb226ff" )


# do the legends for the 95% CI 
# legend(x = -1,
#        y = 0.95,
#        legend = rep("",3),
#        pch = 22,
#        pt.bg = scales::alpha(mc[c(3,2,1)], alp+0.1),
#        pt.cex = 3,
#        y.intersp = 1.5,
#        bty = "n", cex = 0.9)
# legend(x = -0.35,
#        y = 0.95,
#        legend = c("day use", "night use", "night & day use"),
#        lty = c(3,2,1),
#        lwd = 3,
#        seg.len = 2,
#        y.intersp = 1.5,
#        bty = "n", cex = 0.9)
# par(xpd = NA)
# text(x = -0.85, y = 1, labels = "95% CI", pos = 1, cex = 0.9)
# text(y = 1, x = 0.17, labels = "E(X)", pos = 1, cex = 0.9)

x1 <- xp[,2]
x2 <- rev(x1)
y1 <- daynight_cpres[1,]
y2 <- rev(daynight_cpres[3,])


polygon(
  c(x1, x2), c(y1, y2),
  col = scales::alpha(mc[1], alp),
  border = NA
)




# polygon(
#   c(x1, x2), c(y1, y2),
#   col = scales::alpha(mc[1], alp),
#   border = NA
# )

y1 <- night_cpres[1,]
y2 <- rev(night_cpres[3,])

polygon(
  c(x1, x2), c(y1, y2),
  col = scales::alpha(mc[2], alp),
  border = NA
)



y1 <- day_cpres[1,]
y2 <- rev(day_cpres[3,])

polygon(
  c(x1, x2), c(y1, y2),
  col = scales::alpha(mc[3], alp),
  border = NA
)
lines(daynight_cpres[2,] ~ xp[,2], lwd = 3)
lines(night_cpres[2,] ~ xp[,2], lwd = 3, lty = 2)
lines(day_cpres[2,] ~ xp[,2], lwd = 3, lty = 3)
}
dev.off()

# figure out max state
tmp_det <- apply(y, c(1,2), paste, collapse = "-")
highest_obs <- matrix(NA, ncol = 13, nrow = 105)

has2 <- grep("2", tmp_det)
has3 <- grep("3", tmp_det)
has4 <- grep("4", tmp_det)
# only 2 observed
highest_obs[has2[which(!has2 %in% c(has3, has4))]] <- 2
# only 3 observed
highest_obs[has3[which(!has3 %in% c(has2, has4))]] <- 3
# 4 observed
highest_obs[has4] <- 4
highest_obs[has2[which(has2 %in% has3)]] <- 4

# get our urbanization gradients 
curb <- data_list$gam_cov[,2,1]
ccut <- cut(curb, breaks = seq(-3.5, 3.5, length.out = 7))

cdt <- read.csv("../aly_data_request/datetime/datetime_data_chicago.csv")

# get only coyote
cdt <- cdt[which(cdt$Species == "coyote"),]

# get only seasons in our study
cdt <- cdt[cdt$Season %in% unique(coyote$Season),]

# get only sites in our study
cdt <- cdt[cdt$Site %in% unique(coyote$Site),]

# convert to 24 hour time

cdt$dt <- lubridate::ymd_hms(cdt$Datetime)
cdt$time <- format(cdt$dt, format = "%H:%M:%S")
# convert to 24 hour time
ack <- as.numeric(lubridate::hms(cdt$time)) / 86400

cdt$rads <- ack * 2 *pi



# figure out where those lines split
xp[min(which((daynight_cpres[2,] - night_cpres[2,]) < 0)),2]

above1.56 <- covs$Site[covs$urb>0]

cdt <- cdt[-which(duplicated(cdt)),]
library(overlap)
overlapPlot(
  cdt$rads[cdt$Site %in% above1.56],
  cdt$rads[!cdt$Site %in% above1.56],
  xaxs ="i",
  yaxs = "i",
  bty = "l",
  main= ""
)

legend('top',
       c("Activity > 1.56 urb", 
         "Activity < 1.56 urb"),
       lty=c(1,2),
       col=c(1,4), bty='n')
ack <- (lubridate::hour(cdt$dt)) + 
  (lubridate::minute(cdt$dt)/60) +
  (lubridate::second(cdt$dt)/(60*2))
ack <- ack



library(overlap)

ob_state <- matrix(0, ncol = 3, nrow = 6)
for(i in 1:6){
  tmp_dat <- highest_obs[which(ccut == levels(ccut)[i]),]
  ob_state[i,] <- prop.table(tabulate(tmp_dat,4)[2:4])
}

# get average prob of each state across these values


# break it into a few pieces
round(quants,2)
# get naive occupancy
nest <- matrix(NA, ncol = 4, nrow = model_list$nseason)

for(i in 1:nrow(nest)){
  tmp <- model_list$y[,i,]
  tmp <- apply(tmp, 1, max, na.rm = TRUE)
  tmp <- tmp[!is.infinite(tmp)]
  nest[i,] <- prop.table(tabulate(tmp,4))
}

apply(nest, 2, mean)


aa <- c(0, -1.39, -0.32, -1.39-0.32)
aa <- exp(aa)
aa <- aa/ sum(aa)

round(aa, 2)

round(quants,2)

tpm <- predict_latent_conditional(
  apply(mm,2, median),
  model_list
)

tpm <- round(tpm[1,,], 2)

colnames(tpm) <- c("U", "D", "N", "DN")
row.names(tpm) <- colnames(tpm)
