
# Plotting for Ranomafana & Makira---------------------------------------------- 
#
# Written by Dr., Mason Fidino, Dr. Brian Gerber, and Kimberly Rivera
#
# Rivera, K., et al. (In Review) Rethinking mammal habitat occupancy modeling 
# and the role of diel activity in an anthropogenic world. Preprint Available: 
# https://doi.org/10.1101/2021.06.30.450589

# Clean Up ---------------------------------------------------------------------
rm(list=ls())
gc() #garbage collection

# Ranomafana Plotting-----------------------------------------------------------
## Plotting parameters for fosa-------------------------------------------------

# M1 is the most supported model
# Load the model and plots results

setwd("C:/Users/Kim Rivera/Documents/GitHub/multi.state.temporal.activity/")

# libraries for bayes plotting
library(rjags)
library(runjags)
library(coda)
library(bayesplot)
library(ggplot2)
library(tidyverse)
library("modeest")
library("bayestestR")
library("HDInterval")


# load best fit model and data
load("RNP Fosa/M1.full.fit")
load("RNP Fosa/RNP2.data")
covs=RNP2.data[[2]]
cov=covs$DistTown
cov1.unscaled=as.numeric(cov)
cov1.scaled=as.numeric(scale(cov,center = TRUE,scale = TRUE))
fit.matrix=as.matrix(fit)

# Plot the alpha parameters, e.g. the slopes and intercepts for each state
# Y-axis labels with "Int" indicate an intercept for the states (Day, Night, ND = Night and Day) 
# and "Dist.Vill" indicates the slope parameter associated to the variable distance to village.
png(file="C:/Users/Kim Rivera/Documents/Summer Project/Play with Git/RNP/RNP.fosa.parms.png",res=200,units = "in",height=8,width=8)
color_scheme_set("blue")
mcmc_areas(fit.matrix,
           pars = c("alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]","alpha[5]","alpha[6]"),
           prob = 0.5,
           prob_outer=0.99) + geom_vline(xintercept=0, linetype="solid", 
                                         color = "black", size=1)+
  labs(x = "Logit Value", y="Occupancy State Parameters")+ #working
  scale_y_discrete(labels = c(expression(alpha^"Day,Int"),expression(alpha^"Day,Dist.Vill"), expression(alpha^"Night,Int"), expression(alpha^"Night,Dist.Vill"), expression(alpha^"ND,Int"), expression(alpha^"ND,Dist.Vill"))) +
  theme(text = element_text(family = "URWTimes", size=26), axis.text.y = element_text(family = "serif")) 
dev.off()

# Plot probability of use for use and probability of use for fosa
# for each diel category given fosa are present---------------------------------

library(runjags)
library(markovchain)
library(coda)

#Need to predict site occupancy using new covariates
y <- seq(0, 1, length=1000)

n.mcmc=dim(fit.matrix)[1]


# assigning output to alphas (slopes and intercepts)
alpha1=fit[,which(grepl("alpha[1]",colnames(fit), fixed = TRUE))]
alpha2=fit[,which(grepl("alpha[2]",colnames(fit), fixed = TRUE))]
alpha3=fit[,which(grepl("alpha[3]",colnames(fit), fixed = TRUE))]
alpha4=fit[,which(grepl("alpha[4]",colnames(fit), fixed = TRUE))]
alpha5=fit[,which(grepl("alpha[5]",colnames(fit), fixed = TRUE))]
alpha6=fit[,which(grepl("alpha[6]",colnames(fit), fixed = TRUE))]


x.pred=seq(500,4000,length.out = 100) #distance to nearest village interval
x.pred.scaled=as.numeric(scale(x.pred,center = TRUE,scale = TRUE))
psi.state.mcmc=array(NA, dim=c(length(x.pred),n.mcmc,4))
dim(psi.state.mcmc)

#Get psi estimates for each state by site (for all n.mcmc)
for (s in 1:length(x.pred)){
  phi1 <- 1
  phi2 <- exp(alpha1+alpha2*x.pred.scaled[s])
  phi3 <- exp(alpha3+alpha4*x.pred.scaled[s])
  phi4 <- exp(alpha1+alpha2*x.pred.scaled[s]+
                alpha3+alpha4*x.pred.scaled[s]+
                alpha5+alpha6*x.pred.scaled[s])
  constant=phi1+phi2+phi3+phi4
  psi.state.mcmc[s,,1]     <- phi1/constant
  psi.state.mcmc[s,,2]     <- phi2/constant
  psi.state.mcmc[s,,3]     <- phi3/constant
  psi.state.mcmc[s,,4]     <- phi4/constant
  
}

#example, this is site 1 and state 4
dim(psi.state.mcmc)
hist(psi.state.mcmc[1,,4])

hm2=psi.state.mcmc

day_given_night2 <- apply(psi.state.mcmc[,,c(2,4)], c(1,2), sum)
night_given_day2 <- apply(psi.state.mcmc[,,c(3,4)], c(1,2), sum)
day_not_night2 <- apply(psi.state.mcmc[,,c(1,2)], c(1,2), sum)
night_not_day2 <- apply(psi.state.mcmc[,,c(1,3)], c(1,2), sum)
fosa_present <- apply(psi.state.mcmc[,,2:4], c(1,2), sum)

# wrapper to get quantiles
dq <- function(x,y){
  apply(x, y, quantile, probs = c(0.025,0.5,0.975))
}

# get the quantiles from the steady state object.
quants <- dq(hm2, c(1,3))

# and for plotting, the urbanization metric
xp <- matrix(1, ncol = 2, nrow = 100)
xp[,2] <- seq(500, 4000, length.out = 100)

# state 1.
plot(quants[2,,1] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "unoccupied", xlab= "Distance to Nearest Village (m)",
     ylab = "Probability of occupancy", bty = 'l')
lines(quants[1,,1] ~ xp[,2], lty = 2)
lines(quants[3,,1] ~ xp[,2], lty = 2)

# marginal prob of day (i.e., day + day & night)
day_quants <- dq(day_given_night2,1)

plot(day_quants[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "day use", xlab= "Distance to Nearest Village (m)",
     ylab = "Probability of occupancy", bty = 'l')
lines(day_quants[1,] ~ xp[,2], lty = 2)
lines(day_quants[3,] ~ xp[,2], lty = 2)

# marginal prob of night (i.e., night + day & night)
night_quants <- dq(night_given_day2, 1)
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
fosa_occupancy <- dq(
  apply(hm2[,,2:4], c(1,2), sum),
  1
)

plot(fosa_occupancy[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "Coyote, any state", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(fosa_occupancy[1,] ~ xp[,2], lty = 2)
lines(fosa_occupancy[3,] ~ xp[,2], lty = 2)


#calculate the probability of day use given no night use
cond_dnn <- dq(
  hm2[,,2] / day_not_night2,
  1
)

plot(cond_dnn[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "day use | no night use", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(cond_dnn[1,] ~ xp[,2], lty = 2)
lines(cond_dnn[3,] ~ xp[,2], lty = 2)

# conditional night use given no day use
cond_nnd <- dq(
  hm2[,,3] / night_not_day2,
  1
)

plot(cond_nnd[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "night use | no day use", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(cond_nnd[1,] ~ xp[,2], lty = 2)
lines(cond_nnd[3,] ~ xp[,2], lty = 2)


# The conditional probabilities GIVEN coyote presence 
daynight_cpres <- dq(
  hm2[,,4] / fosa_present,
  1
)

plot(daynight_cpres[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "day night use | coyote present", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(daynight_cpres[1,] ~ xp[,2], lty = 2)
lines(daynight_cpres[3,] ~ xp[,2], lty = 2)

# day use given coyote present
day_cpres  <- dq(
  hm2[,,2] / fosa_present,
  1
)

plot(day_cpres[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "day use | coyote present", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(day_cpres[1,] ~ xp[,2], lty = 2)
lines(day_cpres[3,] ~ xp[,2], lty = 2)


night_cpres <- dq(
  hm2[,,3] / fosa_present,
  1
)

plot(night_cpres[2,] ~ xp[,2], type = 'l', ylim = c(0,1),
     main = "night use | coyote present", xlab= "urban intensity",
     ylab = "Probability of occupancy", bty = 'l')
lines(night_cpres[1,] ~ xp[,2], lty = 2)
lines(night_cpres[3,] ~ xp[,2], lty = 2)

# Plot both marginal and conditional probabilities together---------------------
windows(5.5,11) #creates plotting window
tiff("C:/Users/Kim Rivera/Documents/Summer Project/fosa.rano.occ.4.state.tiff", height = 11, width = 5.5, units = "in",
     res = 600, compression = "lzw")
{
  par(mar = c(5,5,0.5,0.5)+1, mfrow = c(2,1))
  
  
  plot(1~1, ylim = c(0,1), type = "n", xlim = c(500,4000),
       bty = "l", lwd = 2, xlab = "", 
       ylab = "", las = 1, xaxs = "i", yaxs = "i", cex.axis = 1.2) # i allows perfect btw 0,1 & 500&4000
  axis(1, seq(500,4000, 500), labels = FALSE,  tck = -0.015)
  axis(2, seq(0,1, 0.1), labels = FALSE,  tck = -0.015)
  mtext(text = "Distance to Nearest Village (m)", 1 , at = 2250, line = 3.1, cex = 1.5) # at = get center point (500,4000)
  mtext(text = "Probability of Occupancy",2, las = 0, at = 0.5, line = 3.1,
        cex = 1.5)
  
  # get the plotting margins to place sub-figure symbol
  u <- par("usr") # plotting window x left x right, x bottom, y top
  
  text(x = u[1]+300, u[4]-.05, labels = "(a)", cex = 1.5) # check the movement on x-axis (like 50)
  # alpha shading
  alp <- 0.2
  mc <- c("#24d5f7ff", "gray50" ,"#5ee38bff")
  
  
  # do the legends for the 95% CI 
  legend(x = 1500,
         y = .999,
         legend = rep("",4),
         pch = 22,
         pt.bg = scales::alpha(c(mc[c(3,2,1)],"navy"), alp+0.1),
         pt.cex = 3,
         y.intersp = 1.5,
         bty = "n", cex = 1.1)
  legend(x = 1700,
         y = .999,
         legend = c("day use", "night use", "day & night use", "no use"),
         lty = c(3,2,1,4),
         lwd = 3,
         seg.len = 2,
         y.intersp = 1.5,
         bty = "n", cex = 1.1)
  par(xpd = NA)
  text(x = 1400, y = 1.06, labels = "95% CI", pos = 1, cex = 1.1)
  text(y = 1.06, x = 2200, labels = "E(X)", pos = 1, cex = 1.1)
  
  x1 <- xp[,2]
  x2 <- rev(x1)
  
  
  
  y1 <- quants[1,,4]
  y2 <- rev(quants[3,,4])
  
  
  polygon(
    c(x1, x2), c(y1, y2),
    col = scales::alpha(mc[1], alp),
    border = NA
  )
  
  

  #create for no use 
  y1 <- quants[1,,3]
  y2 <- rev(quants[3,,3])
  
  polygon(
    c(x1, x2), c(y1, y2),
    col = scales::alpha(mc[2], alp),
    border = NA
  )
  
  y1 <- quants[1,,1]
  y2 <- rev(quants[3,,1])
  
  polygon(
    c(x1, x2), c(y1, y2),
    col = scales::alpha("navy", alp),
    border = NA
  ) 
  
  # day use 
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
  lines(quants[2,,1] ~ xp[,2], lwd = 3, lty = 4)
  
  #start the conditional plot here
  plot(1~1, ylim = c(0,1), type = "n", xlim = c(500,4000),
       bty = "l", lwd = 2, xlab = "", 
       ylab = "", las = 1, xaxs = "i", yaxs = "i", cex.axis = 1.2)
  axis(1, seq(500,4000, 500), labels = FALSE,  tck = -0.015)
  axis(2, seq(0,1, 0.1), labels = FALSE,  tck = -0.015)
  mtext(text = "Distance to Nearest Village (m)", 1 , at = 2250, line = 3.1, cex = 1.5) # at = get center point (500,4000)
  mtext(text = "Conditional probability of use",2, las = 0, at = 0.5, line = 3.1,
        cex = 1.5)
  text(x = u[1]+300, u[4]-.05, labels = "(b)", cex = 1.5)
  # alpha shadeing
  alp <- 0.2
  mc <- c("#24d5f7ff", "gray50" ,"#5ee38bff")
  
  x1 <- xp[,2]
  x2 <- rev(x1)
  y1 <- daynight_cpres[1,]
  y2 <- rev(daynight_cpres[3,])
  
  
  polygon(
    c(x1, x2), c(y1, y2),
    col = scales::alpha(mc[1], alp),
    border = NA
  )
  
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

# Makira Plotting---------------------------------------------------------------
load("Makira Fosa2/M1.full.no.covs.out")
fit <- combine.mcmc(M1.full.no.covs)
#load the prepared data file
load("Makira Fosa2/Makira.data2")
y=Makira.data2[[2]] #detection history, sites x occs x survey area

#assign covariate data to object
cov=Makira.data2[[3]]

#These are the mean (across survey areas- i.e., mean of random effect) for each logit scaled parameter for states 2, 3, and 4
#alpha 1 applies to state 2 and so on
mu.alpha.1=fit[,which(grepl("mu.alpha1",colnames(fit)))]
mu.alpha.2=fit[,which(grepl("mu.alpha2",colnames(fit)))]
mu.alpha.3=fit[,which(grepl("mu.alpha3",colnames(fit)))]

#Since there are no covariates, these are simply the intercepts on the logit scale
hist(mu.alpha.1)
hist(mu.alpha.2)

#The alpha 3 is how much different state 4 is from simply the combinations of states 2 and 3
#The distribution being a bit negative suggests that we are observing fosa a bit less than expected
#based on simply the combinations of states 2 and 3
hist(mu.alpha.3)

#Derive the mean-level (across all survey areas) state specific occupancy probability
denominator=1+exp(mu.alpha.1)+exp(mu.alpha.2)+exp(mu.alpha.1+mu.alpha.2+mu.alpha.3)
psi.day=exp(mu.alpha.1)/denominator
psi.night=exp(mu.alpha.2)/denominator
psi.ND=exp(mu.alpha.1+mu.alpha.2+mu.alpha.3)/denominator

occ.matrix=cbind(psi.day,
                 psi.night,
                 psi.ND)


# Plot the posterior distributions of fosa state occupancy model parameters-----
# (mean across all sites)
png(file="C:/Users/Kim Rivera/Documents/Summer Project/Play with Git/Makira/Makira.fosa.occ.parms.png",res=200,units = "in",height=8,width=8)
color_scheme_set("blue")
mcmc_areas(occ.matrix, pars = colnames(occ.matrix))+ 
  labs(x = "Probability",y="State Occupancy Parameter")+
  scale_y_discrete("State Occupancy Parameter", labels = c(expression(psi^2),expression(psi^3), expression(psi^4))) +
  theme(text = element_text(family = "URWTimes", size=26), axis.text.y = element_text(family = "serif")) +
  xlim(0,.6) 
dev.off()