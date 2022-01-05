##############################
#
# Making model predictions
#
# written by: M. Fidino
#
##############################

# Here we are evaulating the best fit-model to make some figures.

library(runjags)
library(coda)
library(png)
library(grid)


source("./Chicago coyote/plot_utilities.R")
mout <- readRDS("./Chicago coyote/simpler_models/full_model_urb.RDS")


yo <-summary(mout, vars = c("a", "a_inxs", "b", "d", "f", "g", "h"))

# get the full posterior
mcmc <- as.matrix(as.mcmc.list(mout))

# drop z's
mcmc <- mcmc[,-grep("z", colnames(mcmc))]
if(!file.exists("./Chicago coyote/simpler_models/tpm_preds.RDS")){

my_output <- array(NA, dim = c(40000,250,4,4))
pb <- txtProgressBar(0, 40000)

for(o in 1:40000){
setTxtProgressBar(pb, o)
  my_output[o,,,] <- predict_latent_conditional(
    mcmc = mcmc[o,],
    model_list
  )
}

latents <- apply(
  my_output,
  c(2,3,4),
  quantile,
  probs = c(0.025,0.5,0.975)
)
saveRDS(latents, "./Chicago coyote/simpler_models/tpm_preds.RDS")
} else {
  latents <- readRDS("./Chicago coyote/simpler_models/tpm_preds.RDS") 
}


#windows(8,8)

tiff("./Chicago coyote/figures/transition_plot.tiff",
     height = 8, width = 8, units = "in", res = 900, 
     compression = "lzw")
nr <- 2
hm <- matrix(1:16, ncol = 4, nrow = 4, byrow = TRUE)
m <- matrix(
  c(rep(0,( nr * 4) + 2),
  0, rep(hm[1:4], each = nr),0,
  0, rep(hm[1:4], each = nr),0,
  0, rep(hm[5:8], each = nr),0,
  0, rep(hm[5:8], each = nr),0,
  0, rep(hm[9:12], each = nr),0,
  0, rep(hm[9:12], each = nr),0,
  0, rep(hm[13:16], each = nr),0,
  0, rep(hm[13:16], each = nr),0,
  rep(0,( nr * 4) + 2)),
  ncol = (nr * 4) +2 ,
  nrow = 10,
  byrow = TRUE
)
layout(m)

par(mar = c(1.25,1,0.5,1.25), xpd = NA)


urb <- seq(-4,4, length.out = 250)
locs <- expand.grid(1:4, 1:4)
add_png <- TRUE
for(i in 1:16){
  plot(1~1, ylim = c(0,1), xlim = c(-4,4),
       xlab = "", ylab = "",
       bty = "l", type = "l", xaxs = "i", yaxs = "i",
       xaxt = "n", yaxt = "n")
  axis(1, seq(-4,4, length.out = 5), labels = FALSE,  tck = -0.05)
  axis(1, seq(-4,4, 0.5), labels = FALSE,  tck = -0.05/2)

  if(i %in% c(1:4)){
    mtext(
      sprintf("%0.2f", seq(0, 1, 0.25)), 2, at = seq(0,1,0.25), 
      cex = 1, line = 1 ,las = 1)
  }
  if(i %in% c(4,8,12,16)){
  mtext(seq(-4,4, length.out = 5), 1, at = seq(-4,4, length.out = 5),
        line = 1)
  }
  axis(2, seq(0,1,0.25), labels = FALSE,  tck = -0.05)
  axis(2, seq(0,1,0.25/2), labels = FALSE,  tck = -0.05/2)
  
  
  x1 <- urb
  x2 <- rev(x1)
  y1 <- latents[1,,locs[i,1], locs[i,2]]
  y2 <- rev(latents[3,,locs[i,1], locs[i,2]])
  mc <- c("#24d5f7ff", "gray50" ,"#5ee38bff")
  
  polygon(
    c(x1, x2), c(y1, y2),
    col = scales::alpha(mc[1], 0.25),
    border = NA
  )

  lines(c(-4,-4), y = c(0,1))
  lines(latents[2,,locs[i,1], locs[i,2]] ~ urb, lwd = 2, col = mc[1])
  lines(c(-4,4), y = c(0,0))
  if(i == 1){
    par(xpd = NA)
    ty <- 1.5
    tc <- 1.6
    text(0, y = ty, "'no use' to...", cex = tc)
    mtext("Probability of transition", 2, line = 4.75, at = -1.25 ,cex = 1.6)
    if(add_png){
    dpng <- readPNG("./Chicago coyote/pngs/no_coyote.png")
    grid.raster(dpng,x = 0.195, y = 0.92, width = 0.064)
    }
  }
  if(i == 5){
    par(xpd = NA)
    text(0, y = ty, "'day use' to...", cex = tc)
    if(add_png){
    dpng <- readPNG("./Chicago coyote/pngs/day_coyote.png")
    grid.raster(dpng,x = 0.395, y = 0.92, width = 0.08)
    }
  }
  if(i == 9){
    par(xpd = NA)
    text(0, y = ty, "'night use' to...", cex = tc)
    if(add_png){
    dpng <- readPNG("./Chicago coyote/pngs/night_coyote.png")
    grid.raster(dpng,x = 0.595, y = 0.92, width = 0.064)
    }
  }
  if(i == 13){
    par(xpd = NA)
    text(0, y = ty, "'night & day use' to...", cex = tc)
    text(8.85, 0.5, "'no use'", srt = 270, cex = tc)
    if(add_png){
    dpng <- readPNG("./Chicago coyote/pngs/day_night_coyote.png")
    grid.raster(dpng,x = 0.795, y = 0.92, width = 0.08)
    
    dpng <- readPNG("./Chicago coyote/pngs/no_coyote.png")
    grid.raster(dpng,x = 0.925, y = 0.81, width = 0.064)
    }
  }
  if(i == 14){
    text(8.85, 0.5, "'day use'", srt = 270, cex = tc)
    if(add_png){
    dpng <- readPNG("./Chicago coyote/pngs/day_coyote.png")
    grid.raster(dpng,x = 0.925, y = 0.61, width = 0.08)
    }
  }
  if(i == 15){
    text(8.85, 0.5, "'night use'", srt = 270, cex = tc)
    if(add_png){
    dpng <- readPNG("./Chicago coyote/pngs/night_coyote.png")
    grid.raster(dpng,x = 0.925, y = 0.41, width = 0.064)
    }
  }
  if(i == 16){
    text(8.85, 0.5, "'night and day use'", srt = 270, cex = tc)
    mtext("Urban intensity", 1, line = 4.5, at = -14.75 ,cex = 1.6)
    if(add_png){
    dpng <- readPNG("./Chicago coyote/pngs/day_night_coyote.png")
    grid.raster(dpng,x = 0.925, y = 0.21, width = 0.08)
    }
  }
  
  
}

dev.off()
