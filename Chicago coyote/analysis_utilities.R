

make_model_matrix <- function(formula, df, rho = FALSE, comp_winter = FALSE){
  
  nsite <- length(unique(df$Site))
  nseason <- length(unique(df$Season))
  form <- as.formula(formula)
  tmp <- rep(1, nrow(df))
  
  if(comp_winter){
    seas <- rep(0, nrow(df))
    seas[grep("JA", df$Season)] <- 1
  } else {
    seas <- factor(
      substr( df$Season, 1, 2 ),
      levels = c("JA", "AP", "JU", "OC")
    )
  }
  if(length(grep("urb", formula)) > 0){
    mod_df <- data.frame(
      y = rep(1, nrow(df)),
      season = seas,
      urb = df$urb
    )
  } else {
    mod_df <- data.frame(
      y = tmp,
      season = seas
    )
  }
  if(!rho){
    mod_df <- mod_df[-c(1:nsite),]
  }
  mmat <- model.matrix(
    as.formula(form),
    mod_df
  )
  cnames <- colnames(mmat)
  mmat <- t(mmat)
  dim(mmat) <- c(nrow(mmat), nsite, ifelse(rho, nseason, nseason - 1))
  
  mmat <- aperm(mmat, c(2,1,3))
  dimnames(mmat)[[2]] <- cnames
  return(mmat)
  
}
