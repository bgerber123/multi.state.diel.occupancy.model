

make_model_matrix <- function(formula, df, rho = FALSE){
  
  nsite <- length(unique(df$Site))
  nseason <- length(unique(df$Season))
  form <- as.formula(formula)
  tmp <- rep(1, nrow(df))
  

  if(length(grep("urb", formula)) > 0){
    mod_df <- data.frame(
      y = rep(1, nrow(df)),
      urb = df$urb
    )
  } else {
    mod_df <- data.frame(
      y = tmp
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
