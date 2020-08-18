#Script to calculate CPO from a model

CPO.function.RE=function(fit,y){ #model.type
  
  #To do model selection with CPO, we need to get site-level occurence and detection
  #probabilities.
  
  #these are survey-level state-occurence probabilities
  PSI.index=which(apply(as.matrix(colnames(fit)), 1,FUN=function(x){grepl(paste("PSI[",sep=""),x, fixed = TRUE)}))
  
  #the order goes site, state, survey
  fit.PSI=fit[,PSI.index]
  
  #colnames(fit.PSI)
  n.mcmc=dim(fit.PSI)[1]

  #the number of surveys
  Ns=dim(y)[3]
  
  #create an object to save the likelihood that is site by state by survey
  lik.save=array(NA, dim=c(dim(y)[1],n.mcmc,Ns))

#Loop over the number of surveys  
  for(u in 1:Ns){
    
    
    #Need to get survey-level state detection. Assuming there is no variation in detection by site
    pday.index=which(apply(as.matrix(colnames(fit)), 1,FUN=function(x){grepl(paste("pDay[",u,"]",sep=""), x, fixed = TRUE)}))
    pnight.index=which(apply(as.matrix(colnames(fit)), 1,FUN=function(x){grepl(paste("pNight[",u,"]",sep=""),x, fixed = TRUE)}))
    pND.ND.index=which(apply(as.matrix(colnames(fit)), 1,FUN=function(x){grepl(paste("pND.ND[",u,"]",sep=""),x, fixed = TRUE)}))
    pND.0.index=which(apply(as.matrix(colnames(fit)), 1,FUN=function(x){grepl(paste("pND.0[",u,"]",sep=""),x, fixed = TRUE)}))
    pND.D.index=which(apply(as.matrix(colnames(fit)), 1,FUN=function(x){grepl(paste("pND.D[",u,"]",sep=""),x, fixed = TRUE)}))
    pND.N.index=which(apply(as.matrix(colnames(fit)), 1,FUN=function(x){grepl(paste("pND.N[",u,"]",sep=""),x, fixed = TRUE)}))
    
    #package detection paramters to be used in the multi.state.likelihood function
    det.matrix=cbind(fit[,pday.index],
                     fit[,pnight.index],
                     fit[,pND.0.index],
                     fit[,pND.D.index],
                     fit[,pND.N.index],
                     fit[,pND.ND.index])
    
    #grab state occurenecs for surevy u
    survey.index=which(apply(as.matrix(colnames(fit.PSI)), 1,FUN=function(x){grepl(paste(",",u,"]",sep=""),x, fixed = TRUE)}))
    fit.PSI.survey=fit.PSI[,survey.index]
    
   
    #loop through each site k
    for(k in 1:dim(lik.save)[1]){
      
      site.num=k
      index.occ=which(apply(as.matrix(colnames(fit.PSI.survey)), 1, FUN=function(x){grepl(paste("PSI[",site.num,",",sep=""), x, fixed = TRUE)}))
      
      #For site k and survey u in order of states by 1,2,3,4
      psi.matrix=fit.PSI.survey[,index.occ]
      
      #calcualte site level likelihood
      if(all(is.na(y[k,,u]))==FALSE)(
        lik.save[k,,u]=multi.state.likelihood(psi.matrix,det.matrix,y[k,,u])
      )
      
    }
   # print(u)
  }
  
  #Need combine across all surveys
  lik.save.sites=NULL
  for(u in 1:Ns){
    lik.save.sites=c(lik.save.sites,n.mcmc/apply(1/lik.save[,,u],1,sum,na.rm=TRUE))
  }
  
  
  CPO.sites <- lik.save.sites[is.finite(lik.save.sites)]
  length(CPO.sites)
  CPO=(-1)*sum(log(CPO.sites),na.rm = TRUE)
  CPO
  
  
}