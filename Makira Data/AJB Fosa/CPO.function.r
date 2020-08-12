#Script to calculate CPO from a model

CPO.function=function(fit,y,model.type){
  
  #To do model selection with CPO, we need to get site-level occurence and detection
  #probabilities.
  
  pDay=fit[,which(grepl("pDay",colnames(fit)))]
  pNight=fit[,which(grepl("pNight",colnames(fit)))]
  
  if(model.type=="null"){
    pDay=fit[,which(grepl("p.overall",colnames(fit)))]
    pNight=pDay
    
    mat=cbind(log(3),qlogis(pDay),qlogis(pDay),qlogis(pDay))
    mat=exp(mat)
    omega.det=t(apply(mat,1,FUN=function(x){x/sum(x)}))
    
    pND.0=omega.det[,1]
    pND.D=omega.det[,2]
    pND.N=omega.det[,3]
    pND.ND=omega.det[,4]
  }
  
  if(model.type=="full"){
  pND.0=fit[,which(grepl("pND",colnames(fit)))[1]]
  pND.D=fit[,which(grepl("pND",colnames(fit)))[2]]
  pND.N=fit[,which(grepl("pND",colnames(fit)))[3]]
  pND.ND=fit[,which(grepl("pND",colnames(fit)))[4]]
  }
  if(model.type=="reduced"){
    pND.0=(1-pNight)*(1-pDay)
    pND.D=pDay*(1-pNight)
    pND.N=pNight*(1-pDay)
    pND.ND=pDay*pNight 
  }
  
  det.matrix=cbind(pDay,
                   pNight,
                   pND.0,
                   pND.D,
                   pND.N,
                   pND.ND)
  
  #the likelihood is calcaulted for each site k
  n.mcmc=length(pDay)
  
  lik.save=matrix(NA, nrow=dim(y)[1],ncol=n.mcmc)
  
  #loop through each site k
  for(k in 1:dim(lik.save)[1]){
    
    site.num=k
    index.occ=which(apply(as.matrix(colnames(fit)), 1, FUN=function(x){grepl(paste("PSI[",site.num,",",sep=""), x, fixed = TRUE)}))
    
    #For site q in order of states by 1,2,3,4
    
    psi.matrix=fit[,index.occ]
    
    #calcualte site level likelihood
    if(all(is.na(y[k,]))==FALSE)(
    lik.save[k,]=multi.state.likelihood(psi.matrix,det.matrix,y[k,])
    )
  }
  
  
  #CPO calculation from site-level liklihood
  CPO.sites=n.mcmc/(apply(1/lik.save,1,sum,na.rm=TRUE))
  CPO.sites <- CPO.sites[is.finite(CPO.sites)]
  length(CPO.sites)
  CPO=(-1)*sum(log(CPO.sites),na.rm = TRUE)
  CPO
  
}