#Calculates a Bayesian p-value as a measure of goodness-of-fit

GOF=function(fit,y,model.type){

  cl <- parallel::makeCluster(detectCores())
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
  
  #the likelihood is calculated for each site k
  n.mcmc=length(pDay)
  
  predicted.lik.save=lik.save=matrix(NA, nrow=dim(y)[1],ncol=n.mcmc)
  
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
    
    #Now calcaulte the predicted deviance for site k
    
    #Draw multinomial occurences based on psi matrix.
    #This will be M x n.mcmc
    z.pred=apply(psi.matrix,1,FUN=function(x){rmultinom(1,1,prob=x)})

    #Repackage the detection parameters to be used to create detection matrices
    det.matrix2=cbind(pNight,pDay,pND.ND,pND.N,pND.D,pND.0)

    #Create a list and use lapply to get each s mcmc iteration into a detection matrix
    xy.list <- as.list(as.data.frame(t(det.matrix2)))
    det.out=lapply(xy.list, FUN=function(x){det.matrix.func(x[1],x[2],x[3],x[4],x[5],x[6])})
  

    #Multiple the occurences for each s by detection matrix s to get probabilities to be used for sampling
    linear.pred=list()
    for(s in 1:n.mcmc){
      linear.pred[[s]]=t(z.pred[,s])%*%det.out[[s]]
    }


    #Predict samples for each s by each k site
    y.pred=lapply(linear.pred,FUN=function(x){rmultinom(length(y[k,][!is.na(y[k,])]),1,prob=x)})
    state.occurence=lapply(y.pred,FUN=function(x){apply(x,2,FUN=function(y){which(y==1)})})


 #brute-force for loop 
   for(s in 1:n.mcmc){
     predicted.lik.save[k,s]=multi.state.likelihood(t(matrix(psi.matrix[s,])),t(matrix(det.matrix[s,])),state.occurence[[s]])
   }


  }
  
  Deviance.Observed=-2*apply(log(lik.save),2,sum)
  
  Deviance.Predicted=-2*apply(log(predicted.lik.save),2,sum)
  
  list(Deviance.Observed=Deviance.Observed,Deviance.Predicted=Deviance.Predicted)
  
}  #End function
   
  