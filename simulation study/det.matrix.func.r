det.matrix.func=function(pNight,pDay,pND.N,pND.D,pND.notseen){
  det.matrix=matrix(0, ncol=4,nrow=4)
  colnames(det.matrix)=c("not seen","seen day","seen night","seen ND")
  rownames(det.matrix)=c("not used","day use","night use","ND use")

  
  #Define the detection matrix (4 states)
  det.matrix[1,1]=1
  det.matrix[2,1]=1-pDay
  det.matrix[2,2]=pDay
  
  det.matrix[3,1]=1-pNight
  det.matrix[3,3]=pNight
  
  det.matrix[4,1]=pND.notseen
  det.matrix[4,2]=pND.D
  det.matrix[4,3]=pND.N
  det.matrix[4,4]=pND.ND
  det.matrix
}




